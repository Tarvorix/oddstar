"""Narrowband search of Breakthrough Listen fine-resolution scans held on local disk.

The radio channel measures one small time x frequency block per scan. A fine-resolution scan is one block of 16
spectra by several hundred million channels, so this module does the part the channel cannot: it walks the ON
scan in blocks of `live_block_channels`, finds every `live_window_channels`-wide window whose drift-summed power
stands out (`prefilter_block`, the channel's own statistic and drift grid on a per-window noise estimate), cuts the
same frequencies from the OFF scans, and hands each such window to the channel as an ordinary ON/OFF cadence.
Every number that decides a flag comes from `RadioNarrowbandChannel.measure` and `.score`; nothing here re-decides.

What a star's record holds after `enrich_sample_radio`:

- `radio_fetch.status`: "ok" (a complete scan set was on disk and searched), "not_downloaded" (the star is inside a
  beam and a complete ON/OFF set exists in the archive, but its files are not in the data directory),
  "no_comparison" (inside a beam, but the archive has no OFF scans for it), "no_coverage" (inside no half-power
  beam), "failed: ..." (the index or a file could not be read), "not_fetched" (beyond `max_stars`).
- For "ok": `scans`, `dt_s`, `df_hz`, `center_freq_mhz` of one window (the one hardest to explain: flagged outside
  the allocated bands, else flagged, else the strongest ON-only window, else the strongest window), `in_rfi_band`
  for it, and `radio_search` with the counts and the list of evaluated windows.

Limits that are stated, not hidden: only frequencies where the star is inside the half-power beam of the pointing
are searched; the drift grid is the channel's (161 trial rates, 0.05 Hz/s apart at the default range), so a signal
between two trial rates is smeared over a few channels; one window yields one hit per scan (its strongest), so a
weaker signal next to a stronger one is not seen; a set with a single ON scan cannot test persistence, so
interference that happens to be on during the ON scan only is flagged and has to be removed by vetting or by a
second epoch. A flag is a lead, not a detection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import struct
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io import bl_opendata as bl
from oddstar.io.cache import RawCache

SEARCH_VERSION = "bl_local_search_v1"
# Trial drift rates per search; must equal the grid inside `channels.radio_narrowband.dedoppler_search_scan`
N_DRIFT_TRIALS = 161
# Bands allocated to satellite downlinks, radionavigation satellites and aeronautical radionavigation (ITU Radio
# Regulations Article 5, Region 2), where a narrowband line toward any sky position is expected from transmitters.
# (low MHz, high MHz, service)
RFI_ALLOCATIONS_MHZ: tuple[tuple[float, float, str], ...] = (
    (960.0, 1164.0, "aeronautical radionavigation (DME, TACAN, ADS-B)"),
    (1164.0, 1300.0, "radionavigation satellites (GPS L5/L2, GLONASS L3/L2, Galileo E5/E6)"),
    (1525.0, 1559.0, "mobile satellite downlinks (Inmarsat)"),
    (1559.0, 1610.0, "radionavigation satellites (GPS L1, GLONASS L1, Galileo E1, BeiDou B1)"),
    (1610.0, 1626.5, "mobile satellites (Iridium, Globalstar)"),
    (1670.0, 1710.0, "meteorological satellite downlinks (GOES)"),
    (2180.0, 2200.0, "mobile satellite downlinks"),
    (2320.0, 2345.0, "satellite radio downlinks"),
    (2483.5, 2500.0, "mobile satellite downlinks (Globalstar)"),
    (3700.0, 4200.0, "fixed satellite downlinks (C band)"),
    (7250.0, 7750.0, "fixed satellite downlinks (X band)"),
    (10700.0, 12750.0, "fixed and broadcasting satellite downlinks (Ku band)"),
)
RADIO_QUERY = "breakthrough_listen"
COVERAGE_CATEGORIES = ("no_coverage", "not_downloaded", "no_comparison", "failed", "not_fetched")

ProgressCallback = Callable[[int, int, str, str], None]

_SIGPROC_INT = {"machine_id", "telescope_id", "data_type", "nchans", "nbits", "nifs", "nbeams", "ibeam", "barycentric",
                "pulsarcentric", "nsamples"}
_SIGPROC_DOUBLE = {"tstart", "tsamp", "fch1", "foff", "refdm", "az_start", "za_start", "src_raj", "src_dej", "period"}
_SIGPROC_STR = {"source_name", "rawdatafile"}


def rfi_allocation(freq_mhz: float) -> str | None:
    """Name of the transmitter allocation a frequency falls in, or None."""
    for lo, hi, name in RFI_ALLOCATIONS_MHZ:
        if lo <= freq_mhz <= hi:
            return name
    return None


def max_in_beam_freq_mhz(telescope: str, sep_arcsec: float) -> float | None:
    """Highest frequency at which a position `sep_arcsec` off the pointing is still inside the half-power beam."""
    dish = bl.DISH_DIAMETER_M.get(str(telescope).upper())
    if dish is None:
        return None
    if sep_arcsec <= 0.0:
        return float("inf")
    sep_rad = np.radians(sep_arcsec / 3600.0)
    return float(0.5 * bl.BEAM_FACTOR * bl.SPEED_OF_LIGHT_M_S / (dish * sep_rad) / 1e6)


# ---------------------------------------------------------------------------------------------------- scan files


def read_sigproc_header(path: Path) -> tuple[dict[str, Any], int]:
    """Header keywords of a SIGPROC filterbank file and the byte offset of its data."""
    with open(path, "rb") as fh:
        def _string() -> str:
            raw = fh.read(4)
            if len(raw) != 4:
                raise ArchiveAccessError(f"{path.name}: filterbank header ends early", retryable=False)
            (n,) = struct.unpack("<i", raw)
            if not 1 <= n <= 80:
                raise ArchiveAccessError(f"{path.name}: not a SIGPROC filterbank header", retryable=False)
            return fh.read(n).decode("ascii", "replace")

        if _string() != "HEADER_START":
            raise ArchiveAccessError(f"{path.name}: not a SIGPROC filterbank file", retryable=False)
        hdr: dict[str, Any] = {}
        while True:
            key = _string()
            if key == "HEADER_END":
                break
            if key in _SIGPROC_INT:
                hdr[key] = struct.unpack("<i", fh.read(4))[0]
            elif key in _SIGPROC_DOUBLE:
                hdr[key] = struct.unpack("<d", fh.read(8))[0]
            elif key in _SIGPROC_STR:
                hdr[key] = _string()
            else:
                raise ArchiveAccessError(f"{path.name}: unknown filterbank header keyword '{key}'", retryable=False)
        return hdr, fh.tell()


def _plain(v: Any) -> Any:
    if isinstance(v, bytes):
        return v.decode("ascii", "replace")
    if isinstance(v, np.ndarray):
        return _plain(v.item()) if v.size == 1 else v.tolist()
    return v.item() if isinstance(v, np.generic) else v


class ScanFile:
    """One fine-resolution scan (HDF5 or SIGPROC filterbank), read by channel range in ascending frequency."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self._h5: Any = None
        name = self.path.name.lower()
        if not self.path.is_file():
            raise ArchiveAccessError(f"{self.path.name}: file not found", retryable=False)
        if name.endswith((".h5", ".hdf", ".hdf5")):
            try:
                import h5py
            except ImportError as exc:
                raise ArchiveAccessError("reading Breakthrough Listen HDF5 files needs h5py and hdf5plugin (pip install 'oddstar[radio]')",
                                         retryable=False) from exc
            try:
                import hdf5plugin  # noqa: F401  (registers the bitshuffle filter the archive compresses with)
            except ImportError:
                pass  # an uncompressed file still reads; a compressed one fails below with the filter named
            try:
                self._h5 = h5py.File(self.path, "r")
                self._data = self._h5["data"]
                hdr = {k: _plain(v) for k, v in self._data.attrs.items()}
            except (OSError, KeyError) as exc:
                raise ArchiveAccessError(f"{self.path.name}: not a readable filterbank HDF5 file ({str(exc)[:120]})", retryable=False) from exc
        elif name.endswith(".fil"):
            hdr, offset = read_sigproc_header(self.path)
            if int(hdr.get("nbits", 0)) != 32:
                raise ArchiveAccessError(f"{self.path.name}: {hdr.get('nbits')}-bit filterbank data is not supported (32-bit only)", retryable=False)
            nifs, nchans = int(hdr.get("nifs", 1)), int(hdr.get("nchans", 0))
            if nifs < 1 or nchans < 1:
                raise ArchiveAccessError(f"{self.path.name}: filterbank header has no channels", retryable=False)
            n_ints = (self.path.stat().st_size - offset) // (4 * nifs * nchans)
            if n_ints < 1:
                raise ArchiveAccessError(f"{self.path.name}: filterbank file holds no complete spectrum", retryable=False)
            self._data = np.memmap(self.path, dtype="<f4", mode="r", offset=offset, shape=(n_ints, nifs, nchans))
        else:
            raise ArchiveAccessError(f"{self.path.name}: unknown scan file type (expected .h5, .hdf or .fil)", retryable=False)
        try:
            self.fch1_mhz, self.foff_mhz = float(hdr["fch1"]), float(hdr["foff"])
            self.tsamp_s, self.tstart_mjd = float(hdr["tsamp"]), float(hdr["tstart"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ArchiveAccessError(f"{self.path.name}: header lacks fch1, foff, tsamp or tstart", retryable=False) from exc
        self.source_name = str(hdr.get("source_name", ""))
        if self._data.ndim != 3:
            raise ArchiveAccessError(f"{self.path.name}: data has {self._data.ndim} dimensions, expected (time, feed, channel)", retryable=False)
        self.n_ints, self.nchans = int(self._data.shape[0]), int(self._data.shape[2])
        if self.foff_mhz == 0.0 or self.tsamp_s <= 0.0 or self.n_ints < 2:
            raise ArchiveAccessError(f"{self.path.name}: header values cannot describe a spectrogram", retryable=False)
        self.df_hz = abs(self.foff_mhz) * 1e6
        self.f_low_mhz = min(self.fch1_mhz, self.fch1_mhz + (self.nchans - 1) * self.foff_mhz)

    def freq_of_channel(self, a: float) -> float:
        """Centre frequency (MHz) of ascending channel index `a`."""
        return self.f_low_mhz + a * abs(self.foff_mhz)

    def channel_of_freq(self, freq_mhz: float) -> int:
        return int(round((freq_mhz - self.f_low_mhz) / abs(self.foff_mhz)))

    def read_ascending(self, a0: int, a1: int) -> np.ndarray:
        """Power of ascending channels [a0, a1) as a (time, channel) float32 array."""
        if not 0 <= a0 < a1 <= self.nchans:
            raise ArchiveAccessError(f"{self.path.name}: channel range {a0}:{a1} is outside 0:{self.nchans}", retryable=False)
        try:
            if self.foff_mhz > 0:
                arr = self._data[:, 0, a0:a1]
            else:
                arr = self._data[:, 0, self.nchans - a1:self.nchans - a0][:, ::-1]
        except OSError as exc:  # a missing decompression filter or a truncated file
            raise ArchiveAccessError(f"{self.path.name}: cannot read channels {a0}:{a1} ({str(exc)[:160]})", retryable=False) from exc
        return np.ascontiguousarray(arr, dtype=np.float32)

    def close(self) -> None:
        if self._h5 is not None:
            self._h5.close()
            self._h5 = None
        self._data = None


# ---------------------------------------------------------------------------------------------------- the search


def _settings(cfg: PipelineConfig) -> dict[str, float]:
    ch = cfg.channels["radio_narrowband"]
    return {"window": int(ch["live_window_channels"]), "block": int(ch["live_block_channels"]), "prefilter": float(ch["live_prefilter_snr"]),
            "max_listed": int(ch["live_max_listed_windows"]), "max_drift": float(ch["max_abs_drift_hz_s"]), "min_snr": float(ch["min_snr"])}


def prefilter_block(block: np.ndarray, dt_s: float, df_hz: float, max_abs_drift: float, window: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Strongest drift-summed excess of every `window`-wide stretch of a (time, channel) block.

    The statistic is the channel's: power summed along each trial drift, less the median of the time-averaged
    spectrum, over 1.4826 times its median absolute deviation, both taken per window. Returns per window
    (snr, drift_hz_s, channel index in the block). Channels beyond the last whole window are not used.
    """
    n_t, n_f = block.shape
    n_sub = n_f // window
    n_use = n_sub * window
    if n_sub < 1:
        raise ValueError(f"block of {n_f} channels is narrower than one window of {window}")
    b = block[:, :n_use]
    sub = b.mean(axis=0, dtype=np.float64).reshape(n_sub, window)
    med = np.median(sub, axis=1)
    sig = np.maximum(1e-6, 1.4826 * np.median(np.abs(sub - med[:, None]), axis=1))
    rows = np.arange(n_sub)
    best_snr = np.zeros(n_sub)
    best_drift = np.zeros(n_sub)
    best_ch = rows * window + window // 2
    t_rel = np.arange(n_t, dtype=float) * dt_s
    acc = np.empty(n_use, dtype=np.float64)
    for d in np.linspace(-max_abs_drift, max_abs_drift, N_DRIFT_TRIALS):
        shifts = np.round((d * t_rel) / df_hz).astype(int)
        acc[:] = 0.0
        for ti, sh in enumerate(shifts):
            s = int(sh) % n_use  # acc[i] += row[(i + sh) mod n], the channel's np.roll(row, -sh)
            if s == 0:
                acc += b[ti]
            else:
                acc[:n_use - s] += b[ti, s:]
                acc[n_use - s:] += b[ti, :s]
        snr = (acc.reshape(n_sub, window) / n_t - med[:, None]) / sig[:, None]
        idx = np.argmax(snr, axis=1)
        val = snr[rows, idx]
        better = val > best_snr
        best_snr[better] = val[better]
        best_drift[better] = d
        best_ch[better] = rows[better] * window + idx[better]
    return best_snr, best_drift, best_ch


def _window_key(hit: dict[str, Any]) -> tuple[int, float]:
    """Order of 'hardest to explain': flagged outside the allocated bands, flagged, ON-only excess, any excess."""
    if hit["flagged"]:
        return (3 if not hit["in_rfi_band"] else 2, hit["min_on_snr"])
    return (1, hit["raw_statistic"]) if hit["raw_statistic"] > 0 else (0, hit["max_on_snr"])


def search_scan_set(
    source_id: str,
    on: ScanFile,
    offs: list[ScanFile],
    channel: Any,
    cfg: PipelineConfig,
    *,
    max_freq_mhz: float | None = None,
    progress: Callable[[int, int], None] | None = None,
) -> dict[str, Any]:
    """Search one ON scan against its OFF scans. Returns the counts, the evaluated windows and the selected cadence.

    `channel` is the `RadioNarrowbandChannel`; each candidate window is measured and scored by it. Frequencies above
    `max_freq_mhz` (where the star leaves the half-power beam) are not searched.
    """
    st = _settings(cfg)
    window, block_n = int(st["window"]), int(st["block"])
    if not offs:
        raise ArchiveAccessError("a search needs at least one OFF scan", retryable=False)
    for o in offs:
        if abs(o.df_hz - on.df_hz) > 1e-6 * on.df_hz or abs(o.tsamp_s - on.tsamp_s) > 1e-6 * on.tsamp_s:
            raise ArchiveAccessError(f"{o.path.name} and {on.path.name} differ in channel width or sampling time", retryable=False)
    t_start = time.monotonic()
    order = sorted([(on, "ON"), *[(o, "OFF") for o in offs]], key=lambda p: p[0].tstart_mjd)
    t0 = order[0][0].tstart_mjd
    a_hi = on.nchans
    if max_freq_mhz is not None:
        a_hi = max(0, min(on.nchans, on.channel_of_freq(max_freq_mhz) + 1))
    starts = [a for a in range(0, a_hi, block_n) if min(a + block_n, a_hi) - a >= window]
    counts = {"n_candidate_windows": 0, "n_candidates_merged": 0, "n_windows_evaluated": 0, "n_windows_flagged": 0,
              "n_windows_flagged_outside_rfi_bands": 0, "n_blocks_without_off": 0, "n_channels_searched": 0}
    windows: list[dict[str, Any]] = []
    best: tuple[tuple[int, float], dict[str, Any], dict[str, Any]] | None = None
    strongest: tuple[float, int] | None = None  # (prefilter snr, ascending channel) over the whole scan

    def _evaluate(a_lo: int, powers: list[np.ndarray], pre_snr: float, pre_ch: int) -> None:
        nonlocal best
        cadence = {
            "source_id": source_id, "dt_s": on.tsamp_s, "df_hz": on.df_hz, "center_freq_mhz": round(on.freq_of_channel(a_lo + window / 2.0), 6),
            "scans": [{"scan_index": i, "kind": kind, "t_start_s": round((sf.tstart_mjd - t0) * 86400.0, 3), "power": p}
                      for i, ((sf, kind), p) in enumerate(zip(order, powers, strict=True))],
        }
        m = channel.measure(cadence, {})
        sc = channel.score(m)
        freq = on.freq_of_channel(pre_ch)
        allocation = rfi_allocation(freq)
        hit = {"freq_mhz": round(freq, 6), "window_center_freq_mhz": cadence["center_freq_mhz"], "prefilter_snr": round(float(pre_snr), 2),
               "max_on_snr": round(float(m.values["max_on_snr"]), 2), "min_on_snr": round(float(m.values["min_on_snr"]), 2),
               "off_scan_hits": int(m.values["off_scan_hits"]), "on_persistence": float(m.values["on_persistence"]),
               "drift_hz_s": float(m.values["drift_hz_s"]), "raw_statistic": float(sc.raw_statistic), "flagged": bool(sc.flagged),
               "in_rfi_band": allocation is not None, "rfi_allocation": allocation}
        counts["n_windows_evaluated"] += 1
        if hit["flagged"]:
            counts["n_windows_flagged"] += 1
            counts["n_windows_flagged_outside_rfi_bands"] += 0 if hit["in_rfi_band"] else 1
        windows.append(hit)
        key = _window_key(hit)
        if best is None or key > best[0]:
            for sc_rec in cadence["scans"]:
                sc_rec["power"] = np.asarray(sc_rec["power"], dtype=float).tolist()
            best = (key, hit, cadence)

    def _off_range(o: ScanFile, a0: int, n: int) -> tuple[int, int] | None:
        o0 = o.channel_of_freq(on.freq_of_channel(a0))
        return (o0, o0 + n) if 0 <= o0 and o0 + n <= o.nchans else None

    for bi, a0 in enumerate(starts):
        a1 = min(a0 + block_n, a_hi)
        on_block = on.read_ascending(a0, a1)
        snr, _, ch = prefilter_block(on_block, on.tsamp_s, on.df_hz, st["max_drift"], window)
        n_use = (on_block.shape[1] // window) * window
        top = int(np.argmax(snr))
        if strongest is None or snr[top] > strongest[0]:
            strongest = (float(snr[top]), a0 + int(ch[top]))
        cand = [int(i) for i in np.argsort(-snr) if snr[i] >= st["prefilter"]]
        counts["n_candidate_windows"] += len(cand)
        ranges = [_off_range(o, a0, on_block.shape[1]) for o in offs]
        if any(r is None for r in ranges):
            counts["n_blocks_without_off"] += 1  # an OFF scan does not hold these frequencies: nothing here can be compared
        else:
            counts["n_channels_searched"] += n_use
            if cand:
                by_file = {id(on): on_block}
                for o, r in zip(offs, ranges, strict=True):
                    by_file[id(o)] = o.read_ascending(r[0], r[1])  # type: ignore[index]
                done: list[int] = []
                for i in cand:
                    c = int(ch[i])
                    if any(lo <= c < lo + window for lo in done):
                        counts["n_candidates_merged"] += 1
                        continue
                    lo = max(0, min(c - window // 2, n_use - window))
                    done.append(lo)
                    _evaluate(a0 + lo, [by_file[id(sf)][:, lo:lo + window] for sf, _ in order], float(snr[i]), a0 + c)
        if progress is not None:
            progress(bi + 1, len(starts))

    if best is None and strongest is not None and counts["n_channels_searched"]:
        # Nothing reached the prefilter level: measure the strongest window anyway so the star has a real number
        c = strongest[1]
        lo = max(0, min(c - window // 2, a_hi - window))
        ranges = [_off_range(o, lo, window) for o in offs]
        if all(r is not None for r in ranges):
            by_file = {id(on): on.read_ascending(lo, lo + window)}
            for o, r in zip(offs, ranges, strict=True):
                by_file[id(o)] = o.read_ascending(r[0], r[1])  # type: ignore[index]
            _evaluate(lo, [by_file[id(sf)] for sf, _ in order], strongest[0], c)

    windows.sort(key=lambda h: _window_key(h), reverse=True)
    return {
        "version": SEARCH_VERSION, "source_id": source_id, "on_file": on.path.name, "off_files": [o.path.name for o in offs],
        "on_target": on.source_name, "off_targets": [o.source_name for o in offs], "df_hz": on.df_hz, "dt_s": on.tsamp_s, "n_ints": on.n_ints,
        "freq_lo_mhz": round(on.freq_of_channel(0), 6), "freq_hi_mhz": round(on.freq_of_channel(max(0, a_hi - 1)), 6),
        "max_in_beam_freq_mhz": max_freq_mhz, "n_blocks": len(starts), **counts,
        "n_windows_listed": min(len(windows), int(st["max_listed"])), "windows": windows[:int(st["max_listed"])],
        "selected": None if best is None else {"hit": best[1], "cadence": best[2]},
        "elapsed_s": round(time.monotonic() - t_start, 1),
    }


# ---------------------------------------------------------------------------------------------------- the survey


def local_path(data_dir: Path, url: str) -> Path:
    return Path(data_dir) / url.rsplit("/", 1)[-1]


def plan_files(plan: dict[str, Any]) -> tuple[list[dict[str, Any]], list[list[dict[str, Any]]]]:
    """(ON files, [files of each OFF scan]) of the comparison a plan names."""
    if plan["comparison"] == "simultaneous_beams":
        offs = plan["off_simultaneous"]
    else:
        offs = [e for e in (plan["off_before"], plan["off_after"]) if e is not None]
    return plan["on_files"], [e["files"] for e in offs]


def plan_on_disk(plan: dict[str, Any], data_dir: Path) -> bool:
    """Every file of the plan's comparison is in `data_dir` with the size the index lists."""
    on_files, off_sets = plan_files(plan)
    for f in [*on_files, *[f for s in off_sets for f in s]]:
        p = local_path(data_dir, f["url"])
        if not p.is_file() or (f["size_bytes"] and p.stat().st_size != f["size_bytes"]):
            return False
    return True


def radio_search_status(star: dict[str, Any]) -> str:
    """ "searched" or the reason a star was not: one of `COVERAGE_CATEGORIES`."""
    status = str((star.get("radio_fetch") or {}).get("status") or "")
    if status == "ok" and star.get("scans"):
        return "searched"
    if status.startswith("failed"):
        return "failed"
    return status if status in COVERAGE_CATEGORIES else "no_coverage"


def search_coverage_summary(stars: list[dict[str, Any]], cfg: PipelineConfig) -> dict[str, Any]:
    """Counts of searched and not-searched stars, in the layout the other survey channels use."""
    ids: dict[str, list[str]] = {cat: [] for cat in COVERAGE_CATEGORIES}
    n_searched = 0
    for s in stars:
        cat = radio_search_status(s)
        if cat == "searched":
            n_searched += 1
        else:
            ids[cat].append(str(s["source_id"]))
    out: dict[str, Any] = {"n_input": len(stars), "n_searched": n_searched, "n_not_searched": len(stars) - n_searched}
    for cat, lst in ids.items():
        out[f"n_{cat}"] = len(lst)
    out["not_searched_ids"] = {cat: lst for cat, lst in ids.items() if lst}
    return out


def _search_plan(source_id: str, plan: dict[str, Any], data_dir: Path, channel: Any, cfg: PipelineConfig, cache: RawCache,
                 progress: Callable[[int, int], None] | None) -> list[dict[str, Any]]:
    """Search results of one plan, one per ON file (frequency band). Answers are cached by files and settings."""
    on_files, off_sets = plan_files(plan)
    results: list[dict[str, Any]] = []
    for on_f in on_files:
        off_fs = [next(f for f in s if abs(f["center_freq_mhz"] - on_f["center_freq_mhz"]) <= 1.0) for s in off_sets]
        max_freq = max_in_beam_freq_mhz(plan["telescope"], plan["sep_arcsec"])
        desc = json.dumps({"version": SEARCH_VERSION, "source_id": source_id, "on": [on_f["url"], on_f["size_bytes"]],
                           "offs": [[f["url"], f["size_bytes"]] for f in off_fs], "max_freq_mhz": max_freq, "settings": _settings(cfg)},
                          sort_keys=True)
        sig = hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]
        cached = cache.lookup(RADIO_QUERY, SEARCH_VERSION, sig)
        if cached is not None and not cached.offline_fixture:
            try:
                results.append(json.loads(cache.read_bytes(cached).decode("utf-8")))
                continue
            except (FileNotFoundError, OSError, ValueError):
                pass
        on = ScanFile(local_path(data_dir, on_f["url"]))
        offs: list[ScanFile] = []
        try:
            offs = [ScanFile(local_path(data_dir, f["url"])) for f in off_fs]
            res = search_scan_set(source_id, on, offs, channel, cfg, max_freq_mhz=max_freq, progress=progress)
        finally:
            on.close()
            for o in offs:
                o.close()
        res.update({"telescope": plan["telescope"], "sep_arcsec": plan["sep_arcsec"], "comparison": plan["comparison"], "on_mjd": plan["mjd"]})
        cache.put(archive=RADIO_QUERY, release=SEARCH_VERSION, request_signature=sig, query_text=desc,
                  payload=json.dumps(res, sort_keys=True).encode("utf-8"), row_count=int(res["n_windows_evaluated"]), offline_fixture=False)
        results.append(res)
    return results


def enrich_sample_radio(
    stars: list[dict[str, Any]],
    cfg: PipelineConfig,
    cache: RawCache,
    channel: Any,
    data_dir: Path,
    *,
    max_stars: int | None = None,
    progress: ProgressCallback | None = None,
    block_progress: Callable[[str, int, int], None] | None = None,
    sleep: Callable[[float], Any] = time.sleep,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Radio cadence of every star that has a complete scan set on disk. Returns (stars, summary).

    The index is asked for every star (cached, so a repeat is fast); a star inside a beam gets its scan plan; every
    complete plan whose files are all in `data_dir` is searched. No radio data is downloaded here.
    """
    t_start = time.monotonic()
    n_limit = len(stars) if max_stars is None else max(0, min(int(max_stars), len(stars)))
    out: list[dict[str, Any]] = []
    for i, star in enumerate(stars):
        s = dict(star)
        sid = str(s["source_id"])
        if i >= n_limit:
            s["radio_fetch"] = {"status": "not_fetched"}
            out.append(s)
            continue
        ra, dec = float(s["ra_deg"]), float(s["dec_deg"])
        qs: dict[str, str] = s.setdefault("vetting_archive", {}).setdefault("query_status", {})
        rec: dict[str, Any] = {}
        t_req = time.monotonic()
        hit = True
        try:
            rows, _, hit = bl.query_bl_files(ra, dec, cache, sleep=sleep)
            cov = bl.beam_coverage(ra, dec, rows)
            if cov["status"] != "in_beam":
                status = "no_coverage"
                rec["index_status"] = cov["status"]
            else:
                plans = [p for p in bl.plan_star_scans(ra, dec, rows, cache, sleep=sleep) if p["complete"]]
                rec["n_complete_scan_sets"] = len(plans)
                ready = [p for p in plans if plan_on_disk(p, data_dir)]
                if not plans:
                    status = "no_comparison"
                elif not ready:
                    status = "not_downloaded"
                    rec["smallest_scan_set_bytes"] = min(p["bytes_on_plus_offs"] for p in plans)
                else:
                    results: list[dict[str, Any]] = []
                    for p in ready:
                        results.extend(_search_plan(sid, p, data_dir, channel, cfg, cache,
                                                    (lambda d, t, _sid=sid: block_progress(_sid, d, t)) if block_progress else None))
                    chosen = max((r for r in results if r["selected"] is not None), key=lambda r: _window_key(r["selected"]["hit"]), default=None)
                    if chosen is None:
                        status = "failed: no window of any scan set could be compared with its OFF scans"
                    else:
                        status = "ok"
                        cad = chosen["selected"]["cadence"]
                        s.update({"scans": cad["scans"], "dt_s": cad["dt_s"], "df_hz": cad["df_hz"], "center_freq_mhz": cad["center_freq_mhz"],
                                  "in_rfi_band": bool(chosen["selected"]["hit"]["in_rfi_band"])})
                        s["radio_search"] = {
                            "n_scan_sets_searched": len(ready), "n_scan_sets_not_downloaded": len(plans) - len(ready),
                            "selected_hit": chosen["selected"]["hit"], "selected_on_file": chosen["on_file"],
                            "scan_sets": [{k: v for k, v in r.items() if k != "selected"} for r in results],
                        }
        except ArchiveAccessError as exc:
            status = f"failed: {str(exc)[:200]}"
        rec["status"] = status
        s["radio_fetch"] = rec
        qs[RADIO_QUERY] = "ok" if status == "ok" else status
        out.append(s)
        if progress is not None:
            progress(i + 1, n_limit, sid, status)
        if not hit:
            wait = bl.MIN_REQUEST_INTERVAL_S - (time.monotonic() - t_req)
            if wait > 0:
                sleep(wait)
    cov_sum = search_coverage_summary(out, cfg)
    searched = [s for s in out if radio_search_status(s) == "searched"]
    summary = {
        "survey": "Breakthrough Listen open data, fine-resolution scans on local disk", "data_dir": str(data_dir),
        "n_requested": n_limit, "n_ok": cov_sum["n_searched"], **{f"n_{c}": cov_sum[f"n_{c}"] for c in COVERAGE_CATEGORIES},
        "failed_ids": cov_sum["not_searched_ids"].get("failed", []), "not_downloaded_ids": cov_sum["not_searched_ids"].get("not_downloaded", []),
        "n_windows_evaluated": sum(r["n_windows_evaluated"] for s in searched for r in s["radio_search"]["scan_sets"]),
        "n_windows_flagged": sum(r["n_windows_flagged"] for s in searched for r in s["radio_search"]["scan_sets"]),
        "n_windows_flagged_outside_rfi_bands": sum(r["n_windows_flagged_outside_rfi_bands"] for s in searched for r in s["radio_search"]["scan_sets"]),
        "elapsed_s": round(time.monotonic() - t_start, 1),
    }
    return out, summary


# ---------------------------------------------------------------------------------------------------- download


def download_file(url: str, dest: Path, size_bytes: int, md5sum: str, *, chunk_bytes: int = 8 << 20,
                  progress: Callable[[int, int], None] | None = None, timeout_s: float = 60.0) -> str:
    """Fetch one archive file to `dest`, resuming a partial download, and check its size and MD5.

    Returns "present" (already complete), or "downloaded". The partial file is `dest` + ".part".
    """
    import requests

    dest = Path(dest)
    if dest.is_file() and (not size_bytes or dest.stat().st_size == size_bytes):
        return "present"
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    have = part.stat().st_size if part.is_file() else 0
    if size_bytes and have > size_bytes:
        part.unlink()
        have = 0
    if not size_bytes or have < size_bytes:
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with requests.get(url, headers=headers, stream=True, timeout=timeout_s) as resp:
                if resp.status_code == 200 and have:
                    have = 0  # the server ignored the range: start again
                elif resp.status_code not in (200, 206):
                    raise ArchiveAccessError(f"{url}: HTTP {resp.status_code}", retryable=resp.status_code >= 500)
                with open(part, "ab" if have else "wb") as fh:
                    for chunk in resp.iter_content(chunk_size=chunk_bytes):
                        fh.write(chunk)
                        have += len(chunk)
                        if progress is not None:
                            progress(have, size_bytes)
        except requests.exceptions.RequestException as exc:
            raise ArchiveAccessError(f"{url}: {type(exc).__name__}: {str(exc)[:160]} (rerun to resume from {have} bytes)", retryable=True) from exc
    if size_bytes and part.stat().st_size != size_bytes:
        raise ArchiveAccessError(f"{url}: got {part.stat().st_size} of {size_bytes} bytes (rerun to resume)", retryable=True)
    if md5sum and set(md5sum) != {"0"}:
        h = hashlib.md5()  # noqa: S324  (the archive publishes MD5; this is a transfer check, not security)
        with open(part, "rb") as fh:
            for chunk in iter(lambda: fh.read(chunk_bytes), b""):
                h.update(chunk)
        if h.hexdigest() != md5sum.lower():
            raise ArchiveAccessError(f"{url}: MD5 {h.hexdigest()} differs from the index value {md5sum}; the partial file is kept as {part.name}",
                                     retryable=False)
    part.rename(dest)
    return "downloaded"


def main(argv: list[str] | None = None) -> int:
    """Download the files of one star's scan set, as listed by `python -m oddstar.io.bl_opendata --plan`."""
    from oddstar.config import load_config, repo_root

    ap = argparse.ArgumentParser(prog="python -m oddstar.io.bl_radio", description=main.__doc__)
    ap.add_argument("--plan", required=True, help="bl_scan_plan.json written by `python -m oddstar.io.bl_opendata --plan`")
    ap.add_argument("--source-id", required=True, help="star whose scan set to fetch")
    ap.add_argument("--scan-index", type=int, default=None,
                    help="index of the ON scan in the star's plan (default: the complete scan set with the fewest bytes)")
    ap.add_argument("--data-dir", default=None, help="where the files go (default: breakthrough_listen_scans beside the raw cache)")
    ap.add_argument("--dry-run", action="store_true", help="list the files and sizes, download nothing")
    args = ap.parse_args(argv)
    root = repo_root()
    cfg, _, _ = load_config(root=root)
    data_dir = Path(args.data_dir) if args.data_dir else default_data_dir(root, cfg)
    stars = {str(e["source_id"]): e for e in json.loads(Path(args.plan).read_text(encoding="utf-8"))["stars"]}
    entry = stars.get(str(args.source_id))
    if entry is None or not entry.get("scans"):
        print(f"no scan plan for {args.source_id} in {args.plan}", file=sys.stderr)
        return 2
    if args.scan_index is not None:
        plan = entry["scans"][args.scan_index]
    else:
        complete = [p for p in entry["scans"] if p["complete"]]
        if not complete:
            print(f"{args.source_id} has no complete scan set", file=sys.stderr)
            return 2
        plan = min(complete, key=lambda p: p["bytes_on_plus_offs"])
    on_files, off_sets = plan_files(plan)
    files = [*on_files, *[f for fs in off_sets for f in fs]]
    total = sum(f["size_bytes"] for f in files)
    listing = [{"file": local_path(data_dir, f["url"]).name, "size_bytes": f["size_bytes"], "url": f["url"]} for f in files]
    if args.dry_run:
        print(json.dumps({"data_dir": str(data_dir), "target": plan["target"], "comparison": plan["comparison"], "total_bytes": total,
                          "files": listing}, indent=2))
        return 0
    done: list[dict[str, Any]] = []
    for f, row in zip(files, listing, strict=True):
        last = [0.0]

        def _progress(have: int, size: int, _name: str = row["file"], _last: list[float] = last) -> None:
            if time.monotonic() - _last[0] >= 30.0:
                _last[0] = time.monotonic()
                print(f"  {_name}: {have / 1e9:.2f} of {size / 1e9:.2f} GB", file=sys.stderr)

        outcome = download_file(f["url"], local_path(data_dir, f["url"]), f["size_bytes"], f["md5sum"], progress=_progress)
        print(f"  {row['file']}: {outcome}", file=sys.stderr)
        done.append({**row, "outcome": outcome})
    print(json.dumps({"data_dir": str(data_dir), "target": plan["target"], "comparison": plan["comparison"], "total_bytes": total,
                      "files": done}, indent=2))
    return 0


def default_data_dir(root: Path, cfg: PipelineConfig) -> Path:
    """Where downloaded scans live: beside the raw cache, not inside it (the cache evicts files above its size limit)."""
    return (Path(root) / cfg.storage.raw_dir).parent / "breakthrough_listen_scans"


__all__ = [
    "COVERAGE_CATEGORIES", "RFI_ALLOCATIONS_MHZ", "ScanFile", "default_data_dir", "download_file", "enrich_sample_radio", "local_path",
    "max_in_beam_freq_mhz", "plan_files", "plan_on_disk", "prefilter_block", "radio_search_status", "read_sigproc_header", "rfi_allocation",
    "search_coverage_summary", "search_scan_set",
]

if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
