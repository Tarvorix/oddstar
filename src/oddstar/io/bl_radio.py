"""Narrowband search of Breakthrough Listen fine-resolution scans held on local disk.

The radio channel measures one small time x frequency block per scan. A fine-resolution scan is one block of 16
spectra by several hundred million channels, so this module does the part the channel cannot: it walks the ON
scan in blocks of `live_block_channels`, finds every `live_window_channels`-wide window whose drift-summed power
stands out (`prefilter_block`, the channel's own statistic and drift grid on a per-window noise estimate), cuts the
same frequencies from the OFF scans, and hands each such window to the channel as an ordinary ON/OFF cadence.
Every number that decides a flag comes from `RadioNarrowbandChannel.measure` and `.score`; nothing here re-decides.

What a star's record holds after `enrich_sample_radio`:

- `radio_fetch.status`: "ok" (a complete scan set was on disk, or already searched and its result stored, and was
  searched), "not_downloaded" (the star is inside a beam and a complete ON/OFF set exists in the archive, but its
  files are not in the data directory and no stored result exists),
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
import shutil
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

SEARCH_VERSION = "bl_local_search_v4"
# Version of the cutting phases (candidate windows of the ON scan, the same windows of the other scans). The windows
# saved under `<data_dir>/_work` are keyed by it, not by SEARCH_VERSION, so a change that only touches the scoring
# (the allocation table, the listing) scores the saved windows again without fetching a scan.
CUT_VERSION = "bl_cut_v1"
# Trial drift rates per search; must equal the grid inside `channels.radio_narrowband.dedoppler_search_scan`
N_DRIFT_TRIALS = 161
# Bands allocated to transmitters, where a narrowband line toward any sky position is expected from them.
# From 698 to 4200 MHz the table is complete: every band the United States table of frequency allocations (47 CFR
# 2.106, federal and non-federal columns) gives to a transmitting service is listed, so the only frequencies in that
# range that are not listed are the passive bands 1400 to 1427, 1660.5 to 1668.4 and 2690 to 2700 MHz (radio
# astronomy, passive sensing); from 2900 to 4200 MHz no band is passive. Above 4200 MHz only two satellite downlink
# bands are listed; the table is not complete there and must be completed before scans of those receivers are
# searched. History of the table and the reason for the rule: configs/CHANGELOG.md.
# (low MHz, high MHz, service)
ALLOCATIONS_COMPLETE_MHZ = (698.0, 4200.0)
PASSIVE_BANDS_MHZ: tuple[tuple[float, float], ...] = ((1400.0, 1427.0), (1660.5, 1668.4), (2690.0, 2700.0))
RFI_ALLOCATIONS_MHZ: tuple[tuple[float, float, str], ...] = (
    (698.0, 806.0, "mobile telephone and public safety (700 MHz band)"),
    (806.0, 902.0, "land mobile, mobile telephone and air-ground (800 MHz bands)"),
    (902.0, 928.0, "industrial, scientific and medical band; radiolocation"),
    (928.0, 960.0, "paging, fixed and land mobile"),
    (960.0, 1164.0, "aeronautical radionavigation (DME, TACAN, ADS-B)"),
    (1164.0, 1300.0, "radionavigation satellites (GPS L5/L2, GLONASS L3/L2, Galileo E5/E6)"),
    (1300.0, 1350.0, "aeronautical radionavigation (air route surveillance radar)"),
    (1350.0, 1390.0, "radiolocation and GPS L3 (1381.05 MHz)"),
    (1390.0, 1400.0, "fixed, mobile and medical telemetry"),
    (1427.0, 1435.0, "telemetry, fixed and mobile"),
    (1435.0, 1525.0, "aeronautical mobile telemetry"),
    (1525.0, 1559.0, "mobile satellite downlinks (Inmarsat)"),
    (1559.0, 1610.0, "radionavigation satellites (GPS L1, GLONASS L1, Galileo E1, BeiDou B1)"),
    (1610.0, 1626.5, "mobile satellites (Iridium, Globalstar)"),
    (1626.5, 1660.5, "mobile satellite uplinks (Inmarsat)"),
    (1668.4, 1670.0, "meteorological aids (radiosondes)"),
    (1670.0, 1710.0, "meteorological satellite downlinks (GOES)"),
    (1710.0, 1780.0, "mobile telephone uplinks (AWS)"),
    (1780.0, 1850.0, "fixed, mobile and space operation uplinks"),
    (1850.0, 1915.0, "mobile telephone uplinks (PCS)"),
    (1915.0, 1930.0, "unlicensed personal communications (cordless telephones)"),
    (1930.0, 2000.0, "mobile telephone downlinks (PCS)"),
    (2000.0, 2025.0, "mobile and mobile satellite uplinks"),
    (2025.0, 2110.0, "space operation uplinks and television auxiliary links"),
    (2110.0, 2180.0, "mobile telephone downlinks (AWS)"),
    (2180.0, 2200.0, "mobile satellite downlinks"),
    (2200.0, 2290.0, "space operation and satellite telemetry downlinks"),
    (2290.0, 2300.0, "deep space downlinks; fixed and mobile"),
    (2300.0, 2320.0, "amateur and wireless communications"),
    (2320.0, 2345.0, "satellite radio downlinks"),
    (2345.0, 2360.0, "wireless communications"),
    (2360.0, 2400.0, "aeronautical telemetry, medical devices and amateur"),
    (2400.0, 2483.5, "industrial, scientific and medical band (wireless networks)"),
    (2483.5, 2500.0, "mobile satellite downlinks (Globalstar)"),
    (2500.0, 2690.0, "broadband radio service (mobile)"),
    (2700.0, 2900.0, "airport surveillance and weather radar"),
    (2900.0, 3100.0, "maritime radionavigation and radiolocation (ship radar)"),
    (3100.0, 3300.0, "radiolocation (military radar)"),
    (3300.0, 3550.0, "radiolocation, amateur and mobile broadband"),
    (3550.0, 3700.0, "citizens broadband radio service and radiolocation"),
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


def _scan_meta(sf: ScanFile) -> dict[str, Any]:
    return {"file": sf.path.name, "source_name": sf.source_name, "df_hz": sf.df_hz, "tsamp_s": sf.tsamp_s, "tstart_mjd": sf.tstart_mjd,
            "n_ints": sf.n_ints, "nchans": sf.nchans, "f_low_mhz": sf.f_low_mhz, "foff_abs_mhz": abs(sf.foff_mhz)}


def _load_windows(path: Path, n_windows: int, n_ints: int, window: int) -> np.ndarray:
    if n_windows == 0:
        return np.zeros((0, n_ints, window), dtype=np.float32)
    return np.memmap(path, dtype="<f4", mode="r", shape=(n_windows, n_ints, window))


def cut_channels(window: int, max_drift_hz_s: float, span_s: float | None, df_hz: float) -> int:
    """Channels to cut around a hit: `window`, doubled until a signal at the largest drift stays inside over `span_s`."""
    cut = int(window)
    while span_s and max_drift_hz_s * span_s > (cut / 2.0) * df_hz and cut < (1 << 20):
        cut *= 2
    return cut


def plan_offset_s(plan: dict[str, Any]) -> float:
    """Largest time between the start of the reference ON scan and the start of another scan of the plan's comparison."""
    return max((abs(float(e["mjd"]) - float(plan["mjd"])) * 86400.0 for _, e in plan_other_scans(plan)), default=0.0)


def find_candidate_windows(
    on: ScanFile, cfg: PipelineConfig, store_path: Path, *, max_freq_mhz: float | None = None,
    progress: Callable[[int, int], None] | None = None, span_s: float | None = None,
) -> dict[str, Any]:
    """Phase 1, ON scan only: every window that reaches the prefilter level, plus the strongest window of each block.

    The power of each listed window is appended to `store_path` (float32, time x channel). The returned record
    names each window by its lowest ascending channel and frequency, so the same sky frequencies can be cut from
    the other scans later without the ON file. `span_s` is the time from the start of this scan to the end of the
    farthest scan it will be compared with: the cut (`cut_channels`) is `live_window_channels` wide or, when a
    signal at the largest drift would leave that over `span_s`, the next power-of-two multiple that holds it.
    """
    st = _settings(cfg)
    window, block_n = int(st["window"]), int(st["block"])
    cut = cut_channels(window, float(st["max_drift"]), span_s, on.df_hz)
    a_hi = on.nchans
    if max_freq_mhz is not None:
        a_hi = max(0, min(on.nchans, on.channel_of_freq(max_freq_mhz) + 1))
    starts = [a for a in range(0, a_hi, block_n) if min(a + block_n, a_hi) - a >= cut]
    entries: list[dict[str, Any]] = []
    blocks: list[dict[str, Any]] = []
    n_candidates = n_merged = 0
    store_path = Path(store_path)
    store_path.parent.mkdir(parents=True, exist_ok=True)
    with open(store_path, "wb") as fh:
        for bi, a0 in enumerate(starts):
            a1 = min(a0 + block_n, a_hi)
            block = on.read_ascending(a0, a1)
            snr, _, ch = prefilter_block(block, on.tsamp_s, on.df_hz, st["max_drift"], window)
            n_use = (block.shape[1] // window) * window
            cand = [int(i) for i in np.argsort(-snr) if snr[i] >= st["prefilter"]]
            n_candidates += len(cand)
            picks = [(i, True) for i in cand] or [(int(np.argmax(snr)), False)]  # a quiet block still lists its strongest window
            done: list[int] = []
            for i, is_cand in picks:
                c = int(ch[i])
                if any(lo <= c < lo + cut for lo in done):
                    n_merged += 1
                    continue
                lo = max(0, min(c - cut // 2, n_use - cut))
                done.append(lo)
                fh.write(np.ascontiguousarray(block[:, lo:lo + cut], dtype="<f4").tobytes())
                entries.append({"a_lo": a0 + lo, "f_lo_mhz": on.freq_of_channel(a0 + lo), "block": bi, "candidate": is_cand,
                                "prefilter_snr": float(snr[i]), "hit_channel": a0 + c, "hit_freq_mhz": on.freq_of_channel(a0 + c)})
            blocks.append({"a0": a0, "n_use": n_use, "f_lo_mhz": on.freq_of_channel(a0), "f_hi_mhz": on.freq_of_channel(a0 + n_use - 1)})
            if progress is not None:
                progress(bi + 1, len(starts))
    return {"version": CUT_VERSION, "on": _scan_meta(on), "window": cut, "prefilter_window": window, "block_channels": block_n,
            "max_in_beam_freq_mhz": max_freq_mhz,
            "freq_lo_mhz": round(on.freq_of_channel(0), 6), "freq_hi_mhz": round(on.freq_of_channel(max(0, a_hi - 1)), 6),
            "blocks": blocks, "entries": entries, "n_candidate_windows": n_candidates, "n_candidates_merged": n_merged}


def extract_windows(scan: ScanFile, found: dict[str, Any], store_path: Path, *, progress: Callable[[int, int], None] | None = None) -> dict[str, Any]:
    """Phase 2, one other scan (OFF, or another ON scan of a cadence): the power at the sky frequencies of every listed
    window of the reference ON scan, appended to `store_path`.

    A window the scan does not hold is written as zeros and marked not covered; so is a whole block.
    """
    on = found["on"]
    if abs(scan.df_hz - on["df_hz"]) > 1e-6 * on["df_hz"] or abs(scan.tsamp_s - on["tsamp_s"]) > 1e-6 * on["tsamp_s"]:
        raise ArchiveAccessError(f"{scan.path.name} and {on['file']} differ in channel width or sampling time", retryable=False)
    window = int(found["window"])
    covered: list[bool] = []
    blocks_covered: list[bool] = []
    by_block: dict[int, list[dict[str, Any]]] = {}
    for e in found["entries"]:
        by_block.setdefault(int(e["block"]), []).append(e)
    store_path = Path(store_path)
    store_path.parent.mkdir(parents=True, exist_ok=True)
    with open(store_path, "wb") as fh:
        for bi, blk in enumerate(found["blocks"]):
            o0 = scan.channel_of_freq(blk["f_lo_mhz"])
            ok = 0 <= o0 and o0 + int(blk["n_use"]) <= scan.nchans
            blocks_covered.append(bool(ok))
            data = scan.read_ascending(o0, o0 + int(blk["n_use"])) if ok and by_block.get(bi) else None
            for e in by_block.get(bi, []):
                if data is None:
                    fh.write(np.zeros((scan.n_ints, window), dtype="<f4").tobytes())
                    covered.append(False)
                else:
                    lo = int(e["a_lo"]) - int(blk["a0"])
                    fh.write(np.ascontiguousarray(data[:, lo:lo + window], dtype="<f4").tobytes())
                    covered.append(True)
            if progress is not None:
                progress(bi + 1, len(found["blocks"]))
    return {**_scan_meta(scan), "covered": covered, "blocks_covered": blocks_covered}


def evaluate_windows(
    source_id: str, found: dict[str, Any], on_store: Path, offs: list[tuple[dict[str, Any], Path]], channel: Any, cfg: PipelineConfig,
    *, progress: Callable[[int, int], None] | None = None, kinds: list[str] | None = None,
) -> dict[str, Any]:
    """Phase 3: measure and score every candidate window with the channel. Needs only the stored windows.

    `offs` are the other scans of the comparison; `kinds` names each "OFF" (the default) or "ON" (another ON scan
    of a cadence, in which the channel then requires the signal too). `progress(done, total)` is called after each
    scored window.
    """
    st = _settings(cfg)
    kinds = list(kinds) if kinds is not None else ["OFF"] * len(offs)
    if len(kinds) != len(offs) or "OFF" not in kinds:
        raise ArchiveAccessError("a search needs at least one OFF scan", retryable=False)
    on, window, entries = found["on"], int(found["window"]), found["entries"]
    on_w = _load_windows(on_store, len(entries), int(on["n_ints"]), window)
    off_w = [_load_windows(path, len(entries), int(meta["n_ints"]), window) for meta, path in offs]
    order = sorted([(on["tstart_mjd"], "ON", -1), *[(meta["tstart_mjd"], kinds[k], k) for k, (meta, _) in enumerate(offs)]], key=lambda x: x[0])
    t0 = order[0][0]
    block_ok = [all(meta["blocks_covered"][bi] for meta, _ in offs) for bi in range(len(found["blocks"]))]
    counts = {"n_candidate_windows": int(found["n_candidate_windows"]), "n_candidates_merged": int(found["n_candidates_merged"]),
              "n_windows_evaluated": 0, "n_windows_flagged": 0, "n_windows_flagged_outside_rfi_bands": 0,
              "n_blocks_without_off": sum(1 for ok in block_ok if not ok),
              "n_channels_searched": sum(int(b["n_use"]) for b, ok in zip(found["blocks"], block_ok, strict=True) if ok)}
    windows: list[dict[str, Any]] = []
    best: tuple[tuple[int, float], dict[str, Any], dict[str, Any]] | None = None

    def _evaluate(k: int) -> None:
        nonlocal best
        e = entries[k]
        cadence = {
            "source_id": source_id, "dt_s": on["tsamp_s"], "df_hz": on["df_hz"],
            "center_freq_mhz": round(e["f_lo_mhz"] + (window / 2.0) * on["foff_abs_mhz"], 6),
            "scans": [{"scan_index": i, "kind": kind, "t_start_s": round((t - t0) * 86400.0, 3), "power": on_w[k] if j < 0 else off_w[j][k]}
                      for i, (t, kind, j) in enumerate(order)],
        }
        m = channel.measure(cadence, {})
        sc = channel.score(m)
        allocation = rfi_allocation(e["hit_freq_mhz"])
        hit = {"freq_mhz": round(e["hit_freq_mhz"], 6), "window_center_freq_mhz": cadence["center_freq_mhz"], "prefilter_snr": round(e["prefilter_snr"], 2),
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

    comparable = [k for k, e in enumerate(entries) if block_ok[int(e["block"])] and all(meta["covered"][k] for meta, _ in offs)]
    todo = [k for k in comparable if entries[k]["candidate"]]
    for i, k in enumerate(todo, start=1):
        _evaluate(k)
        if progress is not None:
            progress(i, len(todo))
    if best is None and comparable:
        # Nothing reached the prefilter level: measure the strongest window anyway so the star has a real number
        _evaluate(max(comparable, key=lambda k: entries[k]["prefilter_snr"]))

    windows.sort(key=lambda h: _window_key(h), reverse=True)
    off_metas = [meta for (meta, _), kind in zip(offs, kinds, strict=True) if kind == "OFF"]
    on_other = [meta["file"] for (meta, _), kind in zip(offs, kinds, strict=True) if kind == "ON"]
    return {
        "version": SEARCH_VERSION, "source_id": source_id, "on_file": on["file"], "off_files": [meta["file"] for meta in off_metas],
        "on_target": on["source_name"], "off_targets": [meta["source_name"] for meta in off_metas], "df_hz": on["df_hz"], "dt_s": on["tsamp_s"],
        **({"on_other_files": on_other, "n_on_scans": 1 + len(on_other)} if on_other else {}), "window_channels": window,
        "n_ints": on["n_ints"], "freq_lo_mhz": found["freq_lo_mhz"], "freq_hi_mhz": found["freq_hi_mhz"],
        "max_in_beam_freq_mhz": found["max_in_beam_freq_mhz"], "n_blocks": len(found["blocks"]), **counts,
        "n_windows_listed": min(len(windows), int(st["max_listed"])), "windows": windows[:int(st["max_listed"])],
        "selected": None if best is None else {"hit": best[1], "cadence": best[2]},
    }


def search_scan_set(
    source_id: str,
    on: ScanFile,
    offs: list[ScanFile],
    channel: Any,
    cfg: PipelineConfig,
    *,
    max_freq_mhz: float | None = None,
    progress: Callable[[int, int], None] | None = None,
    work_dir: Path | None = None,
    other_ons: list[ScanFile] | None = None,
    offset_s: float | None = None,
) -> dict[str, Any]:
    """Search one ON scan against its OFF scans. Returns the counts, the evaluated windows and the selected cadence.

    `other_ons` are further ON scans of the same cadence: the signal is then required in each of them too.
    `offset_s` is the largest time between the start of `on` and the start of another scan (from the scan headers
    when not given); with the length of `on` it sets how wide the windows are cut.

    `channel` is the `RadioNarrowbandChannel`; each candidate window is measured and scored by it. Frequencies above
    `max_freq_mhz` (where the star leaves the half-power beam) are not searched. The three phases
    (`find_candidate_windows`, `extract_windows`, `evaluate_windows`) can also be run one file at a time when the
    disk cannot hold a whole scan set.
    """
    import tempfile

    if not offs:
        raise ArchiveAccessError("a search needs at least one OFF scan", retryable=False)
    t_start = time.monotonic()
    with tempfile.TemporaryDirectory(dir=work_dir) as tmp:
        on_store = Path(tmp) / "on.f32"
        others = [*[("ON", o) for o in (other_ons or [])], *[("OFF", o) for o in offs]]
        if offset_s is None:
            offset_s = max(abs(o.tstart_mjd - on.tstart_mjd) * 86400.0 for _, o in others)
        found = find_candidate_windows(on, cfg, on_store, max_freq_mhz=max_freq_mhz, progress=progress,
                                       span_s=offset_s + on.n_ints * on.tsamp_s)
        recs = []
        for k, (_, o) in enumerate(others):
            path = Path(tmp) / f"other{k}.f32"
            recs.append((extract_windows(o, found, path), path))
        res = evaluate_windows(source_id, found, on_store, recs, channel, cfg, kinds=[kind for kind, _ in others])
    res["elapsed_s"] = round(time.monotonic() - t_start, 1)
    return res


# ---------------------------------------------------------------------------------------------------- the survey


def local_path(data_dir: Path, url: str) -> Path:
    """Where a scan file lives on disk. Green Bank names carry the recording node in front (`blc07_guppi_...`,
    `spliced_blc0001..._guppi_...`). Parkes multibeam names start with `guppi_`: the two nodes of one scan write the
    same file name into different directories, so such a name is prefixed with a short hash of its URL and the two
    files stay apart."""
    name = url.rsplit("/", 1)[-1]
    if name.startswith("guppi_"):
        name = f"{hashlib.sha1(_file_key(url).encode('utf-8')).hexdigest()[:8]}_{name}"
    return Path(data_dir) / name


def plan_other_scans(plan: dict[str, Any]) -> list[tuple[str, dict[str, Any]]]:
    """(kind, scan) of every scan the plan's ON scan is compared with: "OFF" scans and, in a cadence, other "ON" scans."""
    if plan["comparison"] == "simultaneous_beams":
        return [("OFF", e) for e in plan["off_simultaneous"]]
    if plan["comparison"] == "cadence":
        return sorted([*[("ON", e) for e in plan["on_other"]], *[("OFF", e) for e in plan["off_cadence"]]], key=lambda x: x[1]["mjd"])
    return [("OFF", e) for e in (plan["off_before"], plan["off_after"]) if e is not None]


def plan_files(plan: dict[str, Any]) -> tuple[list[dict[str, Any]], list[list[dict[str, Any]]]]:
    """(files of the reference ON scan, [files of each other scan]) of the comparison a plan names."""
    return plan["on_files"], [e["files"] for _, e in plan_other_scans(plan)]


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


def _band_sets(plan: dict[str, Any]) -> list[tuple[dict[str, Any], list[dict[str, Any]], list[str]]]:
    """(ON file, the file of each other scan, the kind of each other scan) for every frequency band of a plan."""
    on_files, other_sets = plan_files(plan)
    kinds = [kind for kind, _ in plan_other_scans(plan)]
    return [(on_f, [next(f for f in s if abs(f["center_freq_mhz"] - on_f["center_freq_mhz"]) <= 1.0) for s in other_sets], kinds)
            for on_f in on_files]


def _file_key(url: str) -> str:
    """A file's address without the scheme: the index lists the same host as http and as https."""
    return url.split("://", 1)[-1]


def _band_signature(source_id: str, plan: dict[str, Any], on_f: dict[str, Any], other_fs: list[dict[str, Any]], cfg: PipelineConfig,
                    kinds: list[str] | None = None) -> tuple[str, str, float | None]:
    """Key of one band's stored result. A plan without other ON scans gets the key it had before cadences existed."""
    kinds = kinds if kinds is not None else ["OFF"] * len(other_fs)
    max_freq = max_in_beam_freq_mhz(plan["telescope"], plan["sep_arcsec"])
    listed = {kind: [[_file_key(f["url"]), f["size_bytes"]] for f, k in zip(other_fs, kinds, strict=True) if k == kind] for kind in ("ON", "OFF")}
    desc = json.dumps({"version": SEARCH_VERSION, "source_id": source_id, "on": [_file_key(on_f["url"]), on_f["size_bytes"]],
                       "offs": listed["OFF"], **({"on_other": listed["ON"]} if listed["ON"] else {}),
                       "max_freq_mhz": max_freq, "settings": _settings(cfg)},
                      sort_keys=True)
    return hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24], desc, max_freq


def _work_key(desc: str) -> str:
    """Name of the folder that holds a band's saved windows: the band's description with CUT_VERSION for SEARCH_VERSION."""
    d = json.loads(desc)
    d["version"] = CUT_VERSION
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode("utf-8")).hexdigest()[:24]


def _cached_band(sig: str, cache: RawCache) -> dict[str, Any] | None:
    cached = cache.lookup(RADIO_QUERY, SEARCH_VERSION, sig)
    if cached is None or cached.offline_fixture:
        return None
    try:
        return json.loads(cache.read_bytes(cached).decode("utf-8"))
    except (FileNotFoundError, OSError, ValueError):
        return None


def _store_band(res: dict[str, Any], plan: dict[str, Any], sig: str, desc: str, cache: RawCache) -> dict[str, Any]:
    res.update({"telescope": plan["telescope"], "sep_arcsec": plan["sep_arcsec"], "comparison": plan["comparison"], "on_mjd": plan["mjd"]})
    cache.put(archive=RADIO_QUERY, release=SEARCH_VERSION, request_signature=sig, query_text=desc,
              payload=json.dumps(res, sort_keys=True).encode("utf-8"), row_count=int(res["n_windows_evaluated"]), offline_fixture=False)
    return res


def plan_searched(source_id: str, plan: dict[str, Any], cfg: PipelineConfig, cache: RawCache) -> bool:
    """Every band of the plan already has a stored search result (its scan files are then not needed)."""
    return all(_cached_band(_band_signature(source_id, plan, on_f, other_fs, cfg, kinds)[0], cache) is not None
               for on_f, other_fs, kinds in _band_sets(plan))


def _search_plan(source_id: str, plan: dict[str, Any], data_dir: Path, channel: Any, cfg: PipelineConfig, cache: RawCache,
                 progress: Callable[[int, int], None] | None) -> list[dict[str, Any]]:
    """Search results of one plan, one per ON file (frequency band). Answers are cached by files and settings."""
    results: list[dict[str, Any]] = []
    for on_f, other_fs, kinds in _band_sets(plan):
        sig, desc, max_freq = _band_signature(source_id, plan, on_f, other_fs, cfg, kinds)
        res = _cached_band(sig, cache)
        if res is None:
            on = ScanFile(local_path(data_dir, on_f["url"]))
            others: list[ScanFile] = []
            try:
                others = [ScanFile(local_path(data_dir, f["url"])) for f in other_fs]
                res = search_scan_set(source_id, on, [o for o, k in zip(others, kinds, strict=True) if k == "OFF"], channel, cfg,
                                      max_freq_mhz=max_freq, progress=progress, offset_s=plan_offset_s(plan),
                                      other_ons=[o for o, k in zip(others, kinds, strict=True) if k == "ON"])
            finally:
                on.close()
                for o in others:
                    o.close()
            res = _store_band(res, plan, sig, desc, cache)
        results.append(res)
    return results


def search_plan_one_file_at_a_time(
    source_id: str, plan: dict[str, Any], data_dir: Path, channel: Any, cfg: PipelineConfig, cache: RawCache, *,
    keep_files: bool = False, log: Callable[[str], None] = lambda msg: None,
) -> list[dict[str, Any]]:
    """Download, use and delete each scan of a plan in turn, so the disk never holds more than one scan file.

    The windows cut from each scan are kept in `<data_dir>/_work`, and a finished phase is not repeated, so an
    interrupted run resumes. The result is stored in the cache under the same key `enrich_sample_radio` looks up.
    """
    results: list[dict[str, Any]] = []
    for on_f, other_fs, kinds in _band_sets(plan):
        sig, desc, max_freq = _band_signature(source_id, plan, on_f, other_fs, cfg, kinds)
        res = _cached_band(sig, cache)
        if res is not None:
            results.append(res)
            continue
        t_start = time.monotonic()
        work = Path(data_dir) / "_work" / f"{source_id}_{_work_key(desc)}"
        work.mkdir(parents=True, exist_ok=True)

        def _phase(tag: str, f: dict[str, Any], run: Callable[[ScanFile, Path], dict[str, Any]]) -> dict[str, Any]:
            meta_path, store = work / f"{tag}.json", work / f"{tag}.f32"
            if meta_path.is_file() and store.is_file():
                return json.loads(meta_path.read_text(encoding="utf-8"))
            dest = local_path(data_dir, f["url"])
            last = [time.monotonic()]

            def _dl(have: int, size: int) -> None:
                if time.monotonic() - last[0] >= 60.0:
                    last[0] = time.monotonic()
                    log(f"{dest.name}: {have / 1e9:.1f} of {size / 1e9:.1f} GB")

            log(f"{tag}: {dest.name}: {download_file(f['url'], dest, f['size_bytes'], f['md5sum'], progress=_dl)}")
            scan = ScanFile(dest)
            try:
                rec = run(scan, store)
            finally:
                scan.close()
            meta_path.write_text(json.dumps(rec) + "\n", encoding="utf-8")
            if not keep_files:
                dest.unlink()
            log(f"{tag}: done after {time.monotonic() - t_start:.0f} s")
            return rec

        def _blocks(tag: str) -> Callable[[int, int], None]:
            return lambda d, t: log(f"{tag}: block {d}/{t}") if d == t or d % 32 == 0 else None

        offset_s = plan_offset_s(plan)
        found = _phase("on", on_f, lambda scan, store: find_candidate_windows(
            scan, cfg, store, max_freq_mhz=max_freq, progress=_blocks("on"), span_s=offset_s + scan.n_ints * scan.tsamp_s))
        # OFF scans keep the names they had before cadences existed (off0, off1, ...); other ON scans are on1, on2, ...
        n_kind = {"ON": 0, "OFF": 0}
        tags = []
        for kind in kinds:
            tags.append(f"off{n_kind['OFF']}" if kind == "OFF" else f"on{n_kind['ON'] + 1}")
            n_kind[kind] += 1
        recs = [(_phase(tag, f, lambda scan, store, _tag=tag: extract_windows(scan, found, store, progress=_blocks(_tag))), work / f"{tag}.f32")
                for tag, f in zip(tags, other_fs, strict=True)]
        res = evaluate_windows(source_id, found, work / "on.f32", recs, channel, cfg, kinds=kinds,
                               progress=lambda d, t: log(f"scoring: window {d}/{t}") if d == t or d % 2000 == 0 else None)
        res["elapsed_s"] = round(time.monotonic() - t_start, 1)
        results.append(_store_band(res, plan, sig, desc, cache))
        log(f"evaluated {res['n_windows_evaluated']} windows, {res['n_windows_flagged']} flagged")
    return results


def plan_work_dirs(source_id: str, plan: dict[str, Any], data_dir: Path, cfg: PipelineConfig) -> list[Path]:
    """The folders under `<data_dir>/_work` that hold the windows cut for a plan (one per frequency band)."""
    return [Path(data_dir) / "_work" / f"{source_id}_{_work_key(_band_signature(source_id, plan, on_f, other_fs, cfg, kinds)[1])}"
            for on_f, other_fs, kinds in _band_sets(plan)]


def _band_summary(res: dict[str, Any]) -> dict[str, Any]:
    keep = ("version", "on_file", "on_target", "off_targets", "freq_lo_mhz", "freq_hi_mhz", "max_in_beam_freq_mhz", "n_on_scans",
            "n_channels_searched", "n_windows_evaluated", "n_windows_flagged", "n_windows_flagged_outside_rfi_bands", "elapsed_s")
    return {k: res[k] for k in keep if k in res} | {"selected_hit": (res.get("selected") or {}).get("hit")}


def followup_star(
    source_id: str, ra: float, dec: float, cfg: PipelineConfig, cache: RawCache, channel: Any, data_dir: Path, *,
    all_scan_sets: bool = False, keep_windows: bool = False, min_free_bytes: int = 15_000_000_000, dry_run: bool = False,
    attempts: int = 3, retry_wait_s: float = 60.0, sleep: Callable[[float], Any] = time.sleep,
    log: Callable[[str], None] = lambda msg: None,
) -> dict[str, Any]:
    """The recordings of one star, searched on request (`oddstar radio followup`).

    The index is asked whether the star lies in a beam; if it does, the complete scan set with the fewest bytes
    (or every complete scan set) is fetched and searched one file at a time. A scan set that already has a stored
    result is not fetched again. Before a scan set is fetched the disk must have `min_free_bytes` free. The windows
    cut from the scans are deleted once the result is stored, unless `keep_windows`. With `dry_run` nothing is
    fetched: the record lists what would be. The stored results are the ones `search radio_narrowband --live` reads.

    Status: "no_coverage", "no_comparison", "searched", "would_search" (dry run), "disk_full" or "failed: <error>".
    """
    sid = str(source_id)
    rec: dict[str, Any] = {"source_id": sid, "ra_deg": ra, "dec_deg": dec, "scan_sets": []}
    try:
        rows, _, _ = bl.query_star_rows(ra, dec, cache, sleep=sleep)
        cov = bl.beam_coverage(ra, dec, rows)
        rec["index_status"] = cov["status"]
        if cov["status"] != "in_beam":
            return rec | {"status": "no_coverage"}
        plans = [p for p in bl.plan_star_scans(ra, dec, rows, cache, sleep=sleep) if p["complete"]]
    except ArchiveAccessError as exc:
        return rec | {"status": f"failed: {str(exc)[:300]}"}
    rec["n_complete_scan_sets"] = len(plans)
    if not plans:
        return rec | {"status": "no_comparison"}
    chosen = sorted(plans, key=lambda p: p["mjd"]) if all_scan_sets else [min(plans, key=lambda p: p["bytes_on_plus_offs"])]
    statuses: list[str] = []
    for plan in chosen:
        entry: dict[str, Any] = {"target": plan["target"], "telescope": plan["telescope"], "mjd": plan["mjd"], "sep_arcsec": plan["sep_arcsec"],
                                 "comparison": plan["comparison"], "bytes": plan["bytes_on_plus_offs"]}
        rec["scan_sets"].append(entry)
        stored = plan_searched(sid, plan, cfg, cache)
        if dry_run and not stored:
            entry["status"] = "would_search"
            statuses.append("would_search")
            continue
        if not stored:
            free = shutil.disk_usage(Path(data_dir) if Path(data_dir).exists() else Path(data_dir).parent).free
            if free < min_free_bytes:
                entry["status"] = "disk_full"
                entry["free_bytes"] = free
                statuses.append("disk_full")
                log(f"{sid}: {free / 1e9:.1f} GB free, {min_free_bytes / 1e9:.1f} GB wanted; scan set not fetched")
                continue
        error, retryable = "", True
        results: list[dict[str, Any]] | None = None
        for attempt in range(1, max(1, attempts) + 1):
            try:
                results = search_plan_one_file_at_a_time(sid, plan, data_dir, channel, cfg, cache, log=log)
                break
            except ArchiveAccessError as exc:
                error, retryable = str(exc)[:300], bool(exc.retryable)
                log(f"{sid}: attempt {attempt} failed: {error}")
                if not retryable or attempt >= max(1, attempts):
                    break
                sleep(retry_wait_s)
        if results is None:
            entry["status"] = f"failed: {error}"
            statuses.append("failed")
            # No scan file stays behind; a partial download stays only where a repeat could resume it
            on_files, off_sets = plan_files(plan)
            for f in [*on_files, *[x for fs in off_sets for x in fs]]:
                dest = local_path(data_dir, f["url"])
                dest.unlink(missing_ok=True)
                if not retryable:
                    dest.with_name(dest.name + ".part").unlink(missing_ok=True)
            continue
        entry.update({"status": "searched", "already_stored": stored, "bands": [_band_summary(r) for r in results]})
        statuses.append("searched")
        if not keep_windows:
            for work in plan_work_dirs(sid, plan, data_dir, cfg):
                if work.is_dir():
                    shutil.rmtree(work)
    if "searched" in statuses:
        rec["status"] = "searched"
    elif "would_search" in statuses:
        rec["status"] = "would_search"
    elif "disk_full" in statuses:
        rec["status"] = "disk_full"
    else:
        rec["status"] = next(e["status"] for e in rec["scan_sets"] if str(e["status"]).startswith("failed"))
    rec["n_windows_flagged"] = sum(b["n_windows_flagged"] for e in rec["scan_sets"] for b in e.get("bands", []))
    rec["n_windows_flagged_outside_rfi_bands"] = sum(b["n_windows_flagged_outside_rfi_bands"] for e in rec["scan_sets"] for b in e.get("bands", []))
    return rec


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
            rows, _, hit = bl.query_star_rows(ra, dec, cache, sleep=sleep)
            cov = bl.beam_coverage(ra, dec, rows)
            if cov["status"] != "in_beam":
                status = "no_coverage"
                rec["index_status"] = cov["status"]
            else:
                plans = [p for p in bl.plan_star_scans(ra, dec, rows, cache, sleep=sleep) if p["complete"]]
                rec["n_complete_scan_sets"] = len(plans)
                ready = [p for p in plans if plan_searched(sid, p, cfg, cache) or plan_on_disk(p, data_dir)]
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


HTTPS_HOSTS = ("bldata.berkeley.edu",)


def secure_url(url: str) -> str:
    """The https address of a file the index lists as http on a host known to serve both."""
    scheme, _, rest = url.partition("://")
    return f"https://{rest}" if scheme == "http" and rest.split("/", 1)[0] in HTTPS_HOSTS else url


def has_md5(md5sum: str) -> bool:
    """The index gives a real MD5 (32 hex digits, not all zero); it also writes "-" or zeros for none."""
    m = str(md5sum or "").strip().lower()
    return len(m) == 32 and all(c in "0123456789abcdef" for c in m) and set(m) != {"0"}


def download_file(url: str, dest: Path, size_bytes: int, md5sum: str, *, chunk_bytes: int = 8 << 20,
                  progress: Callable[[int, int], None] | None = None, timeout_s: float = 60.0, attempts: int = 6,
                  sleep: Callable[[float], Any] = time.sleep) -> str:
    """Fetch one archive file to `dest`, resuming a partial download, and check its size and MD5.

    A transfer of many gigabytes is often cut by the server or a proxy. A cut transfer is resumed from the bytes
    already on disk, up to `attempts` requests in a row that add nothing. Returns "present" (already complete) or
    "downloaded". The partial file is `dest` + ".part" and is kept on failure, so a later call resumes too.
    """
    import requests

    dest = Path(dest)
    if dest.is_file() and (not size_bytes or dest.stat().st_size == size_bytes):
        return "present"
    url = secure_url(url)
    part = dest.with_name(dest.name + ".part")
    dest.parent.mkdir(parents=True, exist_ok=True)
    if size_bytes and part.is_file() and part.stat().st_size > size_bytes:
        part.unlink()
    idle = 0
    while not (size_bytes and part.is_file() and part.stat().st_size == size_bytes):
        have = start = part.stat().st_size if part.is_file() else 0
        headers = {"Range": f"bytes={have}-"} if have else {}
        error = ""
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
            error = f"{type(exc).__name__}: {str(exc)[:160]}"
        if not size_bytes and not error:
            break  # no size to check against: the body the server sent is the file
        now = part.stat().st_size if part.is_file() else 0
        if size_bytes and now > size_bytes:
            raise ArchiveAccessError(f"{url}: got {now} bytes, the index lists {size_bytes}", retryable=False)
        idle = 0 if now > start else idle + 1
        if idle >= attempts or (not size_bytes and error):
            raise ArchiveAccessError(f"{url}: {error or 'the transfer ended early'}; {now} of {size_bytes} bytes on disk (rerun to resume)",
                                     retryable=True)
        if now != size_bytes:
            sleep(min(60.0, 5.0 * (idle + 1)))
    if has_md5(md5sum):
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
    ap.add_argument("--search", action="store_true",
                    help="download, search and delete one scan file at a time (needs disk for one file only) and print the result summary")
    ap.add_argument("--keep-files", action="store_true", help="with --search: keep each scan file after it was used")
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
    if args.search:
        from oddstar.channels import build_channel

        cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
        channel = build_channel("radio_narrowband", cfg, root)
        results = search_plan_one_file_at_a_time(str(args.source_id), plan, data_dir, channel, cfg, cache, keep_files=args.keep_files,
                                                 log=lambda msg: print(f"  {msg}", file=sys.stderr, flush=True))
        print(json.dumps([{k: v for k, v in r.items() if k not in ("selected", "windows")} | {
            "selected_hit": (r["selected"] or {}).get("hit"), "top_windows": r["windows"][:20]} for r in results], indent=2))
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
    "evaluate_windows", "extract_windows", "find_candidate_windows", "followup_star", "plan_work_dirs", "max_in_beam_freq_mhz", "plan_files", "plan_offset_s", "plan_on_disk", "plan_other_scans", "plan_searched", "cut_channels",
    "prefilter_block", "radio_search_status", "read_sigproc_header", "rfi_allocation", "search_coverage_summary", "search_plan_one_file_at_a_time",
    "search_scan_set",
]

if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
