"""Live radio search of Breakthrough Listen scans on local disk (offline; scans are small synthetic files)."""

from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oddstar.channels import build_channel
from oddstar.channels.radio_narrowband.channel import dedoppler_search_scan
from oddstar.cli import dispatch
from oddstar.config import load_config, repo_root
from oddstar.errors import ArchiveAccessError
from oddstar.io import bl_opendata as bl
from oddstar.io import bl_radio as br
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.vet.radio_checks import run_radio_vetting

DT = 18.253611008      # s per spectrum of the archive's fine-resolution product
DF_MHZ = 2.7939677e-6  # MHz per channel
N_T, N_F = 16, 8192
F_LOW = 1400.0         # MHz, centre of the lowest channel of the synthetic scans
T_ON = 57585.5083


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


@pytest.fixture(scope="module")
def small_cfg(cfg: Any) -> Any:
    """The real settings with a smaller block and window, so several blocks fit in a small file."""
    c = cfg.model_copy(deep=True)
    c.channels["radio_narrowband"].update({"live_block_channels": 2048, "live_window_channels": 1024})
    return c


def _noise(seed: int, n_f: int = N_F) -> np.ndarray:
    return np.random.default_rng(seed).normal(1000.0, 50.0, size=(N_T, n_f)).astype(np.float32)


def _tone(data: np.ndarray, channel: int, drift_hz_s: float, amp: float = 600.0) -> np.ndarray:
    """Add a drifting line (ascending channel index) to a (time, channel) array."""
    df_hz = DF_MHZ * 1e6
    for ti in range(data.shape[0]):
        ch = channel + int(round(drift_hz_s * ti * DT / df_hz))
        if 0 <= ch < data.shape[1]:
            data[ti, ch] += amp
    return data


def _key(name: str, kind: str, value: Any) -> bytes:
    raw = struct.pack("<i", len(name)) + name.encode("ascii")
    if kind == "i":
        return raw + struct.pack("<i", value)
    if kind == "d":
        return raw + struct.pack("<d", value)
    return raw + struct.pack("<i", len(value)) + value.encode("ascii")


def write_fil(path: Path, data: np.ndarray, *, f_low: float = F_LOW, tstart: float = T_ON, descending: bool = True, source: str = "HIPX") -> Path:
    """`data` is (time, channel) in ascending frequency; the archive stores descending (foff < 0)."""
    n_f = data.shape[1]
    fch1, foff = (f_low + (n_f - 1) * DF_MHZ, -DF_MHZ) if descending else (f_low, DF_MHZ)
    hdr = (struct.pack("<i", 12) + b"HEADER_START" + _key("source_name", "s", source) + _key("machine_id", "i", 20) + _key("telescope_id", "i", 6)
           + _key("data_type", "i", 1) + _key("fch1", "d", fch1) + _key("foff", "d", foff) + _key("nchans", "i", n_f) + _key("nbits", "i", 32)
           + _key("nifs", "i", 1) + _key("tstart", "d", tstart) + _key("tsamp", "d", DT) + _key("src_raj", "d", 3848.5) + _key("src_dej", "d", -132500.0)
           + struct.pack("<i", 10) + b"HEADER_END")
    stored = data[:, ::-1] if descending else data
    path.write_bytes(hdr + np.ascontiguousarray(stored, dtype="<f4").tobytes())
    return path


def write_h5(path: Path, data: np.ndarray, *, f_low: float = F_LOW, tstart: float = T_ON, source: str = "HIPX") -> Path:
    h5py = pytest.importorskip("h5py")
    n_f = data.shape[1]
    with h5py.File(path, "w") as f:
        ds = f.create_dataset("data", data=np.ascontiguousarray(data[:, ::-1], dtype="<f4")[:, None, :], chunks=(1, 1, 1024), compression="gzip")
        for k, v in {"fch1": f_low + (n_f - 1) * DF_MHZ, "foff": -DF_MHZ, "nchans": n_f, "tsamp": DT, "tstart": tstart, "nifs": 1, "nbits": 32}.items():
            ds.attrs[k] = v
        ds.attrs["source_name"] = np.bytes_(source)
    return path


def test_scan_files_read_in_ascending_frequency(tmp_path: Path) -> None:
    data = _noise(1, 4096)
    for i, path in enumerate([write_fil(tmp_path / "d.fil", data), write_fil(tmp_path / "a.fil", data, descending=False),
                              write_h5(tmp_path / "d.h5", data)]):
        sf = br.ScanFile(path)
        assert (sf.n_ints, sf.nchans) == (N_T, 4096) and sf.source_name == "HIPX" and sf.tsamp_s == pytest.approx(DT)
        assert sf.df_hz == pytest.approx(DF_MHZ * 1e6) and sf.f_low_mhz == pytest.approx(F_LOW, abs=1e-9), i
        np.testing.assert_array_equal(sf.read_ascending(100, 612), data[:, 100:612])
        assert sf.channel_of_freq(sf.freq_of_channel(777)) == 777
        with pytest.raises(ArchiveAccessError, match="outside"):
            sf.read_ascending(4000, 5000)
        sf.close()
    (tmp_path / "bad.fil").write_bytes(b"\x00" * 64)
    with pytest.raises(ArchiveAccessError, match="SIGPROC"):
        br.ScanFile(tmp_path / "bad.fil")
    with pytest.raises(ArchiveAccessError, match="not found"):
        br.ScanFile(tmp_path / "missing.fil")
    (tmp_path / "x.txt").write_text("x")
    with pytest.raises(ArchiveAccessError, match="unknown scan file type"):
        br.ScanFile(tmp_path / "x.txt")


def test_prefilter_uses_the_channels_statistic_and_drift_grid() -> None:
    window = 1024
    block = _tone(_noise(2, 4 * window), 2 * window + 300, 0.5)
    snr, drift, ch = br.prefilter_block(block, DT, DF_MHZ * 1e6, 4.0, window)
    assert int(np.argmax(snr)) == 2 and ch[2] == 2 * window + 300 and drift[2] == pytest.approx(0.5)
    c_snr, c_f, c_drift = dedoppler_search_scan(block[:, 2 * window:3 * window].astype(float), DT, DF_MHZ * 1e6, max_abs_drift=4.0)
    assert snr[2] == pytest.approx(c_snr, rel=1e-4) and c_drift == pytest.approx(0.5) and c_f == pytest.approx((300 - window / 2) * DF_MHZ * 1e6)
    assert np.all(np.delete(snr, 2) < 8.0), "noise alone does not reach the prefilter level"
    assert len(np.linspace(-4.0, 4.0, br.N_DRIFT_TRIALS)) == 161
    with pytest.raises(ValueError):
        br.prefilter_block(block[:, :100], DT, DF_MHZ * 1e6, 4.0, window)


def test_search_flags_a_line_present_only_in_the_on_scan(tmp_path: Path, cfg: Any, small_cfg: Any) -> None:
    channel = build_channel("radio_narrowband", cfg, repo_root())
    on_data = _noise(10)
    _tone(on_data, 700, 0.5)                 # toward the target only
    _tone(on_data, 3000, 0.0)                # steady interference, also in both OFF scans
    _tone(on_data, 5000, 0.5)                # drifting, but also in an OFF scan
    _tone(on_data, 7500, 0.0)                # ON only but not drifting
    off1 = _tone(_noise(11), 3000, 0.0)
    off2 = _tone(_tone(_noise(12), 3000, 0.0), 5000 + int(round(0.5 * 684.0 / (DF_MHZ * 1e6))), 0.5)
    on = br.ScanFile(write_fil(tmp_path / "on.fil", on_data))
    offs = [br.ScanFile(write_fil(tmp_path / "off1.fil", off1, tstart=T_ON - 342.0 / 86400.0, source="OFFA")),
            br.ScanFile(write_fil(tmp_path / "off2.fil", off2, tstart=T_ON + 342.0 / 86400.0, source="OFFB"))]
    seen: list[tuple[int, int]] = []
    res = br.search_scan_set("77", on, offs, channel, small_cfg, progress=lambda d, t: seen.append((d, t)))
    assert res["n_blocks"] == 4 and seen[-1] == (4, 4) and res["n_channels_searched"] == N_F and res["n_blocks_without_off"] == 0
    assert res["n_candidate_windows"] == 4 and res["n_windows_evaluated"] == 4 and res["n_windows_flagged"] == 1
    by_ch = {on.channel_of_freq(w["freq_mhz"]): w for w in res["windows"]}
    assert sorted(by_ch) == [700, 3000, 5000, 7500]
    assert by_ch[700]["flagged"] and by_ch[700]["off_scan_hits"] == 0 and by_ch[700]["drift_hz_s"] == pytest.approx(0.5) and by_ch[700]["min_on_snr"] > 20
    assert not by_ch[3000]["flagged"] and by_ch[3000]["off_scan_hits"] == 2
    assert not by_ch[5000]["flagged"] and by_ch[5000]["off_scan_hits"] == 1, "the OFF line is where the drift carries the ON line 342 s later"
    assert not by_ch[7500]["flagged"] and by_ch[7500]["off_scan_hits"] == 0 and abs(by_ch[7500]["drift_hz_s"]) < 0.05
    assert res["windows"][0] is by_ch[700] or res["windows"][0] == by_ch[700], "flagged windows are listed first"
    sel = res["selected"]
    assert sel["hit"]["flagged"] and on.channel_of_freq(sel["hit"]["freq_mhz"]) == 700
    cad = sel["cadence"]
    assert [s["kind"] for s in cad["scans"]] == ["OFF", "ON", "OFF"] and [s["t_start_s"] for s in cad["scans"]] == pytest.approx([0.0, 342.0, 684.0], abs=0.01)
    # The small test window (1024 channels = 2.9 kHz) cannot hold 4 Hz/s over the 634 s of this scan set: the cut is doubled
    assert res["window_channels"] == 2048 and br.cut_channels(1024, 4.0, 342.0 + N_T * DT, DF_MHZ * 1e6) == 2048
    assert np.asarray(cad["scans"][1]["power"]).shape == (N_T, 2048) and cad["dt_s"] == pytest.approx(DT)
    # The selected cadence is an ordinary cadence: the channel alone reproduces the flag from it
    m = channel.measure(json.loads(json.dumps(cad)), {})
    assert channel.score(m).flagged and m.values["off_scan_hits"] == 0.0

    # Above the frequency where the star leaves the beam nothing is searched
    cut = br.search_scan_set("77", on, offs, channel, small_cfg, max_freq_mhz=on.freq_of_channel(2047))
    assert cut["n_blocks"] == 1 and cut["n_windows_flagged"] == 1 and cut["n_channels_searched"] == 2048
    low = br.search_scan_set("77", on, offs, channel, small_cfg, max_freq_mhz=on.freq_of_channel(500))
    assert low["n_blocks"] == 0 and low["selected"] is None and low["n_windows_evaluated"] == 0
    for sf in (on, *offs):
        sf.close()


def test_a_quiet_scan_still_gets_a_measured_window_and_offs_are_matched_by_frequency(tmp_path: Path, cfg: Any, small_cfg: Any) -> None:
    channel = build_channel("radio_narrowband", cfg, repo_root())
    on = br.ScanFile(write_fil(tmp_path / "on.fil", _noise(20)))
    # This OFF scan starts 1024 channels lower in frequency and is wider: the same sky frequencies sit at other indices
    wide = _noise(21, N_F + 2048)
    off = br.ScanFile(write_fil(tmp_path / "off.fil", wide, f_low=F_LOW - 1024 * DF_MHZ, tstart=T_ON + 0.004))
    res = br.search_scan_set("5", on, [off], channel, small_cfg)
    assert res["n_candidate_windows"] == 0 and res["n_windows_evaluated"] == 1 and res["n_windows_flagged"] == 0
    hit = res["selected"]["hit"]
    assert not hit["flagged"] and 3.0 < hit["max_on_snr"] < 10.0
    lo = on.channel_of_freq(res["selected"]["cadence"]["center_freq_mhz"] - 1024 * DF_MHZ)
    np.testing.assert_allclose(np.asarray(res["selected"]["cadence"]["scans"][1]["power"]), wide[:, lo + 1024:lo + 3072], rtol=1e-6)

    # An OFF scan that does not hold the ON frequencies compares nothing, and says so
    short = br.ScanFile(write_fil(tmp_path / "short.fil", _noise(22, 4096), tstart=T_ON + 0.004))
    res2 = br.search_scan_set("5", on, [short], channel, small_cfg)
    assert res2["n_blocks_without_off"] == 2 and res2["n_channels_searched"] == 4096
    other = br.ScanFile(write_fil(tmp_path / "far.fil", _noise(23), f_low=1500.0, tstart=T_ON + 0.004))
    res3 = br.search_scan_set("5", on, [other], channel, small_cfg)
    assert res3["n_blocks_without_off"] == 4 and res3["n_channels_searched"] == 0 and res3["selected"] is None
    with pytest.raises(ArchiveAccessError, match="at least one OFF"):
        br.search_scan_set("5", on, [], channel, small_cfg)
    for sf in (on, off, short, other):
        sf.close()


def test_a_cadence_requires_the_line_in_every_on_scan(tmp_path: Path, cfg: Any) -> None:
    """A, B, A, C, A: a line has to be in all three ON scans; one that shows in a single ON scan is interference."""
    channel = build_channel("radio_narrowband", cfg, repo_root())
    c = cfg.model_copy(deep=True)
    c.channels["radio_narrowband"].update({"live_block_channels": 8192, "live_window_channels": 2048})
    df_hz, step = DF_MHZ * 1e6, 640.0

    def shifted(ch: int, drift: float, dt: float) -> int:
        return ch + int(round(drift * dt / df_hz))

    n_f = 16384
    ons = [_noise(30, n_f), _noise(31, n_f), _noise(32, n_f)]
    for k, data in enumerate(ons):                   # ON scans start at -640 s, 0 s, +640 s
        _tone(data, shifted(3000, 0.5, (k - 1) * step), 0.5)     # in every ON scan, where its drift carries it
        _tone(data, shifted(11000, 3.5, (k - 1) * step), 3.5)    # the same at a drift that leaves a 2048-channel window
    _tone(ons[1], 7000, 0.5)                         # in the middle ON scan only
    scans = [br.ScanFile(write_fil(tmp_path / f"on{k}.fil", data, tstart=T_ON + (k - 1) * step / 86400.0, source="NGCX")) for k, data in enumerate(ons)]
    offs = [br.ScanFile(write_fil(tmp_path / f"off{k}.fil", _noise(40 + k, n_f), tstart=T_ON + (k - 0.5) * step / 86400.0, source=f"OFF{k}"))
            for k in range(3)]
    assert br.cut_channels(2048, 4.0, step + N_T * DT, df_hz) == 4096 and br.cut_channels(2048, 4.0, 342.0 + N_T * DT, df_hz) == 2048
    assert br.cut_channels(2048, 4.0, None, df_hz) == 2048

    alone = br.search_scan_set("9", scans[1], offs[:2], channel, c)
    assert alone["n_windows_flagged"] == 3 and "on_other_files" not in alone and alone["window_channels"] == 2048
    res = br.search_scan_set("9", scans[1], offs, channel, c, other_ons=[scans[0], scans[2]])
    by_ch = {scans[1].channel_of_freq(w["freq_mhz"]): w for w in res["windows"]}
    assert res["n_on_scans"] == 3 and res["on_other_files"] == ["on0.fil", "on2.fil"] and res["off_files"] == ["off0.fil", "off1.fil", "off2.fil"]
    assert res["window_channels"] == 4096 and res["n_windows_evaluated"] == 3 and res["n_windows_flagged"] == 2
    assert by_ch[3000]["flagged"] and by_ch[3000]["on_persistence"] == 1.0 and by_ch[3000]["min_on_snr"] > 20
    assert by_ch[11000]["flagged"], "a fast line stays inside the wider cut"
    assert not by_ch[7000]["flagged"] and by_ch[7000]["min_on_snr"] < 10, "seen in one ON scan of three"
    cad = res["selected"]["cadence"]
    assert [sc["kind"] for sc in cad["scans"]] == ["ON", "OFF", "ON", "OFF", "ON", "OFF"]
    assert np.asarray(cad["scans"][0]["power"]).shape == (N_T, 4096)
    m = channel.measure(json.loads(json.dumps(cad)), {})
    assert channel.score(m).flagged and m.values["on_persistence"] == 1.0
    with pytest.raises(ArchiveAccessError, match="at least one OFF"):
        br.search_scan_set("9", scans[1], [], channel, c, other_ons=[scans[0]])
    for sf in (*scans, *offs):
        sf.close()

    # The plan of a cadence: files, kinds and the key of the stored result
    def entry(name: str, minutes: float) -> dict[str, Any]:
        return {"target": name, "mjd": T_ON + minutes / 1440.0, "files": [{"url": f"http://x/{name}.fil", "size_bytes": 5, "md5sum": "-", "center_freq_mhz": 1400.0}]}

    single = {"telescope": "GBT", "sep_arcsec": 10.0, "mjd": T_ON, "comparison": "before_after", "on_files": entry("on1", 0.0)["files"],
              "off_before": entry("offa", -5.3), "off_after": entry("offb", 5.3), "off_simultaneous": []}
    cadence = dict(single, comparison="cadence", off_before=None, off_after=None, on_other=[entry("on0", -10.6), entry("on2", 10.6)],
                   off_cadence=[entry("offa", -5.3), entry("offb", 5.3)])
    assert [k for k, _ in br.plan_other_scans(single)] == ["OFF", "OFF"] and br.plan_offset_s(single) == pytest.approx(318.0, abs=0.1)
    assert [(k, e["target"]) for k, e in br.plan_other_scans(cadence)] == [("ON", "on0"), ("OFF", "offa"), ("OFF", "offb"), ("ON", "on2")]
    assert br.plan_offset_s(cadence) == pytest.approx(636.0, abs=0.1) and len(br.plan_files(cadence)[1]) == 4
    (on_f, other_fs, kinds), = br._band_sets(cadence)
    assert kinds == ["ON", "OFF", "OFF", "ON"]
    sig_c, desc_c, _ = br._band_signature("9", cadence, on_f, other_fs, c, kinds)
    (on_s, other_s, kinds_s), = br._band_sets(single)
    sig_s, desc_s, _ = br._band_signature("9", single, on_s, other_s, c, kinds_s)
    assert sig_c != sig_s and "on_other" in desc_c and "on_other" not in desc_s
    assert br._band_signature("9", single, on_s, other_s, c)[0] == sig_s, "a plan without other ON scans keeps its key"
    # The saved windows are keyed by the cutting version: a new search version finds them again
    key = br._work_key(desc_s)
    assert key != sig_s and key != br._work_key(desc_c)
    old_version = br.SEARCH_VERSION
    try:
        br.SEARCH_VERSION = "bl_local_search_test"
        sig_new, desc_new, _ = br._band_signature("9", single, on_s, other_s, c, kinds_s)
        assert sig_new != sig_s and br._work_key(desc_new) == key
    finally:
        br.SEARCH_VERSION = old_version


def test_allocations_and_beam_edge() -> None:
    assert "GPS L1" in (br.rfi_allocation(1575.42) or "") and "Iridium" in (br.rfi_allocation(1621.0) or "")
    assert br.rfi_allocation(1420.4) is None and br.rfi_allocation(1700.0) is not None and br.rfi_allocation(6000.0) is None
    for freq, word in ((1340.0, "radar"), (1381.05, "GPS L3"), (1440.0, "telemetry"), (1627.5, "Inmarsat"), (1750.0, "AWS"), (1906.0, "PCS"),
                       (751.5, "700 MHz"), (868.9, "800 MHz"), (915.0, "industrial"), (953.6, "paging"),
                       (1922.9, "cordless"), (1960.0, "PCS"), (2250.0, "telemetry"), (2440.0, "wireless networks"), (2800.0, "radar"),
                       (3000.0, "ship radar"), (3200.0, "radiolocation"), (3500.0, "mobile broadband"), (3600.0, "citizens broadband")):
        assert word in (br.rfi_allocation(freq) or ""), freq
    assert br.SEARCH_VERSION == "bl_local_search_v4"
    # From 698 to 4200 MHz the table is complete: a frequency is in a listed band or in one of the three passive bands
    rows = sorted(br.RFI_ALLOCATIONS_MHZ)
    assert list(br.RFI_ALLOCATIONS_MHZ) == rows and all(lo < hi for lo, hi, _ in rows)
    assert all(a[1] <= b[0] for a, b in zip(rows, rows[1:], strict=False)), "bands do not overlap"
    lo_c, hi_c = br.ALLOCATIONS_COMPLETE_MHZ
    edges = sorted([(lo, hi) for lo, hi, _ in rows if lo_c <= lo < hi_c] + list(br.PASSIVE_BANDS_MHZ))
    assert edges[0][0] == lo_c and edges[-1][1] == hi_c and all(a[1] == b[0] for a, b in zip(edges, edges[1:], strict=False)), "no gap"
    for lo, hi in br.PASSIVE_BANDS_MHZ:
        assert br.rfi_allocation((lo + hi) / 2.0) is None
    assert br.rfi_allocation(1420.4057) is None and br.rfi_allocation(1665.4) is None and br.rfi_allocation(2695.0) is None
    assert br.rfi_allocation(611.0) is None and br.rfi_allocation(680.0) is None and br.rfi_allocation(5000.0) is None
    f_edge = br.max_in_beam_freq_mhz("GBT", 209.3)
    assert f_edge == pytest.approx(1773.0, abs=3.0)
    assert bl.half_power_radius_arcsec("GBT", f_edge) == pytest.approx(209.3, rel=1e-6)
    assert br.max_in_beam_freq_mhz("GBT", 0.0) == float("inf") and br.max_in_beam_freq_mhz("APF", 10.0) is None


# ------------------------------------------------------------------------------------------------ survey level

RA, DEC = 9.7, 13.4


def _index_row(name: str, target: str, ra: float, dec: float, mjd: float, size: int, freq: float = 1400.0114) -> dict[str, Any]:
    return {"target": target, "ra": ra, "decl": dec, "center_freq": freq, "telescope": "GBT", "size": size, "mjd": mjd,
            "file_type": "filterbank", "md5sum": "0" * 32, "url": f"http://blpd0.ssl.berkeley.edu/x/{name}"}


class _Archive:
    """A fake file index over synthetic scans in `data_dir` (position and time-window queries)."""

    def __init__(self, data_dir: Path) -> None:
        self.rows: list[dict[str, Any]] = []
        self.data_dir = data_dir
        self.calls = 0
        self.fail_ra: set[str] = set()

    def add(self, name: str, target: str, ra: float, dec: float, mjd: float, data: np.ndarray | None, size: int | None = None) -> None:
        if data is not None:
            write_fil(self.data_dir / name, data, tstart=mjd, source=target)
            size = (self.data_dir / name).stat().st_size
        self.rows.append(_index_row(name, target, ra, dec, mjd, int(size or 0)))

    def get(self, url: str, params: Any, **kw: Any) -> Any:
        self.calls += 1
        if "time-start" in params:
            hits = [r for r in self.rows if float(params["time-start"]) <= r["mjd"] <= float(params["time-end"])]
        else:
            if params["pos-ra"] in self.fail_ra:
                raise ArchiveAccessError("index unreachable", retryable=True)
            ra, dec, rad = float(params["pos-ra"]), float(params["pos-dec"]), float(params["pos-rad"])
            hits = [r for r in self.rows if ot.angular_separation_arcsec(ra, dec, r["ra"], r["decl"]) <= rad * 3600.0]

        class _R:
            status_code = 200
            content = json.dumps({"data": hits, "result": "success"}).encode("utf-8")

        return _R()


def _sample(tmp_path: Path) -> tuple[list[dict[str, Any]], _Archive, Path]:
    data_dir = tmp_path / "scans"
    data_dir.mkdir()
    arc = _Archive(data_dir)
    step = 5.7 / 1440.0
    # Star 1: 100 arcsec from the HIPON pointing, scans on disk; a drifting line toward the target only
    arc.add("on_0022.gpuspec.0000.fil", "HIPON", RA, DEC + 100.0 / 3600.0, T_ON, _tone(_noise(30), 4000, 0.3))
    arc.add("off_0021.gpuspec.0000.fil", "HIPOFF", RA + 3.0, DEC, T_ON - step, _noise(31))
    arc.add("off_0023.gpuspec.0000.fil", "HIPOFF", RA + 3.0, DEC, T_ON + step, _noise(32))
    # Star 2: inside a beam with a full scan set in the archive, nothing downloaded
    arc.add("b_on.gpuspec.0000.fil", "HIPB", 40.0, 10.0, T_ON + 1.0, None, size=15_000_000_000)
    arc.add("b_off1.gpuspec.0000.fil", "HIPB1", 43.0, 10.0, T_ON + 1.0 - step, None, size=15_000_000_000)
    arc.add("b_off2.gpuspec.0000.fil", "HIPB2", 43.0, 11.0, T_ON + 1.0 + step, None, size=15_000_000_000)
    # Star 3: inside a beam, but the archive holds no other scan near in time
    arc.add("c_on.gpuspec.0000.fil", "HIPC", 80.0, 20.0, T_ON + 2.0, None, size=15_000_000_000)
    # Star 4: a pointing 10 arcmin away; star 5: nothing; star 6: the index fails
    arc.add("d_on.gpuspec.0000.fil", "HIPD", 120.0, 30.0, T_ON + 3.0, None, size=15_000_000_000)
    stars = [{"source_id": str(i + 1), "ra_deg": ra, "dec_deg": dec, "live_enriched": False, "distance_pc": 300.0, "teff_k": 5200.0}
             for i, (ra, dec) in enumerate([(RA, DEC), (40.0, 10.0 + 60.0 / 3600.0), (80.0, 20.0 + 60.0 / 3600.0),
                                            (120.0, 30.0 + 600.0 / 3600.0), (200.0, 45.0), (250.0, -5.0)])]
    arc.fail_ra = {"250.00000"}
    return stars, arc, data_dir


def test_survey_searches_only_stars_with_a_complete_scan_set_on_disk(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Any) -> None:
    stars, arc, data_dir = _sample(tmp_path)
    monkeypatch.setattr(ot, "_get_with_retries", arc.get)
    cache = RawCache(tmp_path / "raw")
    channel = build_channel("radio_narrowband", cfg, repo_root())
    seen: list[str] = []
    out, summary = br.enrich_sample_radio(stars, cfg, cache, channel, data_dir, sleep=lambda s: None,
                                          progress=lambda d, t, sid, st: seen.append(f"{sid}:{st.split(':')[0]}"))
    assert seen == ["1:ok", "2:not_downloaded", "3:no_comparison", "4:no_coverage", "5:no_coverage", "6:failed"]
    assert [br.radio_search_status(s) for s in out] == ["searched", "not_downloaded", "no_comparison", "no_coverage", "no_coverage", "failed"]
    assert summary["n_ok"] == 1 and summary["n_not_downloaded"] == 1 and summary["n_no_comparison"] == 1 and summary["n_no_coverage"] == 2
    assert summary["n_failed"] == 1 and summary["failed_ids"] == ["6"] and summary["not_downloaded_ids"] == ["2"] and summary["n_windows_flagged"] == 1
    cov = br.search_coverage_summary(out, cfg)
    assert cov["n_input"] == 6 and cov["n_searched"] == 1 and cov["n_not_searched"] == 5
    assert cov["not_searched_ids"] == {"no_coverage": ["4", "5"], "not_downloaded": ["2"], "no_comparison": ["3"], "failed": ["6"]}

    s1 = out[0]
    assert s1["vetting_archive"]["query_status"]["breakthrough_listen"] == "ok" and s1["in_rfi_band"] is False
    rs = s1["radio_search"]
    assert rs["n_scan_sets_searched"] == 1 and rs["selected_hit"]["flagged"] and rs["scan_sets"][0]["comparison"] == "before_after"
    assert rs["scan_sets"][0]["on_file"] == "on_0022.gpuspec.0000.fil" and rs["scan_sets"][0]["max_in_beam_freq_mhz"] > 3000.0
    assert "selected" not in rs["scan_sets"][0]
    m = channel.measure(s1, {})
    assert channel.score(m).flagged and m.values["drift_hz_s"] == pytest.approx(0.3)
    vet = {r.check: r.status for r in run_radio_vetting(s1, cfg, tmp_path / "vet")}
    assert vet == {"off_scan": "pass", "zero_drift": "pass", "rfi_band": "pass", "multi_target": "not_evaluated", "reobservation": "not_evaluated"}

    assert out[1]["radio_fetch"]["smallest_scan_set_bytes"] == 45_000_000_000 and "scans" not in out[1]
    assert out[1]["vetting_archive"]["query_status"]["breakthrough_listen"] == "not_downloaded"
    assert {r.status for r in run_radio_vetting(out[1], cfg, tmp_path / "vet2")} == {"not_evaluated"}
    m2 = channel.measure(out[1], {})
    assert not channel.score(m2).flagged and m2.values["max_on_snr"] == 0.0

    # Second pass: the index answers and the search result come from the cache
    n_calls = arc.calls
    monkeypatch.setattr(br, "search_scan_set", lambda *a, **kw: (_ for _ in ()).throw(AssertionError("searched twice")))
    out2, summary2 = br.enrich_sample_radio(stars, cfg, cache, channel, data_dir, sleep=lambda s: None)
    assert arc.calls == n_calls + 1, "only the failed star is asked for again"
    assert out2[0]["radio_search"]["selected_hit"] == rs["selected_hit"] and summary2["n_ok"] == 1
    out3, summary3 = br.enrich_sample_radio(stars, cfg, cache, channel, data_dir, max_stars=1, sleep=lambda s: None)
    assert summary3["n_not_fetched"] == 5 and out3[1]["radio_fetch"] == {"status": "not_fetched"}

    # A scan file of the wrong size is not a downloaded file
    (data_dir / "off_0023.gpuspec.0000.fil").write_bytes(b"short")
    out4, _ = br.enrich_sample_radio(stars[:1], cfg, RawCache(tmp_path / "raw2"), channel, data_dir, sleep=lambda s: None)
    assert out4[0]["radio_fetch"]["status"] == "not_downloaded"


def test_download_resumes_and_checks_the_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import requests

    body = bytes(range(256)) * 40
    md5 = hashlib.md5(body).hexdigest()  # noqa: S324
    calls: list[dict[str, str]] = []

    class _Resp:
        def __init__(self, data: bytes, status: int) -> None:
            self.data, self.status_code = data, status

        def __enter__(self) -> "_Resp":
            return self

        def __exit__(self, *a: Any) -> None:
            return None

        def iter_content(self, chunk_size: int) -> Any:
            for i in range(0, len(self.data), 1000):
                yield self.data[i:i + 1000]

    def fake_get(url: str, headers: dict[str, str], stream: bool, timeout: float) -> _Resp:
        calls.append(dict(headers))
        if "Range" in headers:
            start = int(headers["Range"].split("=")[1].rstrip("-"))
            return _Resp(body[start:], 206)
        return _Resp(body, 200)

    monkeypatch.setattr(requests, "get", fake_get)
    dest = tmp_path / "d" / "scan.fil"
    assert br.download_file("http://x/scan.fil", dest, len(body), md5) == "downloaded" and dest.read_bytes() == body
    assert br.download_file("http://x/scan.fil", dest, len(body), md5) == "present" and len(calls) == 1

    dest2 = tmp_path / "d" / "scan2.fil"
    dest2.with_name("scan2.fil.part").write_bytes(body[:3000])
    assert br.download_file("http://x/scan2.fil", dest2, len(body), md5) == "downloaded" and dest2.read_bytes() == body
    assert calls[-1] == {"Range": "bytes=3000-"}

    dest3 = tmp_path / "d" / "scan3.fil"
    with pytest.raises(ArchiveAccessError, match="MD5"):
        br.download_file("http://x/scan3.fil", dest3, len(body), "f" * 32)

    # A transfer that is cut twice is resumed from the bytes on disk, without the caller asking again
    cuts = [4000, 7000]

    class _Cut(_Resp):
        def iter_content(self, chunk_size: int) -> Any:
            yield self.data[:3000]
            raise requests.exceptions.ChunkedEncodingError("Connection broken: reset by peer")

    def flaky_get(url: str, headers: dict[str, str], stream: bool, timeout: float) -> _Resp:
        calls.append(dict(headers))
        start = int(headers["Range"].split("=")[1].rstrip("-")) if "Range" in headers else 0
        return _Cut(body[start:], 206 if start else 200) if cuts and cuts.pop(0) else _Resp(body[start:], 206 if start else 200)

    monkeypatch.setattr(requests, "get", flaky_get)
    n, slept = len(calls), []
    dest4 = tmp_path / "d" / "scan4.fil"
    assert br.download_file("http://x/scan4.fil", dest4, len(body), md5, sleep=slept.append) == "downloaded" and dest4.read_bytes() == body
    assert calls[n:] == [{}, {"Range": "bytes=3000-"}, {"Range": "bytes=6000-"}] and len(slept) == 2

    # A server that keeps cutting without sending anything gives up after `attempts` requests and keeps the part file
    class _Dead(_Resp):
        def iter_content(self, chunk_size: int) -> Any:
            raise requests.exceptions.ChunkedEncodingError("Connection broken")
            yield b""

    monkeypatch.setattr(requests, "get", lambda url, headers, stream, timeout: _Dead(b"", 200))
    dest5 = tmp_path / "d" / "scan5.fil"
    with pytest.raises(ArchiveAccessError, match="rerun to resume"):
        br.download_file("http://x/scan5.fil", dest5, len(body), md5, attempts=3, sleep=slept.append)
    assert not dest5.exists()
    assert not dest3.exists() and dest3.with_name("scan3.fil.part").exists()
    monkeypatch.setattr(requests, "get", lambda url, headers, stream, timeout: _Resp(body[:500], 200))
    with pytest.raises(ArchiveAccessError, match="rerun to resume"):
        br.download_file("http://x/scan6.fil", tmp_path / "d" / "scan6.fil", len(body), md5, attempts=2, sleep=slept.append)
    monkeypatch.setattr(requests, "get", lambda url, headers, stream, timeout: _Resp(b"", 404))
    with pytest.raises(ArchiveAccessError, match="HTTP 404"):
        br.download_file("http://x/scan5.fil", tmp_path / "d" / "scan5.fil", len(body), md5)

    # The index writes "-" (or zeros) where it has no MD5: the size is then the only check
    assert not br.has_md5("-") and not br.has_md5("0" * 32) and not br.has_md5("") and br.has_md5(md5.upper())
    seen: list[str] = []
    monkeypatch.setattr(requests, "get", lambda url, headers, stream, timeout: (seen.append(url), _Resp(body, 200))[1])
    assert br.download_file("http://bldata.berkeley.edu/a/scan7.h5", tmp_path / "d" / "scan7.h5", len(body), "-") == "downloaded"
    assert br.download_file("http://blpd0.ssl.berkeley.edu/a/scan8.h5", tmp_path / "d" / "scan8.h5", len(body), md5) == "downloaded"
    # a host known to serve https is asked over https; the stored result does not depend on the scheme the index listed
    assert seen == ["https://bldata.berkeley.edu/a/scan7.h5", "http://blpd0.ssl.berkeley.edu/a/scan8.h5"]
    assert br._file_key("http://bldata.berkeley.edu/a/scan7.h5") == br._file_key("https://bldata.berkeley.edu/a/scan7.h5")


def test_live_radio_search_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """search -> vet -> null-report -> run record: only stars whose ON and OFF scans were compared count as searched."""
    import shutil

    from oddstar import cli
    from oddstar.run_records import build_run_record

    stars, arc, data_dir = _sample(tmp_path)
    monkeypatch.setattr(ot, "_get_with_retries", arc.get)
    enriched: list[str] = []
    monkeypatch.setattr(cli, "fetch_live_halo_bulge_sample", lambda cfg, cache, **kw: ([json.loads(json.dumps(s)) for s in stars], []))
    monkeypatch.setattr(cli, "enrich_star_live_for_vetting", lambda star, cfg, cache: (enriched.append(star["source_id"]) or dict(star), []))
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    monkeypatch.setattr("time.sleep", lambda s: None)
    run_dir = repo_root() / "runs" / "test-radio-survey"
    try:
        code, payload = dispatch(["search", "radio_narrowband", "--live", "--region", "halo", "--limit", "6", "--radio-data-dir", str(data_dir),
                                  "--run-id", "test-radio-survey", "--json"])
        assert code == 0, payload
        d = payload["data"]
        assert d["n_input"] == 6 and d["n_searched"] == 1 and d["flagged_ids"] == ["1"] and enriched == []
        assert d["coverage"]["n_not_downloaded"] == 1 and d["coverage"]["n_no_comparison"] == 1 and d["coverage"]["n_no_coverage"] == 2
        assert d["coverage"]["n_failed"] == 1 and d["radio_fetch"]["n_windows_flagged"] == 1
        assert any("5 of 6 stars not searched" in w and "not_downloaded=1" in w for w in payload["warnings"])
        saved = json.loads((run_dir / "search_radio_narrowband.json").read_text(encoding="utf-8"))
        summaries = {s["source_id"]: s["summary"] for s in saved["scores"]}
        assert summaries["2"].startswith("NOT SEARCHED (not_downloaded)") and summaries["5"].startswith("NOT SEARCHED (no_coverage)")
        assert summaries["1"].startswith("Radio min_ON_SNR=") and list(saved["radio_search"]) == ["1"]
        kept = json.loads((run_dir / "live_stars.json").read_text(encoding="utf-8"))
        assert list(kept) == ["1"] and len(kept["1"]["scans"]) == 3

        code, payload = dispatch(["vet", "test-radio-survey", "--json"])
        assert code == 0, payload
        vet = json.loads((run_dir / "vetting.json").read_text(encoding="utf-8"))["vetting"]["1"]
        status = {v["check"]: v["status"] for v in vet}
        assert status["off_scan"] == "pass" and status["zero_drift"] == "pass" and status["rfi_band"] == "pass"
        assert status["reobservation"] == "not_evaluated" and status["multi_target"] == "not_evaluated"

        code, payload = dispatch(["null-report", "test-radio-survey", "--json"])
        assert code == 0, payload
        rep = json.loads((run_dir / "null_report_radio_narrowband.json").read_text(encoding="utf-8"))
        assert rep["n_searched"] == 1 and rep["n_flagged"] == 1

        rec, samples = build_run_record(run_dir)
        s = rec["searches"][0]
        assert s["coverage"]["n_searched"] == 1 and s["coverage"]["n_not_downloaded"] == 1 and s["radio_fetch"]["n_ok"] == 1
        assert len(s["not_searched_ids"]) == 5 and samples[0]["n"] == 6
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def test_one_file_at_a_time_search_gives_the_same_result_and_leaves_no_scan_on_disk(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Any) -> None:
    import shutil

    stars, arc, archive_dir = _sample(tmp_path)
    monkeypatch.setattr(ot, "_get_with_retries", arc.get)
    channel = build_channel("radio_narrowband", cfg, repo_root())
    ref_cache = RawCache(tmp_path / "raw_ref")
    ref, _ = br.enrich_sample_radio(stars[:1], cfg, ref_cache, channel, archive_dir, sleep=lambda s: None)

    cache = RawCache(tmp_path / "raw")
    rows, _, _ = bl.query_bl_files(RA, DEC, cache)
    plan = bl.plan_star_scans(RA, DEC, rows, cache)[0]
    data_dir = tmp_path / "empty"
    data_dir.mkdir()
    fetched: list[str] = []
    on_disk_at_once: list[int] = []

    def fake_download(url: str, dest: Path, size_bytes: int, md5sum: str, **kw: Any) -> str:
        fetched.append(dest.name)
        shutil.copy(archive_dir / dest.name, dest)
        on_disk_at_once.append(len([p for p in data_dir.iterdir() if p.is_file()]))
        return "downloaded"

    monkeypatch.setattr(br, "download_file", fake_download)
    assert not br.plan_searched("1", plan, cfg, cache)
    log: list[str] = []
    res = br.search_plan_one_file_at_a_time("1", plan, data_dir, channel, cfg, cache, log=log.append)
    assert fetched == ["on_0022.gpuspec.0000.fil", "off_0021.gpuspec.0000.fil", "off_0023.gpuspec.0000.fil"] and on_disk_at_once == [1, 1, 1]
    assert [p.name for p in data_dir.iterdir() if p.is_file()] == [], "each scan is deleted after use"
    assert any("evaluated" in line for line in log)
    want = ref[0]["radio_search"]["scan_sets"][0]
    got = res[0]
    assert got["selected"]["hit"] == ref[0]["radio_search"]["selected_hit"]
    for k in ("n_windows_evaluated", "n_windows_flagged", "n_candidate_windows", "n_channels_searched", "windows", "comparison", "on_mjd"):
        assert got[k] == want[k], k

    # The stored result stands in for the files: the survey reports the star as searched with nothing on disk
    assert br.plan_searched("1", plan, cfg, cache)
    out, summary = br.enrich_sample_radio(stars[:2], cfg, cache, channel, data_dir, sleep=lambda s: None)
    assert [br.radio_search_status(s) for s in out] == ["searched", "not_downloaded"] and summary["n_ok"] == 1
    assert out[0]["scans"] == ref[0]["scans"]
    # A second call does nothing; an interrupted run repeats only the unfinished phase
    fetched.clear()
    assert br.search_plan_one_file_at_a_time("1", plan, data_dir, channel, cfg, cache) == res and fetched == []
    cache2 = RawCache(tmp_path / "raw2")
    work = next((data_dir / "_work").iterdir())
    (work / "off1.json").unlink()
    res2 = br.search_plan_one_file_at_a_time("1", plan, data_dir, channel, cfg, cache2)
    assert fetched == ["off_0023.gpuspec.0000.fil"] and res2[0]["selected"]["hit"] == got["selected"]["hit"]


def test_radio_followup_searches_chosen_stars_on_request(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Any) -> None:
    """`oddstar radio followup`: the recordings of chosen stars are fetched, searched and removed; the survey then reads the result."""
    import shutil

    from oddstar import cli

    stars, arc, archive_dir = _sample(tmp_path)
    monkeypatch.setattr(ot, "_get_with_retries", arc.get)
    channel = build_channel("radio_narrowband", cfg, repo_root())
    cache = RawCache(tmp_path / "raw")
    data_dir = tmp_path / "pass_through"
    data_dir.mkdir()
    fetched: list[str] = []

    def fake_download(url: str, dest: Path, size_bytes: int, md5sum: str, **kw: Any) -> str:
        if not (archive_dir / dest.name).is_file():
            raise ArchiveAccessError(f"{url}: HTTP 404", retryable=False)
        fetched.append(dest.name)
        shutil.copy(archive_dir / dest.name, dest)
        return "downloaded"

    monkeypatch.setattr(br, "download_file", fake_download)
    by_id = {s["source_id"]: s for s in stars}

    def follow(sid: str, **kw: Any) -> dict[str, Any]:
        return br.followup_star(sid, by_id[sid]["ra_deg"], by_id[sid]["dec_deg"], cfg, cache, channel, data_dir, sleep=lambda s: None, **kw)

    # A dry run lists the scan set and fetches nothing
    dry = follow("1", dry_run=True)
    assert dry["status"] == "would_search" and dry["scan_sets"][0]["status"] == "would_search" and dry["scan_sets"][0]["bytes"] > 0 and fetched == []
    # Too little free disk: nothing is fetched
    assert follow("1", min_free_bytes=10**18)["status"] == "disk_full" and fetched == []
    # The search: three files pass through, nothing stays on disk, the flagged window is reported
    rec = follow("1")
    assert rec["status"] == "searched" and len(fetched) == 3 and rec["n_windows_flagged"] == 1 and rec["n_windows_flagged_outside_rfi_bands"] == 1
    band = rec["scan_sets"][0]["bands"][0]
    assert band["selected_hit"]["in_rfi_band"] is False and band["n_windows_evaluated"] >= 1 and rec["scan_sets"][0]["already_stored"] is False
    assert [p for p in data_dir.rglob("*") if p.is_file()] == [], "neither a scan nor its cut windows stay on disk"
    # A repeat reads the stored result
    fetched.clear()
    again = follow("1")
    assert again["status"] == "searched" and again["scan_sets"][0]["already_stored"] is True and fetched == []
    # The other cases: recordings that cannot be fetched, no comparison scan, no pointing
    failed = follow("2")
    assert failed["status"].startswith("failed: ") and "HTTP 404" in failed["status"] and [p for p in data_dir.rglob("*") if p.is_file()] == []
    assert follow("3")["status"] == "no_comparison" and follow("5")["status"] == "no_coverage"
    assert follow("6")["status"].startswith("failed: ")
    # The survey reads the stored result of star 1 and downloads nothing
    out, summary = br.enrich_sample_radio(stars[:2], cfg, cache, channel, data_dir, sleep=lambda s: None)
    assert [br.radio_search_status(s) for s in out] == ["searched", "not_downloaded"] and summary["n_ok"] == 1

    # The command: by default the stars vetting left at Tier C or better; --ids and --all-covered choose otherwise
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    monkeypatch.setattr(cli, "fetch_gaia_astrometry", lambda ids, cfg, cache, **kw: {})  # Gaia holds none of the test ids
    runs = repo_root() / "runs"
    src, out_dir = runs / "test-followup-src", runs / "radio-followup-test-followup-src"
    try:
        src.mkdir(parents=True)
        score = {"source_id": "2", "channel": "ir_excess", "raw_statistic": 9.0, "false_alarm_prob": 1e-4, "flagged": True, "rank_score": 9.0,
                     "summary": "x"}
        (src / "search_ir_excess.json").write_text(json.dumps({
            "channel": "ir_excess", "sample_stars": stars, "flagged_ids": ["2"], "scores": [score]}), encoding="utf-8")
        code, payload = dispatch(["radio", "followup", "--from-run", "test-followup-src", "--radio-data-dir", str(data_dir), "--dry-run", "--json"])
        assert code == 0, payload
        d = payload["data"]
        assert d["selection"] == "tier C or better" and [s["source_id"] for s in d["stars"]] == ["2"] and d["stars"][0]["tier_before"] == "C"
        assert d["status_counts"] == {"would_search": 1} and not out_dir.exists(), "a dry run writes nothing"
        code, payload = dispatch(["radio", "followup", "--from-run", "test-followup-src", "--ids", "1,5", "--radio-data-dir", str(data_dir), "--json"])
        assert code == 0, payload
        assert payload["data"]["status_counts"] == {"searched": 1, "no_coverage": 1}
        saved = json.loads((out_dir / "radio_followup.json").read_text(encoding="utf-8"))
        assert [s["status"] for s in saved["stars"]] == ["searched", "no_coverage"] and saved["n_windows_flagged_outside_rfi_bands"] == 1
        code, payload = dispatch(["radio", "followup", "--from-run", "test-followup-src", "--all-covered", "--radio-data-dir", str(data_dir),
                                  "--dry-run", "--json"])
        assert code == 0 and payload["data"]["n_selected"] == 6 and payload["data"]["status_counts"]["searched"] == 1
        code, payload = dispatch(["radio", "followup", "--from-run", "test-followup-src", "--ids", "99", "--json"])
        assert code != 0 and "No position" in json.dumps(payload)
    finally:
        shutil.rmtree(src, ignore_errors=True)
        shutil.rmtree(out_dir, ignore_errors=True)


def test_two_nodes_of_one_parkes_scan_do_not_share_a_local_file(tmp_path: Path) -> None:
    # The two recording nodes of one Parkes multibeam scan write the same file name into different directories
    a = "http://blpd12.ssl.berkeley.edu/collate_mb_2/blc11/guppi_58326_62282_852088_G299.60-41.91_0001.0000.h5"
    b = "http://blpd12.ssl.berkeley.edu/collate_mb_2/blc12/guppi_58326_62282_852088_G299.60-41.91_0001.0000.h5"
    pa, pb = br.local_path(tmp_path, a), br.local_path(tmp_path, b)
    assert pa != pb and pa.name.endswith("_guppi_58326_62282_852088_G299.60-41.91_0001.0000.h5") and pa.suffix == ".h5"
    assert br.local_path(tmp_path, a.replace("http://", "https://")) == pa, "the scheme does not change the name"
    # Names that carry their node keep their name
    gbt = "https://bldata.berkeley.edu/blc07_guppi_59426_77412_NGC4168_0012.rawspec.0000.h5"
    assert br.local_path(tmp_path, gbt).name == "blc07_guppi_59426_77412_NGC4168_0012.rawspec.0000.h5"


# ---------------------------------------------------------------------------------------------------- the beam test


def test_range_reader_serves_reads_from_blocks_and_counts_what_it_fetched() -> None:
    raw = bytes(range(256)) * 40  # 10240 bytes
    asked: list[tuple[int, int]] = []

    def fetch(a: int, b: int) -> bytes:
        asked.append((a, b))
        return raw[a:b]

    r = br.RangeReader(fetch, len(raw), block_bytes=1000, max_blocks=8)
    assert r.readable() and r.seekable() and not r.writable() and r.tell() == 0
    assert r.read(10) == raw[:10] and asked == [(0, 1000)] and r.tell() == 10
    assert r.read(500) == raw[10:510] and asked == [(0, 1000)], "the same block is not fetched twice"
    r.seek(2500)
    assert r.read(3000) == raw[2500:5500] and asked[-1] == (2000, 6000), "neighbouring missing blocks are one request"
    r.seek(900)
    assert r.read(1200) == raw[900:2100] and asked[-1] == (1000, 2000), "only the block in the middle was missing"
    assert r.seek(-240, 2) == 10000 and r.read() == raw[10000:] and asked[-1] == (10000, 10240), "the last block is short"
    assert r.read(5) == b"" and r.seek(0, 2) == len(raw)
    r.seek(4, 0)
    r.seek(2, 1)
    buf = bytearray(8)
    assert r.readinto(buf) == 8 and bytes(buf) == raw[6:14]
    assert r.n_requests == len(asked) and r.bytes_fetched == sum(b - a for a, b in asked)
    with pytest.raises(ArchiveAccessError):
        br.RangeReader(lambda a, b: raw[a:b - 1], len(raw), block_bytes=1000).read(10)
    with pytest.raises(ArchiveAccessError):
        br.RangeReader(fetch, 0)


def test_http_range_fetch_accepts_only_the_asked_range() -> None:
    body = bytes(range(200))

    class _R:
        def __init__(self, status: int, content: bytes) -> None:
            self.status_code, self.content = status, content

        def close(self) -> None:
            pass

    class _Session:
        def __init__(self, answers: list[Any]) -> None:
            self.answers, self.ranges = answers, []

        def get(self, url: str, headers: dict[str, str], **kw: Any) -> Any:
            self.ranges.append(headers["Range"])
            a = self.answers.pop(0)
            if isinstance(a, Exception):
                raise a
            return a

    import requests

    s = _Session([_R(206, body[10:20])])
    assert br.http_range_fetch("http://x/f.h5", session=s)(10, 20) == body[10:20] and s.ranges == ["bytes=10-19"]
    # a server error, a cut transfer and a short body are tried again
    s = _Session([_R(503, b""), requests.exceptions.ConnectionError("cut"), _R(206, body[:5]), _R(206, body[:10])])
    assert br.http_range_fetch("http://x/f.h5", session=s, sleep=lambda t: None)(0, 10) == body[:10] and len(s.ranges) == 4
    # a server that ignores the range would send the whole file: refused, not retried
    s = _Session([_R(200, body)])
    with pytest.raises(ArchiveAccessError, match="does not serve byte ranges") as e:
        br.http_range_fetch("http://x/f.h5", session=s, sleep=lambda t: None)(0, 10)
    assert e.value.retryable is False and len(s.ranges) == 1
    s = _Session([_R(404, b"")])
    with pytest.raises(ArchiveAccessError, match="HTTP 404"):
        br.http_range_fetch("http://x/f.h5", session=s, sleep=lambda t: None)(0, 10)
    s = _Session([_R(500, b"")] * 3)
    with pytest.raises(ArchiveAccessError, match="HTTP 500") as e:
        br.http_range_fetch("http://x/f.h5", session=s, attempts=3, sleep=lambda t: None)(0, 10)
    assert e.value.retryable is True and len(s.ranges) == 3


def _as_remote(path: Path, asked: list[tuple[int, int]] | None = None, block: int = 4096) -> br.ScanFile:
    """A local scan file read the way an archive file is: through a `RangeReader`."""
    raw = path.read_bytes()

    def fetch(a: int, b: int) -> bytes:
        if asked is not None:
            asked.append((a, b))
        return raw[a:b]

    return br.ScanFile(path.name, fileobj=br.RangeReader(fetch, len(raw), block_bytes=block))


def test_a_scan_read_through_byte_ranges_equals_the_local_file_and_fetches_only_its_chunks(tmp_path: Path) -> None:
    pytest.importorskip("h5py")
    data = _noise(70)
    path = write_h5(tmp_path / "a.gpuspec.0000.h5", data)
    local = br.ScanFile(path)
    asked: list[tuple[int, int]] = []
    remote = _as_remote(path, asked, block=512)
    assert remote.nchans == N_F and remote.n_ints == N_T and remote.chunk_channels == 1024 and local.chunk_channels == 1024
    assert remote.freq_of_channel(0) == local.freq_of_channel(0) and remote.tstart_mjd == local.tstart_mjd
    n_open = sum(b - a for a, b in asked)
    # two windows inside one chunk of the file (ascending 3072..4095), one in another: read together, sliced apart
    wins = remote.read_windows([3100, 3600, 7000, -5, N_F - 100], 200)
    assert wins[3] is None and wins[4] is None, "a window the file does not hold is None"
    for a, w in zip([3100, 3600, 7000], wins[:3], strict=True):
        assert w is not None and np.array_equal(w, data[:, a:a + 200]) and np.array_equal(w, local.read_ascending(a, a + 200))
    n_read = sum(b - a for a, b in asked) - n_open
    assert 0 < n_read < 0.45 * path.stat().st_size, "two of eight chunk columns were fetched, not the file"
    assert remote.fileobj.bytes_fetched == sum(b - a for a, b in asked)
    remote.close()
    local.close()
    with pytest.raises(ArchiveAccessError, match="only an HDF5 scan"):
        br.ScanFile("x.gpuspec.0000.fil", fileobj=br.RangeReader(lambda a, b: b"", 10))


MB_RA, MB_DEC = 150.0, 20.0
KB_RA, KB_DEC = 160.0, 25.0


def _multibeam_sample(tmp_path: Path) -> tuple[list[dict[str, Any]], _Archive, Path]:
    """Two stars in a multibeam scan each: every beam of a scan carries the ON beam's target name and start time.

    Star 11 has two lines in its ON beam: line A (channel 1000) only there, line B (channel 6000, stronger) also in a
    near and in a far beam. Star 12 has only line B, also in a far beam; one of its far beam files is not in the archive.
    """
    pytest.importorskip("h5py")
    data_dir = tmp_path / "scans"
    data_dir.mkdir()
    arc = _Archive(data_dir)
    step = 5.7 / 1440.0
    asec = 1.0 / 3600.0

    def add(name: str, target: str, ra: float, dec: float, mjd: float, data: np.ndarray | None) -> None:
        size = 1000
        if data is not None:
            size = write_h5(data_dir / name, data, tstart=mjd, source=target).stat().st_size
        arc.rows.append(_index_row(name, target, ra, dec, mjd, size))

    def line_b(seed: int) -> np.ndarray:
        return _tone(_noise(seed), 6000, 0.2, amp=900.0)

    add("mb1_on.gpuspec.0000.h5", "MB1", MB_RA, MB_DEC + 100 * asec, T_ON, _tone(line_b(40), 1000, 0.3))
    add("mb1_near.gpuspec.0000.h5", "MB1", MB_RA, MB_DEC + 1000 * asec, T_ON, line_b(41))       # 1000": outside the beam, nearer than 6 radii
    add("mb1_far1.gpuspec.0000.h5", "MB1", MB_RA, MB_DEC + 2000 * asec, T_ON, line_b(42))
    add("mb1_far2.gpuspec.0000.h5", "MB1", MB_RA, MB_DEC - 2500 * asec, T_ON, _noise(43))
    add("mb1_far3.gpuspec.0000.h5", "MB1", MB_RA + 1.0, MB_DEC, T_ON, _noise(44))
    add("mb1_off1.gpuspec.0000.h5", "MBOFF", MB_RA + 3.0, MB_DEC, T_ON - step, _noise(45))
    add("mb1_off2.gpuspec.0000.h5", "MBOFF", MB_RA + 3.0, MB_DEC, T_ON + step, _noise(46))

    add("mb2_on.gpuspec.0000.h5", "MB2", KB_RA, KB_DEC + 100 * asec, T_ON + 1.0, line_b(50))
    add("mb2_far1.gpuspec.0000.h5", "MB2", KB_RA, KB_DEC + 2000 * asec, T_ON + 1.0, line_b(51))
    add("mb2_far2.gpuspec.0000.h5", "MB2", KB_RA, KB_DEC - 2500 * asec, T_ON + 1.0, None)       # listed by the index, not served
    add("mb2_far3.gpuspec.0000.h5", "MB2", KB_RA + 1.0, KB_DEC, T_ON + 1.0, _noise(53))
    add("mb2_off1.gpuspec.0000.h5", "MBOFF2", KB_RA + 3.0, KB_DEC, T_ON + 1.0 - step, _noise(54))
    add("mb2_off2.gpuspec.0000.h5", "MBOFF2", KB_RA + 3.0, KB_DEC, T_ON + 1.0 + step, _noise(55))
    stars = [{"source_id": sid, "ra_deg": ra, "dec_deg": dec, "live_enriched": False, "distance_pc": 300.0, "teff_k": 5200.0}
             for sid, ra, dec in (("11", MB_RA, MB_DEC), ("12", KB_RA, KB_DEC))]
    return stars, arc, data_dir


def test_beam_test_reads_the_other_beams_of_the_instant_and_feeds_the_multi_target_check(
        tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Any) -> None:
    import shutil

    from oddstar import cli

    stars, arc, data_dir = _multibeam_sample(tmp_path)
    monkeypatch.setattr(ot, "_get_with_retries", arc.get)
    cache = RawCache(tmp_path / "raw")
    channel = build_channel("radio_narrowband", cfg, repo_root())
    df_hz = DF_MHZ * 1e6

    # The search as it was: the beams carry the ON beam's name, so the comparison is the scans before and after
    out, summary = br.enrich_sample_radio(stars, cfg, cache, channel, data_dir, sleep=lambda s: None)
    assert [br.radio_search_status(s) for s in out] == ["searched", "searched"] and summary["n_windows_flagged_outside_rfi_bands"] == 3
    assert summary["n_windows_tested"] == 0
    s11, s12 = out
    assert s11["radio_search"]["scan_sets"][0]["comparison"] == "before_after" and "multi_target_match" not in s11
    f_a, f_b = F_LOW + 1000 * DF_MHZ, F_LOW + 6000 * DF_MHZ
    assert s11["radio_search"]["selected_hit"]["freq_mhz"] == pytest.approx(f_b, abs=2e-6), "the stronger line is the selected window"
    vet = {r.check: r for r in run_radio_vetting(s11, cfg, tmp_path / "vet0")}
    assert vet["multi_target"].status == "not_evaluated" and "oddstar radio beams" in vet["multi_target"].note

    plan = bl.plan_star_scans(MB_RA, MB_DEC, bl.query_star_rows(MB_RA, MB_DEC, cache, sleep=lambda s: None)[0], cache, sleep=lambda s: None)[0]
    beams = plan["beams_same_instant"]
    assert [b["files"][0]["url"].rsplit("/", 1)[-1][:8] for b in beams] == ["mb1_near", "mb1_far1", "mb1_far2", "mb1_far3"], "nearest first"
    assert [b["far"] for b in beams] == [False, True, True, True] and plan["off_simultaneous"] == []

    asked: list[tuple[int, int]] = []
    opened: list[str] = []

    def remote(f: dict[str, Any]) -> br.ScanFile:
        path = data_dir / f["url"].rsplit("/", 1)[-1]
        opened.append(path.name)
        if not path.is_file():
            raise ArchiveAccessError(f"{f['url']}: HTTP 404", retryable=False)
        return _as_remote(path, asked)

    logs: list[str] = []
    rec = br.beam_test_star("11", MB_RA, MB_DEC, cfg, cache, channel, open_scan=remote, sleep=lambda s: None, log=logs.append)
    assert rec["status"] == "tested" and rec["n_windows_tested"] == 2 and rec["n_windows_in_far_beams"] == 1 and rec["n_windows_in_no_far_beam"] == 1
    band = rec["scan_sets"][0]["bands"][0]
    assert band["already_stored"] is False and rec["scan_sets"][0]["n_beams_same_instant"] == 4
    test = band["test"]
    assert test["n_beams_read"] == 4 and test["n_far_beams_read"] == 3 and test["bytes_fetched"] == sum(b - a for a, b in asked) > 0
    assert test["n_range_requests"] == len(asked) and test["chunk_channels"] == 1024
    # the ON beam, the four other beams, then the two OFF scans for the window that may become the selected one
    assert opened == ["mb1_on.gpuspec.0000.h5", "mb1_near.gpuspec.0000.h5", "mb1_far1.gpuspec.0000.h5", "mb1_far2.gpuspec.0000.h5",
                      "mb1_far3.gpuspec.0000.h5", "mb1_off1.gpuspec.0000.h5", "mb1_off2.gpuspec.0000.h5"]
    w_b, w_a = sorted(band["windows"], key=lambda w: -w["freq_mhz"])
    assert w_a["freq_mhz"] == pytest.approx(f_a, abs=2e-6) and w_b["freq_mhz"] == pytest.approx(f_b, abs=2e-6)
    assert w_a["matches_stored_hit"] and w_b["matches_stored_hit"], "the ON beam read again gives the stored hit"
    assert w_b["multi_target_match"] is True and (w_b["n_far_beams"], w_b["n_far_beams_with_hit"], w_b["n_near_beams"], w_b["n_near_beams_with_hit"]) == (3, 1, 1, 1)
    assert w_a["multi_target_match"] is False and (w_a["n_far_beams"], w_a["n_far_beams_with_hit"], w_a["n_near_beams_with_hit"]) == (3, 0, 0)
    far1 = next(b for b in w_b["beams"] if b["beam"] == 1)
    assert far1["hit"] and far1["snr"] > 30 and far1["snr_at_on_freq"] == pytest.approx(far1["snr"], abs=0.01) and far1["drift_hz_s"] == pytest.approx(0.2)
    assert abs(far1["freq_offset_hz"] - w_b["on"]["freq_offset_hz"]) <= df_hz
    quiet = next(b for b in w_b["beams"] if b["beam"] == 2)
    assert not quiet["hit"] and quiet["snr"] < 10 and quiet["snr_at_on_freq"] < 6
    assert w_a["remeasured"]["flagged"] and w_a["remeasured"]["off_scan_hits"] == 0 and "remeasured" not in w_b
    assert any("line in 1 of 2 window(s)" in m for m in logs)

    # A repeat reads the stored test; --redo would read the beams again
    opened.clear()
    again = br.beam_test_star("11", MB_RA, MB_DEC, cfg, cache, channel, open_scan=remote, sleep=lambda s: None)
    assert again["scan_sets"][0]["bands"][0]["already_stored"] is True and opened == [] and again["n_windows_in_far_beams"] == 1

    # Star 12: its only line is in a far beam; a beam file the archive does not serve is listed with its error
    rec12 = br.beam_test_star("12", KB_RA, KB_DEC, cfg, cache, channel, open_scan=remote, sleep=lambda s: None)
    t12 = rec12["scan_sets"][0]["bands"][0]
    assert rec12["status"] == "tested" and rec12["n_windows_tested"] == 1 and rec12["n_windows_in_far_beams"] == 1
    assert t12["test"]["n_beams_listed"] == 3 and t12["test"]["n_beams_read"] == 2 and "HTTP 404" in t12["test"]["beams"][1]["error"]
    assert t12["windows"][0]["n_far_beams"] == 2 and t12["windows"][0]["n_far_beams_with_hit"] == 1

    # The survey now reads the stored tests: star 11's selected window is the line no far beam shows, star 12 keeps its line
    out, summary = br.enrich_sample_radio(stars, cfg, cache, channel, data_dir, sleep=lambda s: None)
    s11, s12 = out
    assert summary["n_windows_flagged_outside_rfi_bands"] == 3, "the search counts are not changed by the test"
    assert (summary["n_windows_tested"], summary["n_windows_in_far_beams"], summary["n_windows_in_no_far_beam"]) == (3, 2, 1)
    hit = s11["radio_search"]["selected_hit"]
    assert hit["freq_mhz"] == pytest.approx(f_a, abs=2e-6) and hit["multi_target_match"] is False and s11["multi_target_match"] is False
    assert s11["radio_search"]["scan_sets"][0]["beam_test"]["n_windows_in_far_beams"] == 1
    listed = s11["radio_search"]["scan_sets"][0]["windows"]
    assert [w.get("multi_target_match") for w in listed[:2]] == [False, True], "the window a far beam shows ranks below the one none shows"
    m = channel.measure(s11, {})
    assert channel.score(m).flagged and m.values["drift_hz_s"] == pytest.approx(0.3) and len(s11["scans"]) == 3
    vet = {r.check: r for r in run_radio_vetting(s11, cfg, tmp_path / "vet1")}
    assert vet["multi_target"].status == "pass" and "none of the 3 far beams" in vet["multi_target"].note
    assert vet["off_scan"].status == "pass" and vet["rfi_band"].status == "pass" and vet["reobservation"].status == "not_evaluated"
    assert s12["multi_target_match"] is True and s12["radio_search"]["selected_hit"]["freq_mhz"] == pytest.approx(f_b, abs=2e-6)
    assert channel.score(channel.measure(s12, {})).flagged, "still a flagged window: vetting rejects it"
    vet = {r.check: r for r in run_radio_vetting(s12, cfg, tmp_path / "vet2")}
    assert vet["multi_target"].status == "fail" and "in 1 of 2 far beams" in vet["multi_target"].note

    # The command
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    monkeypatch.setattr(br, "open_remote_scan", remote)
    runs = repo_root() / "runs"
    src, out_dir = runs / "test-beams-src", runs / "radio-beams-test-beams-src"
    try:
        src.mkdir(parents=True)
        score = {"source_id": "12", "channel": "radio_narrowband", "raw_statistic": 9.0, "false_alarm_prob": 1e-4, "flagged": True,
                 "rank_score": 9.0, "summary": "x"}
        (src / "search_radio_narrowband.json").write_text(json.dumps({
            "channel": "radio_narrowband", "sample_stars": stars, "flagged_ids": ["12"], "scores": [score]}), encoding="utf-8")
        code, payload = dispatch(["radio", "beams", "--from-run", "test-beams-src", "--json"])
        assert code == 0, payload
        d = payload["data"]
        assert d["selection"] == "tier C or better" and [s["source_id"] for s in d["stars"]] == ["12"] and d["status_counts"] == {"tested": 1}
        assert d["n_windows_tested"] == 1 and d["n_windows_in_far_beams"] == 1 and "1 also in a far beam" in payload["summary"]
        opened.clear()
        code, payload = dispatch(["radio", "beams", "--from-run", "test-beams-src", "--ids", "11,12", "--redo", "--json"])
        assert code == 0, payload
        assert payload["data"]["n_windows_tested"] == 3 and payload["data"]["n_windows_in_no_far_beam"] == 1 and len(opened) == 11
        saved = json.loads((out_dir / "radio_beams.json").read_text(encoding="utf-8"))
        assert [s["status"] for s in saved["stars"]] == ["tested", "tested"] and "cadence" not in json.dumps(saved)
    finally:
        shutil.rmtree(src, ignore_errors=True)
        shutil.rmtree(out_dir, ignore_errors=True)


def test_beam_test_has_nothing_to_do_without_other_beams_or_flagged_windows(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cfg: Any) -> None:
    stars, arc, data_dir = _sample(tmp_path)
    monkeypatch.setattr(ot, "_get_with_retries", arc.get)
    cache = RawCache(tmp_path / "raw")
    channel = build_channel("radio_narrowband", cfg, repo_root())

    def never(f: dict[str, Any]) -> br.ScanFile:
        raise AssertionError("no file is opened")

    def test(i: int) -> dict[str, Any]:
        return br.beam_test_star(stars[i]["source_id"], stars[i]["ra_deg"], stars[i]["dec_deg"], cfg, cache, channel, open_scan=never, sleep=lambda s: None)

    assert test(0)["status"] == "not_searched" and test(4)["status"] == "no_coverage" and test(5)["status"].startswith("failed: ")
    br.enrich_sample_radio(stars[:1], cfg, cache, channel, data_dir, sleep=lambda s: None)
    rec = test(0)
    assert rec["status"] == "no_other_beams" and rec["scan_sets"][0]["bands"][0]["n_windows_flagged_outside_rfi_bands"] == 1
    assert rec["n_windows_tested"] == 0
