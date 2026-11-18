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
    assert np.asarray(cad["scans"][1]["power"]).shape == (N_T, 1024) and cad["dt_s"] == pytest.approx(DT)
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
    lo = on.channel_of_freq(res["selected"]["cadence"]["center_freq_mhz"] - 512 * DF_MHZ)
    np.testing.assert_allclose(np.asarray(res["selected"]["cadence"]["scans"][1]["power"]), wide[:, lo + 1024:lo + 2048], rtol=1e-6)

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


def test_allocations_and_beam_edge() -> None:
    assert "GPS L1" in (br.rfi_allocation(1575.42) or "") and "Iridium" in (br.rfi_allocation(1621.0) or "")
    assert br.rfi_allocation(1420.4) is None and br.rfi_allocation(1700.0) is not None and br.rfi_allocation(6000.0) is None
    # search version 2: radar, GPS L3, air telemetry, satellite and telephone uplinks; the protected bands stay clear
    for freq, word in ((1340.0, "radar"), (1381.05, "GPS L3"), (1440.0, "telemetry"), (1627.5, "Inmarsat"), (1750.0, "AWS"), (1906.0, "PCS")):
        assert word in (br.rfi_allocation(freq) or ""), freq
    assert br.rfi_allocation(1410.0) is None and br.rfi_allocation(1665.0) is None and br.rfi_allocation(1800.0) is None
    assert br.SEARCH_VERSION == "bl_local_search_v2"
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
