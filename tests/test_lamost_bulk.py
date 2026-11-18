"""Survey-scale LAMOST spectra for spectral_lines (offline, fake HTTP) and honest not-searched accounting."""

from __future__ import annotations

import gzip
import io
import json
import re
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oddstar.channels import build_channel
from oddstar.channels.spectral_lines.channel import search_coverage_summary, spectra_search_status
from oddstar.cli import CHANNELS_WITHOUT_LIVE_SOURCE, dispatch
from oddstar.config import load_config, repo_root
from oddstar.errors import ArchiveAccessError
from oddstar.io import lamost_bulk as lb
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.vet.spectral_checks import KNOWN_SKY_LINES_A, LINE_EXCLUSION_LSF_FWHM, run_spectral_vetting

C_KMS = 299_792.458


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


def _air_to_vacuum(air: np.ndarray) -> np.ndarray:
    return air * (np.asarray(air) / lb.vacuum_to_air(air))  # first-order inverse, good to 1e-4 A here


def _lamost_fits(obsid: int, *, z: float = 0.0, helio_rv: float = 0.0, snrr: float = 80.0, lmjd: int = 57000, seed: int = 1,
                 emission: tuple[float, float] | None = None, masked: tuple[int, ...] = (), noise: float = 0.005,
                 gz: bool = True) -> bytes:
    """A LAMOST-like low-resolution file: one row of arrays on a log vacuum grid, observed (heliocentric) frame.

    The star has a continuum falling to the red, photospheric H-alpha at redshift `z`, and optionally an emission
    line `(rest air wavelength, amplitude)` at the line-spread width.
    """
    from astropy.io import fits

    rng = np.random.default_rng(seed)
    wave_vac = 10 ** np.arange(np.log10(6350.0), np.log10(6800.0), 1e-4)
    rest_air = lb.vacuum_to_air(wave_vac) / (1.0 + z)
    lsf_sig = lb.HALPHA_AIR_A / lb.LRS_RESOLVING_POWER / 2.355
    flux = 1.0 - 0.3 * np.exp(-0.5 * ((rest_air - lb.HALPHA_AIR_A) / 2.4) ** 2)
    if emission is not None:
        flux = flux + emission[1] * np.exp(-0.5 * ((rest_air - emission[0]) / lsf_sig) ** 2)
    flux = 3000.0 * (1.0 - 4e-4 * (wave_vac - 6575.0)) * flux * (1.0 + rng.normal(0.0, noise, len(wave_vac)))
    ormask = np.zeros(len(wave_vac), dtype=np.float32)
    for i in masked:
        ormask[i] = 16.0
        flux[i] = 0.0
    n = len(wave_vac)
    cols = [fits.Column(name=name, format=f"{n}E", array=np.asarray([arr], dtype=np.float32))
            for name, arr in (("FLUX", flux), ("IVAR", np.full(n, 1e-3)), ("WAVELENGTH", wave_vac), ("ANDMASK", np.zeros(n)),
                              ("ORMASK", ormask), ("NORMALIZATION", flux / 3000.0))]
    hdr = fits.Header()
    for k, v in (("OBSID", obsid), ("Z", z), ("HELIO_RV", helio_rv), ("SNRR", snrr), ("LMJD", lmjd), ("DATE-OBS", "2015-01-01T12:00:00"),
                 ("NEXP", 3), ("VACUUM", True)):
        hdr[k] = v
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(header=hdr), fits.BinTableHDU.from_columns(cols, name="COADD")]).writeto(buf)
    return gzip.compress(buf.getvalue()) if gz else buf.getvalue()


def _votable(rows: list[dict[str, Any]], status: str = "OK", columns: tuple[str, ...] = lb.CATALOGUE_COLUMNS) -> bytes:
    fields = "".join(f'<FIELD name="{c}" datatype="char" arraysize="*" />' for c in columns)
    trs = "".join("<TR>" + "".join(f"<TD>{r.get(c, '')}</TD>" for c in columns) + "</TR>" for r in rows)
    return (f'<?xml version="1.0" encoding="UTF-8"?><VOTABLE xmlns="http://www.ivoa.net/xml/VOTable/v1.3" version="1.4"><RESOURCE type="results">'
            f'<INFO name="QUERY_STATUS" value="{status}"><![CDATA[ no such table ]]></INFO><TABLE>{fields}<DATA><TABLEDATA>{trs}</TABLEDATA></DATA>'
            "</TABLE></RESOURCE></VOTABLE>").encode("utf-8")


def _row(obsid: int, gaia_id: str, *, snrr: float = 80.0, lmjd: int = 57000, z: float = 0.0) -> dict[str, Any]:
    return {"obsid": obsid, "gaia_source_id": gaia_id, "ra": 150.0, "dec": 30.0, "snrr": snrr, "lmjd": lmjd, "obsdate": "2015-01-01",
            "class": "STAR", "subclass": "G2", "z": z, "fibermask": 0}


class _Resp:
    def __init__(self, content: bytes) -> None:
        self.status_code = 200
        self.content = content


class _Lamost:
    """Fake LAMOST: catalogue rows per Gaia id and spectrum files per obsid. Counts every request."""

    def __init__(self, *, canary: bool = True) -> None:
        self.rows: list[dict[str, Any]] = [_row(14004088, lb.CANARY_GAIA_ID)] if canary else []
        self.files: dict[int, bytes | Exception] = {}
        self.catalogue_requests: list[list[str]] = []
        self.file_requests: list[int] = []
        self.catalogue_error: Exception | None = None

    def add(self, gaia_id: str, epochs: list[tuple[dict[str, Any], bytes | Exception]]) -> None:
        for row, body in epochs:
            self.rows.append(dict(row, gaia_source_id=gaia_id))
            self.files[int(row["obsid"])] = body

    def install(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def get(url: str, params: dict[str, Any] | None, **kw: Any) -> _Resp:
            if url == lb.LAMOST_TAP_URL:
                assert params is not None
                asked = re.findall(r"'(\d+)'", params["QUERY"])
                self.catalogue_requests.append(asked)
                if self.catalogue_error is not None:
                    raise self.catalogue_error
                return _Resp(_votable([r for r in self.rows if r["gaia_source_id"] in asked]))
            obsid = int(url.rsplit("/", 1)[-1])
            self.file_requests.append(obsid)
            body = self.files[obsid]
            if isinstance(body, Exception):
                raise body
            return _Resp(body)

        monkeypatch.setattr(ot, "_get_with_retries", get)


def _star(sid: str) -> dict[str, Any]:
    return {"source_id": sid, "ra_deg": 150.0, "dec_deg": 30.0, "live_enriched": False, "spectra": {}, "vetting_archive": {"query_status": {}}}


def _survey() -> tuple[list[dict[str, Any]], _Lamost]:
    """Six stars: 1 quiet two-epoch, 2 emission line in both epochs, 3 one epoch, 4 no spectrum, 5 fetch fails, 6 same-night pair."""
    b = _Lamost()
    b.add("1", [(_row(101, "1", lmjd=57000), _lamost_fits(101, seed=1)), (_row(102, "1", lmjd=57400, snrr=60.0), _lamost_fits(102, seed=2))])
    b.add("2", [(_row(201, "2", lmjd=57000, z=-0.0008), _lamost_fits(201, z=-0.0008, seed=3, emission=(6626.0, 0.15))),
                (_row(202, "2", lmjd=57900, z=-0.0008), _lamost_fits(202, z=-0.0008, seed=4, emission=(6626.0, 0.15)))])
    b.add("3", [(_row(301, "3"), _lamost_fits(301, seed=5))])
    b.add("5", [(_row(501, "5", lmjd=57000), _lamost_fits(501, seed=6)),
                (_row(502, "5", lmjd=57300), ArchiveAccessError("GET failed after 3 attempt(s): HTTP 503", retryable=True))])
    b.add("6", [(_row(601, "6", lmjd=57000), _lamost_fits(601, seed=7)), (_row(602, "6", lmjd=57000), _lamost_fits(602, seed=8)),
                (_row(603, "6", lmjd=57500, snrr=8.0), _lamost_fits(603, seed=9))])
    return [_star(str(i)) for i in range(1, 7)], b


# ----------------------------------------------------------------------------------------------- catalogue answers


def test_votable_parser_reads_rows_and_refuses_anything_else() -> None:
    rows = lb.parse_lamost_votable(_votable([_row(7, "42", snrr=33.5)]))
    assert rows == [{c: str(v) for c, v in _row(7, "42", snrr=33.5).items()}]
    assert lb.parse_lamost_votable(_votable([])) == []
    with pytest.raises(ArchiveAccessError, match="did not answer OK"):
        lb.parse_lamost_votable(_votable([], status="ERROR"))
    with pytest.raises(ArchiveAccessError, match="lacks column"):
        lb.parse_lamost_votable(_votable([], columns=("obsid", "ra")))
    with pytest.raises(ArchiveAccessError, match="not XML"):
        lb.parse_lamost_votable(b"<html>Service unavailable")


def test_an_empty_catalogue_answer_is_believed_only_with_the_known_star(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cache = RawCache(tmp_path / "raw")
    b = _Lamost(canary=False)  # the service answers OK with no rows at all, as a broken table would
    b.install(monkeypatch)
    epochs, errors, _ = lb.lookup_lamost_epochs(["1", "2"], cache, sleep=lambda s: None)
    assert epochs == {} and set(errors) == {"1", "2"} and "not trusted" in errors["1"]
    # The untrusted answer was not cached: a working service is asked again and now answers
    good = _Lamost()
    good.add("1", [(_row(101, "1"), b"")])
    good.install(monkeypatch)
    epochs, errors, manifests = lb.lookup_lamost_epochs(["1", "2"], cache, sleep=lambda s: None)
    assert errors == {} and [e["obsid"] for e in epochs["1"]] == [101] and epochs["2"] == [] and len(manifests) == 1
    assert lb.CANARY_GAIA_ID in good.catalogue_requests[0] and lb.CANARY_GAIA_ID not in epochs
    # ... and that answer is served from the cache afterwards
    lb.lookup_lamost_epochs(["1", "2"], cache, sleep=lambda s: None)
    assert len(good.catalogue_requests) == 1
    with pytest.raises(ValueError, match="not a Gaia source id"):
        lb.lookup_lamost_epochs(["1' OR 1=1 --"], cache)


# ------------------------------------------------------------------------------------------------ one spectrum


def test_reduce_moves_the_window_to_the_rest_frame_and_flattens_the_continuum() -> None:
    z = -0.0009  # -270 km/s: H-alpha is observed 5.9 A to the blue
    rec = lb.reduce_lamost_fits(_lamost_fits(11, z=z, helio_rv=12.0, snrr=91.0, lmjd=57123, masked=(150, 151)))
    assert rec["usable"] and rec["obsid"] == 11 and rec["z"] == pytest.approx(z) and rec["lmjd"] == 57123 and rec["n_bad_pixels"] <= 2
    w, f = np.asarray(rec["wavelength_a"]), np.asarray(rec["flux"])
    assert w[0] <= lb.GRID_START_A and w[-1] >= lb.GRID_STOP_A
    assert abs(w[int(np.argmin(f))] - lb.HALPHA_AIR_A) < 1.6  # within one pixel of rest H-alpha
    far = np.abs(w - lb.HALPHA_AIR_A) > 15.0
    assert abs(float(np.median(f[far])) - 1.0) < 0.005
    slope = float(np.polyfit(w[far], f[far], 1)[0])
    assert abs(slope) * 150.0 < 0.01  # the 6 per cent tilt of the input is gone
    assert np.all(np.isfinite(f)) and float(f.min()) > 0.5  # the masked pixels were interpolated, not left at zero
    # plain (not gzipped) files are read too
    assert lb.reduce_lamost_fits(_lamost_fits(11, gz=False))["usable"]


def test_reduce_refuses_a_spectrum_without_a_redshift_or_with_many_masked_pixels() -> None:
    rec = lb.reduce_lamost_fits(_lamost_fits(12, z=-9999.0))
    assert not rec["usable"] and "redshift" in rec["reason"] and rec["flux"] == []
    rec = lb.reduce_lamost_fits(_lamost_fits(13, masked=tuple(range(140, 175))))
    assert not rec["usable"] and "masked" in rec["reason"]


def test_sky_lines_follow_each_epoch_into_the_rest_frame() -> None:
    win = {"z": -0.001, "helio_rv_kms": 20.0}
    got = lb.sky_lines_in_rest_frame([win, {"z": 0.0, "helio_rv_kms": 0.0}])
    lab = [*KNOWN_SKY_LINES_A, *lb.LAMOST_EXTRA_SKY_LINES_A]
    shifted = [round(line * (1 + 20.0 / C_KMS) / (1 - 0.001), 3) for line in lab]
    assert got == sorted(set(shifted + [round(x, 3) for x in lab])) and len(lab) == 9
    assert all(6.5 < a - b < 7.5 for a, b in zip(shifted, lab))  # 320 km/s is about 7 A here


# ------------------------------------------------------------------------------------------------ the survey fetch


def test_survey_fetch_statuses_and_resume(cfg: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, b = _survey()
    b.install(monkeypatch)
    cache = RawCache(tmp_path / "raw")
    seen: list[tuple[str, str]] = []
    out, summary = lb.enrich_sample_lamost(stars, cfg, cache, sleep=lambda s: None, progress=lambda d, t, sid, st: seen.append((sid, st.split(":")[0])))
    by = {s["source_id"]: s for s in out}
    assert [s["source_id"] for s in out] == ["1", "2", "3", "4", "5", "6"] and stars[0]["spectra"] == {}  # input not modified
    assert {k: v["spectra_fetch"]["status"].split(":")[0] for k, v in by.items()} == {
        "1": "ok", "2": "ok", "3": "single_epoch", "4": "no_coverage", "5": "failed", "6": "single_epoch"}
    assert summary["n_ok"] == 2 and summary["n_single_epoch"] == 2 and summary["n_no_coverage"] == 1 and summary["n_failed"] == 1
    assert summary["failed_ids"] == ["5"] and summary["n_requested"] == 6 and summary["n_not_fetched"] == 0
    assert seen[-1] == ("6", "single_epoch") and len(seen) == 6

    sp = by["1"]["spectra"]
    grid = np.asarray(sp["wavelength_a"])
    assert grid[0] == lb.GRID_START_A and grid[-1] == lb.GRID_STOP_A and len(sp["flux_exp1"]) == len(sp["flux_exp2"]) == len(grid) == 101
    assert sp["lsf_fwhm_a"] == pytest.approx(3.646) and "sky_flux" not in sp and "detector_artifact" not in sp
    assert by["1"]["spectra_fetch"]["obsids_used"] == [101, 102] and len(by["1"]["spectra_fetch"]["manifest_ids"]) == 2
    assert by["1"]["vetting_archive"]["query_status"] == {"lamost_lrs": "ok", "spectra": "ok"}
    # A star that was not loaded has no spectra and no "spectra" status, so vetting cannot read stale arrays
    for sid in ("3", "4", "5", "6"):
        assert by[sid]["spectra"] == {} and "spectra" not in by[sid]["vetting_archive"]["query_status"]
    assert by["5"]["vetting_archive"]["query_status"]["lamost_lrs"].startswith("failed: ArchiveAccessError")
    # Star 6: the second epoch of the same night and the low-S/N epoch are not used, and are not downloaded
    notes = {e["obsid"]: e["note"] for e in by["6"]["spectra_fetch"]["epochs"]}
    assert notes[602] == "same night as a used epoch" and "S/N below" in notes[603] and 602 not in b.file_requests and 603 not in b.file_requests
    assert b.file_requests.count(502) == cfg.channels["spectral_lines"]["lamost_max_attempts"]  # retried, then reported

    # Second run: everything that succeeded comes from the cache; only the failed file is asked for again
    n_cat, b.file_requests[:] = len(b.catalogue_requests), []
    b.files[502] = _lamost_fits(502, seed=10)
    out2, summary2 = lb.enrich_sample_lamost(stars, cfg, cache, sleep=lambda s: None)
    assert len(b.catalogue_requests) == n_cat and b.file_requests == [502]
    assert summary2["n_ok"] == 3 and summary2["n_failed"] == 0
    assert out2[0]["spectra"] == by["1"]["spectra"]

    out3, summary3 = lb.enrich_sample_lamost(stars, cfg, cache, max_stars=2, sleep=lambda s: None)
    assert summary3["n_requested"] == 2 and summary3["n_not_fetched"] == 4 and out3[4]["spectra_fetch"] == {"status": "not_fetched"}


def test_a_failed_catalogue_lookup_is_a_failure_not_an_absence(cfg: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, b = _survey()
    b.catalogue_error = ArchiveAccessError("GET failed after 3 attempt(s): ConnectionError", retryable=True)
    b.install(monkeypatch)
    out, summary = lb.enrich_sample_lamost(stars, cfg, RawCache(tmp_path / "raw"), sleep=lambda s: None)
    assert summary["n_failed"] == 6 and summary["n_no_coverage"] == 0 and b.file_requests == []
    assert all(s["spectra_fetch"]["status"].startswith("failed: catalogue lookup") for s in out)
    assert search_coverage_summary(out, cfg)["n_failed"] == 6


# ------------------------------------------------------------------------------------------------ search and vetting


def test_channel_flags_a_line_seen_in_both_epochs_and_counts_only_loaded_stars(cfg: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, b = _survey()
    b.install(monkeypatch)
    out, _ = lb.enrich_sample_lamost(stars, cfg, RawCache(tmp_path / "raw"), sleep=lambda s: None)
    ch = build_channel("spectral_lines", cfg, repo_root())
    res = {}
    for s in out:
        m = ch.measure(s, {})
        res[s["source_id"]] = (m, ch.score(m))
    assert [sid for sid, (_, sc) in res.items() if sc.flagged] == ["2"]
    m2 = res["2"][0].values
    assert abs(m2["peak_wavelength_a"] - 6626.0) <= 1.5 and m2["repeat_ratio"] > 0.6
    assert res["1"][0].extra["search_status"] == "searched" and res["1"][1].summary.startswith("Spectral residual")
    assert res["3"][1].summary.startswith("NOT SEARCHED (single_epoch)") and res["4"][1].summary.startswith("NOT SEARCHED (no_coverage)")
    assert [spectra_search_status(s) for s in out] == ["searched", "searched", "single_epoch", "no_coverage", "failed", "single_epoch"]
    cov = search_coverage_summary(out, cfg)
    assert cov["n_input"] == 6 and cov["n_searched"] == 2 and cov["n_not_searched"] == 4 and cov["n_single_epoch"] == 2
    assert cov["not_searched_ids"] == {"single_epoch": ["3", "6"], "no_coverage": ["4"], "failed": ["5"]}
    assert spectra_search_status({"source_id": "9", "spectra": {}}) == "no_spectrum"

    # The same line in one epoch only is not a detection
    one = json.loads(json.dumps(out[1]))
    one["spectra"]["flux_exp2"] = out[0]["spectra"]["flux_exp2"]
    assert not ch.score(ch.measure(one, {})).flagged
    # A line at a listed night-sky position (displaced into this star's rest frame) is excluded
    sky = json.loads(json.dumps(out[1]))
    sky["spectra"]["sky_line_wavelengths_a"] = [6625.0]
    m_sky = ch.measure(sky, {})
    assert not ch.score(m_sky).flagged and abs(m_sky.values["peak_wavelength_a"] - 6625.0) > LINE_EXCLUSION_LSF_FWHM * 3.6
    # With the full night-sky list about half of the window is open to the search, and the measurement says so
    assert 0.3 < res["2"][0].extra["searched_window_fraction"] < 0.6 < m_sky.extra["searched_window_fraction"]
    # The peak's observed wavelength in each epoch is kept, so a sky or detector feature shared by many stars shows
    assert res["2"][0].extra["peak_observed_wavelengths_a"] == [round(m2["peak_wavelength_a"] * (1 - 0.0008), 2)] * 2

    # Vetting reads the loaded arrays: the checks with inputs are evaluated, the two without are named as missing
    vet = {r.check: r for r in run_spectral_vetting(out[1], cfg, tmp_path / "vet")}
    assert vet["cosmic_ray"].status == "pass" and vet["repeat"].status == "pass" and vet["stellar_emission"].status == "pass"
    assert vet["sky_line"].status == "not_evaluated" and "sky_flux" in vet["sky_line"].note
    assert vet["detector"].status == "not_evaluated"
    lines = json.loads(Path(vet["sky_line"].evidence_path).read_text(encoding="utf-8"))["metrics"] if getattr(vet["sky_line"], "evidence_path", None) else None
    if lines is not None:
        assert lines["sky_line_list_source"] == "spectra.sky_line_wavelengths_a"
    # ... and a star whose spectra were not loaded is not evaluated at all
    assert {r.status for r in run_spectral_vetting(out[2], cfg, tmp_path / "vet3")} == {"not_evaluated"}


def test_live_search_refuses_only_channels_without_a_real_data_source() -> None:
    assert set(CHANNELS_WITHOUT_LIVE_SOURCE) == {"radio_narrowband"}


def test_live_spectral_search_records_coverage_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """search -> vet -> null-report -> run record on a faked live sample: only two-epoch stars count as searched."""
    import shutil

    from oddstar import cli
    from oddstar.run_records import build_run_record

    stars, b = _survey()
    for s in stars:
        s.update({"distance_pc": 150.0, "teff_k": 5500.0})
    b.install(monkeypatch)
    enriched: list[str] = []
    monkeypatch.setattr(cli, "fetch_live_halo_bulge_sample", lambda cfg, cache, **kw: ([json.loads(json.dumps(s)) for s in stars], []))
    monkeypatch.setattr(cli, "enrich_star_live_for_vetting", lambda star, cfg, cache: (enriched.append(star["source_id"]) or dict(star), []))
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    monkeypatch.setattr("time.sleep", lambda s: None)
    run_dir = repo_root() / "runs" / "test-lamost-survey"
    try:
        argv = ["search", "spectral_lines", "--live", "--region", "halo", "--limit", "6", "--run-id", "test-lamost-survey", "--json"]
        code, payload = dispatch([*argv, "--defer-enrichment"])
        assert code == 0, payload
        d = payload["data"]
        assert d["n_input"] == 6 and d["n_searched"] == 2 and d["flagged_ids"] == ["2"]
        assert d["coverage"]["n_single_epoch"] == 2 and d["coverage"]["n_no_coverage"] == 1 and d["coverage"]["n_failed"] == 1
        assert d["spectra_fetch"]["n_ok"] == 2 and any("4 of 6 stars not searched" in w for w in payload["warnings"])
        # The flagged star is kept with its spectra although the archive enrichment was deferred
        assert enriched == []
        kept = json.loads((run_dir / "live_stars.json").read_text(encoding="utf-8"))
        assert list(kept) == ["2"] and len(kept["2"]["spectra"]["flux_exp1"]) == 101
        saved = json.loads((run_dir / "search_spectral_lines.json").read_text(encoding="utf-8"))
        assert saved["n_searched"] == 2 and saved["not_searched_ids"] == {"single_epoch": ["3", "6"], "no_coverage": ["4"], "failed": ["5"]}
        assert saved["spectra_fetch"]["failed_ids"] == ["5"] and len(saved["sample_stars"]) == 6

        code, payload = dispatch(["vet", "test-lamost-survey", "--json"])
        assert code == 0, payload
        vet = json.loads((run_dir / "vetting.json").read_text(encoding="utf-8"))["vetting"]["2"]
        status = {v["check"]: v["status"] for v in vet}
        assert status["cosmic_ray"] == "pass" and status["repeat"] == "pass" and status["detector"] == "not_evaluated"

        code, payload = dispatch(["null-report", "test-lamost-survey", "--json"])
        assert code == 0, payload
        rep = json.loads((run_dir / "null_report_spectral_lines.json").read_text(encoding="utf-8"))
        assert rep["n_searched"] == 2 and rep["n_flagged"] == 1

        rec, samples = build_run_record(run_dir)
        s = rec["searches"][0]
        assert s["coverage"]["n_searched"] == 2 and s["coverage"]["n_single_epoch"] == 2 and s["coverage"]["n_no_spectrum"] == 0
        assert s["not_searched_ids"] == ["3", "4", "5", "6"] and s["spectra_fetch"]["n_failed"] == 1 and samples[0]["n"] == 6
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)
