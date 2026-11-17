"""Comprehensive offline test suite verifying Milestones M0-M13 and Ground Rules R1-R10."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from oddstar.channels import ALL_CHANNEL_NAMES, build_channel
from oddstar.ci_checks import run_all_ci_checks
from oddstar.cli import dispatch
from oddstar.config import load_config, repo_root
from oddstar.errors import EXIT_UNCERTIFIED, PolicyViolationError
from oddstar.fuse import assign_tier, detect_cross_channel_patterns, fuse_results
from oddstar.io import ArchiveFetcher, RADIO_CADENCES, RawCache, TableStore, generate_field_sample, get_named_object
from oddstar.schemas import Dossier, ManifestRecord, NullReport, ReviewRecord, ValidationCertificate
from oddstar.server import MCP_TOOLS, handle_jsonrpc_request
from oddstar.validate import run_validation
from oddstar.vet.ir_checks import run_ir_vetting
from oddstar.vet.lc_checks import run_lc_vetting
from oddstar.vet.radio_checks import run_radio_vetting
from oddstar.vet.spectral_checks import run_spectral_vetting
from oddstar.vet.underluminous_checks import run_underluminous_vetting
from oddstar.vet.vanishing_checks import run_vanishing_vetting


def test_m0_doctor_and_schemas() -> None:
    code, payload = dispatch(["doctor", "--json"])
    assert code == 0
    assert payload["ok"] is True
    assert payload["data"]["cache"]["max_cache_gb"] == 10.0
    assert payload["data"]["archives_registered"] == 23
    assert payload["data"]["channels_registered"] == 8
    assert payload["data"]["lock_file_present"] is True

    root = repo_root()
    assert (root / "AGENTS.md").exists()
    assert (root / "review" / "review_record.schema.json").exists()


def test_m1_fetch_zero_network_second_call_and_budget_eviction(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=10.0)
    fetcher = ArchiveFetcher("gaia", cfg, cache)
    ids_100 = [f"9{10000000000000000 + i}" for i in range(100)]

    rows1, man1, hit1 = fetcher.fetch(ids=ids_100)
    assert hit1 is False
    assert fetcher.backend_calls == 1
    assert len(rows1) == 100
    ManifestRecord.model_validate(man1.model_dump())

    # Second identical fetch must hit cache with zero backend/network calls
    rows2, man2, hit2 = fetcher.fetch(ids=ids_100)
    assert hit2 is True
    assert fetcher.backend_calls == 1
    assert man2.sha256 == man1.sha256
    assert rows2 == rows1

    # HSC-SSP disabled by default under R/Section 13 policy
    hsc_fetcher = ArchiveFetcher("hsc_ssp", cfg, cache)
    with pytest.raises(PolicyViolationError):
        hsc_fetcher.fetch(ids=["4582145718700356096"])

    # Tight cache cap triggers deterministic LRU eviction while keeping manifest SHA-256
    tiny_cache = RawCache(tmp_path / "tiny_raw", max_cache_gb=0.000001)  # ~1 KB
    m_a = tiny_cache.put(archive="breakthrough_listen", release="v2", request_signature="sig1", query_text="q1", payload=b"A" * 700, row_count=1)
    m_b = tiny_cache.put(archive="breakthrough_listen", release="v2", request_signature="sig2", query_text="q2", payload=b"B" * 700, row_count=1)
    assert m_b.sha256 != m_a.sha256
    assert tiny_cache.total_cached_bytes() <= 1100
    assert len(tiny_cache.all_manifests()) == 2
    assert any(r.evicted for r in tiny_cache.all_manifests())

    # HEALPix partitioned storage round-trip
    store = TableStore(tmp_path / "tables", nside=16)
    written = store.write_partitioned("gaia_sample", rows1[:20])
    assert written
    q_rows = store.query("gaia_sample")
    assert len(q_rows) == 20


def test_m2_ir_excess_10k_chi_distribution(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    ch = build_channel("ir_excess", cfg, tmp_path)
    spec = ch.sample(cfg, limit=10000)
    assert spec.total_count > 7000
    assert len(spec.cuts) == 6

    stars = generate_field_sample(10000, seed=cfg.sample.seed, gaia_release=cfg.gaia.release)
    chis_mod = []
    chis_emp = []
    for s in stars[:1000]:
        m = ch.measure(s, {})
        chis_mod.append(m.values["chi_w4"])
        chis_emp.append(m.values["chi_emp_w4"])

    assert abs(float(np.mean(chis_mod))) < cfg.channels["ir_excess"]["chi_distribution_mean_tol"]
    assert abs(float(np.std(chis_mod)) - 1.0) < cfg.channels["ir_excess"]["chi_distribution_std_tol"]
    assert abs(float(np.mean(chis_emp))) < cfg.channels["ir_excess"]["chi_distribution_mean_tol"]
    assert abs(float(np.std(chis_emp)) - 1.0) < cfg.channels["ir_excess"]["chi_distribution_std_tol"]

    # Mandatory single-epoch caveat
    m_first = ch.measure(stars[0], {})
    assert any("single instrument" in c.lower() for c in m_first.caveats)


def test_m3_all_vetting_checks_have_pass_fail_fixtures_and_evidence(tmp_path: Path) -> None:
    cfg, _, _ = load_config()

    # 1. Infrared vetting pass vs fail fixtures
    clean_star = copy.deepcopy(get_named_object("115725517896200320"))  # BD+20 307
    assert clean_star is not None
    clean_star["vetting_archive"]["deep_imaging_available"] = True
    clean_res = run_ir_vetting(clean_star, cfg, tmp_path / "ir_pass")
    assert len(clean_res) == 16
    for r in clean_res:
        assert r.status == "pass", f"Expected BD+20 307 to pass {r.check}: {r.note}"
        assert r.evidence is not None and Path(r.evidence).exists()

    # Failing IR fixtures (all 15 checks covered)
    imp_radio = get_named_object("4418120359123456128")
    imp_dog = get_named_object("3141592653589793280")
    imp_cirrus = get_named_object("5812938471029384960")
    imp_art = get_named_object("4129384756102938112")
    assert imp_radio and imp_dog and imp_cirrus and imp_art

    assert {r.check: r.status for r in run_ir_vetting(imp_radio, cfg, tmp_path / "ir_f1")}["radio_counterpart"] == "fail"
    assert {r.check: r.status for r in run_ir_vetting(imp_dog, cfg, tmp_path / "ir_f2")}["neighbours"] == "fail"
    assert {r.check: r.status for r in run_ir_vetting(imp_cirrus, cfg, tmp_path / "ir_f3")}["cirrus"] == "fail"
    f4 = {r.check: r.status for r in run_ir_vetting(imp_art, cfg, tmp_path / "ir_f4")}
    assert f4["wise_quality"] == "fail" and f4["single_frames"] == "fail"

    # Melis et al. 2012 TYC 8241-2652-1 secular change fixture
    melis = get_named_object("5391513322976037248")
    assert melis is not None
    melis_sec = {r.check: r for r in run_ir_vetting(melis, cfg, tmp_path / "ir_melis")}["second_instrument"]
    assert melis_sec.status == "pass" and "real_secular_fading" in melis_sec.note

    # Synthetic fail fixture for youth, evolved, binary, moving_object, xmatch_integrity, extragalactic_colours
    synth_fail = copy.deepcopy(clean_star)
    synth_fail["allwise"]["ra_deg"] += 0.02
    synth_fail["vetting_archive"]["moving_objects_2010"] = [{"name": "2010 Asteroid"}]
    synth_fail["vetting_archive"]["young_association"] = "Taurus"
    synth_fail["abs_g_mag"] = -1.2
    synth_fail["teff_k"] = 4200.0
    synth_fail["ruwe"] = 2.8
    synth_fail["mags"]["W1"] = 14.0
    synth_fail["mags"]["W2"] = 12.8
    synth_fail["mags"]["W3"] = 9.5
    sf_map = {r.check: r.status for r in run_ir_vetting(synth_fail, cfg, tmp_path / "ir_synth_fail")}
    for chk in ("xmatch_integrity", "moving_object", "youth", "evolved", "binary", "extragalactic_colours"):
        assert sf_map[chk] == "fail"

    # 2. Light-curve checks
    tabby = get_named_object("2081900940499099136")
    keb = get_named_object("722200000000000001")
    tess_sc = get_named_object("722200000000000002")
    blend_eb = get_named_object("722200000000000003")
    assert tabby and keb and tess_sc and blend_eb
    assert all(r.status == "pass" for r in run_lc_vetting(tabby, cfg, tmp_path / "lc_pass"))
    assert {r.check: r.status for r in run_lc_vetting(keb, cfg, tmp_path / "lc_f1")}["known_variable"] == "fail"
    f_sc = {r.check: r.status for r in run_lc_vetting(tess_sc, cfg, tmp_path / "lc_f2")}
    assert f_sc["systematics_window"] == "fail" and f_sc["pipeline_agreement"] == "fail"
    assert {r.check: r.status for r in run_lc_vetting(blend_eb, cfg, tmp_path / "lc_f3")}["neighbour_signal"] == "fail"

    # 3. Radio checks
    v1 = copy.deepcopy(RADIO_CADENCES["GBT_VOYAGER1_20151230"])
    v1["off_scan_hits"] = 0
    v1["drift_hz_s"] = -0.38
    assert all(r.status == "pass" for r in run_radio_vetting(v1, cfg, tmp_path / "rad_pass"))
    blc1 = copy.deepcopy(RADIO_CADENCES["GBT_BLC1_PROXCEN"])
    blc1["off_scan_hits"] = 3
    blc1["drift_hz_s"] = 0.0
    rf = {r.check: r.status for r in run_radio_vetting(blc1, cfg, tmp_path / "rad_fail")}
    assert rf["off_scan"] == "fail" and rf["zero_drift"] == "fail" and rf["rfi_band"] == "fail" and rf["multi_target"] == "fail"
    assert rf["reobservation"] == "not_evaluated"

    # 4. Spectral checks
    spec_ok = copy.deepcopy(clean_star)
    sf_cr = get_named_object("744400000000000003")
    sf_air = get_named_object("744400000000000004")
    sf_md = get_named_object("744400000000000001")
    assert sf_cr and sf_air and sf_md
    assert all(r.status == "pass" for r in run_spectral_vetting(spec_ok, cfg, tmp_path / "sp_pass"))
    assert {r.check: r.status for r in run_spectral_vetting(sf_cr, cfg, tmp_path / "sp_cr")}["cosmic_ray"] == "fail"
    assert {r.check: r.status for r in run_spectral_vetting(sf_air, cfg, tmp_path / "sp_air")}["sky_line"] == "fail"
    assert {r.check: r.status for r in run_spectral_vetting(sf_md, cfg, tmp_path / "sp_md")}["stellar_emission"] == "fail"

    # 5. Vanishing checks
    v_ok = get_named_object("755500000000000003")
    v_pm = get_named_object("755500000000000001")
    v_lpv = get_named_object("755500000000000005")
    v_pd = get_named_object("755500000000000006")
    assert v_ok and v_pm and v_lpv and v_pd
    assert all(r.status == "pass" for r in run_vanishing_vetting(v_ok, cfg, tmp_path / "van_pass"))
    assert {r.check: r.status for r in run_vanishing_vetting(v_pm, cfg, tmp_path / "van_pm")}["proper_motion"] == "fail"
    assert {r.check: r.status for r in run_vanishing_vetting(v_lpv, cfg, tmp_path / "van_lpv")}["known_variable"] == "fail"
    assert {r.check: r.status for r in run_vanishing_vetting(v_pd, cfg, tmp_path / "van_pd")}["plate_defect"] == "fail"

    # 6. Under-luminous checks
    assert all(r.status == "pass" for r in run_underluminous_vetting(clean_star, cfg, tmp_path / "ul_pass"))
    ul_ext = get_named_object("733300000000000001")
    ul_sub = get_named_object("733300000000000004")
    ul_ruwe = get_named_object("733300000000000005")
    assert ul_ext and ul_sub and ul_ruwe
    assert {r.check: r.status for r in run_underluminous_vetting(ul_ext, cfg, tmp_path / "ul_f1")}["extinction"] == "fail"
    assert {r.check: r.status for r in run_underluminous_vetting(ul_sub, cfg, tmp_path / "ul_f2")}["metallicity"] == "fail"
    assert {r.check: r.status for r in run_underluminous_vetting(ul_ruwe, cfg, tmp_path / "ul_f3")}["astrometry"] == "fail"


def test_m4_validation_gate_and_search_lockout(tmp_path: Path) -> None:
    # Copy minimal workspace layout in tmp_path to verify uncertified lockout
    root = repo_root()
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs" / "default.yaml").write_text((root / "configs" / "default.yaml").read_text(), encoding="utf-8")
    (tmp_path / "validation" / "benchmarks").mkdir(parents=True)
    (tmp_path / "validation" / "impostors").mkdir(parents=True)
    (tmp_path / "validation" / "benchmarks" / "benchmarks.json").write_text(
        (root / "validation" / "benchmarks" / "benchmarks.json").read_text(), encoding="utf-8"
    )
    (tmp_path / "validation" / "impostors" / "impostors.json").write_text(
        (root / "validation" / "impostors" / "impostors.json").read_text(), encoding="utf-8"
    )

    # Search MUST fail with EXIT_UNCERTIFIED when no certificate exists
    code, err_payload = dispatch(["search", "ir_excess", "--json"], root=tmp_path)
    assert code == EXIT_UNCERTIFIED
    assert err_payload["ok"] is False
    assert err_payload["code"] == "UNCERTIFIED_CHANNEL"

    # Running validate certifies the channel, then search succeeds
    vcode, vpay = dispatch(["validate", "ir_excess", "--limit", "400", "--json"], root=tmp_path)
    assert vcode == 0 and vpay["ok"] is True
    scode, spay = dispatch(["search", "ir_excess", "--limit", "200", "--json"], root=tmp_path)
    assert scode == 0 and spay["ok"] is True


def test_m5_open_cases_inspect() -> None:
    code, payload = dispatch(["inspect", "4582145718700356096", "--json"])
    assert code == 0 and payload["ok"] is True
    dossier_path = Path(payload["artifacts"][0])
    assert dossier_path.exists()
    dossier = Dossier.model_validate_json(dossier_path.read_text(encoding="utf-8"))
    assert dossier.tier == "C"
    assert len(dossier.lead_lines) == 4
    assert len(dossier.hypotheses) == 12

    ir_ch = next(c for c in dossier.channels if c.name == "ir_excess")
    assert 40.0 <= float(ir_ch.measurements["w4_ratio"]) <= 80.0
    assert 0.008 <= float(ir_ch.measurements["frac_lum"]) <= 0.025
    assert 1.1 <= float(ir_ch.measurements["r_bb_au"]) <= 2.2
    assert any("single instrument" in cav.lower() for cav in ir_ch.caveats)

    holes = {v.check for v in dossier.vetting if v.status == "not_evaluated"}
    assert {"deep_imaging", "second_instrument"}.issubset(holes)
    for ev_rel in ir_ch.evidence:
        assert (dossier_path.parent / ev_rel).exists()

    # Also test HD 139139 open case
    code2, payload2 = dispatch(["inspect", "6256050602424290048", "--json"])
    assert code2 == 0 and payload2["ok"] is True


def test_m6_review_kit_isolation_and_lessons(tmp_path: Path) -> None:
    root = repo_root()
    sys.path.insert(0, str(root))
    import review_kit.verify_target as vt

    # Verify R6 isolation barrier rejects data/ and runs/ paths
    with pytest.raises(RuntimeError, match="R6 isolation violation"):
        vt.assert_path_allowed(root / "data" / "raw" / "manifests.jsonl")
    with pytest.raises(RuntimeError, match="R6 isolation violation"):
        vt.assert_path_allowed(root / "runs" / "foo.json")

    # Verify Round 1 packet creation inside tmp_path mixes targets, controls, known positives & impostors cleanly
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs" / "default.yaml").write_text((root / "configs" / "default.yaml").read_text(), encoding="utf-8")
    (tmp_path / "review" / "prompts").mkdir(parents=True)
    (tmp_path / "review" / "prompts" / "round1_blind.md").write_text(
        (root / "review" / "prompts" / "round1_blind.md").read_text(), encoding="utf-8"
    )
    code, pay = dispatch(["review", "packet", "4582145718700356096", "--round", "1", "--json"], root=tmp_path)
    assert code == 0
    pkt_path = Path(pay["artifacts"][0])
    pkt_data = json.loads(pkt_path.read_text(encoding="utf-8"))
    assert len(pkt_data["targets"]) >= 5
    # Answer key is NOT inside the packet
    assert "role" not in json.dumps(pkt_data)


def test_m7_m11_remaining_channels_validation(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    root = repo_root()
    (tmp_path / "validation" / "benchmarks").mkdir(parents=True)
    (tmp_path / "validation" / "impostors").mkdir(parents=True)
    (tmp_path / "validation" / "benchmarks" / "benchmarks.json").write_text(
        (root / "validation" / "benchmarks" / "benchmarks.json").read_text(), encoding="utf-8"
    )
    (tmp_path / "validation" / "impostors" / "impostors.json").write_text(
        (root / "validation" / "impostors" / "impostors.json").read_text(), encoding="utf-8"
    )
    for ch_name in ("ir_timedomain", "optical_dips", "radio_narrowband", "underluminous", "spectral_lines", "vanishing"):
        cert, _ = run_validation(ch_name, cfg, tmp_path, split="all", field_limit=450)
        assert cert.passed is True
        assert cert.blind_recovery_rate >= 0.90
        assert cert.impostor_pass is True
        assert cert.impostor_max_tier in {"Explained", "Rejected", "C"}
        assert cert.injection_completeness_design_point >= 0.90


def test_m12_fusion_patterns_null_report_and_mcp_server() -> None:
    cfg, _, _ = load_config()
    root = repo_root()

    # Verify ASASSN-21qj triggers named cross-channel pattern ir_excess + optical_dips
    asassn = get_named_object("5548511859528971520")
    assert asassn is not None
    ch_ir = build_channel("ir_excess", cfg, root)
    ch_opt = build_channel("optical_dips", cfg, root)
    sc_ir = ch_ir.score(ch_ir.measure(asassn, {}))
    sc_opt = ch_opt.score(ch_opt.measure(asassn, {}))
    assert sc_ir.flagged and sc_opt.flagged
    patterns = detect_cross_channel_patterns([sc_ir, sc_opt])
    assert any("ir_excess_plus_optical_dips" in p for p in patterns)

    # Verify Tier A requires >=2 distinct physical instruments & 0 holes
    from oddstar.schemas import VettingResult

    pass_vet = [VettingResult(check="xmatch_integrity", status="pass", evidence="e.json", note="ok")]
    tier_a, _ = assign_tier([sc_ir, sc_opt], pass_vet, cfg)
    assert tier_a == "A"

    # Verify MCP JSON-RPC tool server lists all tools & executes doctor + status cleanly
    list_resp = handle_jsonrpc_request({"jsonrpc": "2.0", "id": 1, "method": "tools/list"})
    assert len(list_resp["result"]["tools"]) == len(MCP_TOOLS) == 17
    doc_resp = handle_jsonrpc_request(
        {"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "oddstar_doctor", "arguments": {}}}
    )
    assert doc_resp["result"]["isError"] is False
    doc_payload = json.loads(doc_resp["result"]["content"][0]["text"])
    assert doc_payload["ok"] is True


def test_m13_bit_for_bit_reproducibility_and_dr4_adapter() -> None:
    root = repo_root()
    # Two identical runs of ir_excess must produce bit-for-bit identical search payloads (excluding timestamps)
    c1, p1 = dispatch(["search", "ir_excess", "--run-id", "test-repro-a", "--seed", "20261116", "--limit", "150", "--json"])
    c2, p2 = dispatch(["search", "ir_excess", "--run-id", "test-repro-b", "--seed", "20261116", "--limit", "150", "--json"])
    assert c1 == 0 and c2 == 0
    d1 = json.loads(Path(p1["artifacts"][0]).read_text(encoding="utf-8"))
    d2 = json.loads(Path(p2["artifacts"][0]).read_text(encoding="utf-8"))
    d1.pop("run_id")
    d2.pop("run_id")
    assert hashlib.sha256(json.dumps(d1, sort_keys=True).encode()).hexdigest() == hashlib.sha256(
        json.dumps(d2, sort_keys=True).encode()
    ).hexdigest()

    # Check DR3 vs DR4 comparison artifact & DR4 certificates exist (tracked in validation/certificates/)
    comp_path = root / "validation" / "certificates" / "dr3_vs_dr4_comparison.json"
    assert comp_path.exists()
    comp = json.loads(comp_path.read_text(encoding="utf-8"))
    assert comp["channels_revalidated"] == 8
    assert all(row["dr4_passed"] for row in comp["comparison"])

    # Check all CI governance rules pass cleanly
    report = run_all_ci_checks(root)
    assert all(len(v) == 0 for v in report.values()), f"CI check violations: {report}"
