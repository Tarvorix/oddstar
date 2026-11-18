"""Radio narrowband vetting checks (Section 6).

Input rule (R9, lessons L007/L008): a check returns `pass` or `fail` only from inputs that were actually obtained.

- A live star (`live_enriched` key) counts an input only when `vetting_archive.query_status["breakthrough_listen"]
  == "ok"`. The live radio search (`oddstar.io.bl_radio`, config 1.5.0) writes that for a star whose ON and OFF
  fine-resolution scans were on disk and compared; it attaches the cadence of one window and `in_rfi_band` for its
  frequency. It does not produce a session cross-check or a second-epoch record, so those two checks stay
  `not_evaluated` for a live star. For every other live star (no beam contains it, the scans were not downloaded,
  no OFF scan exists, the search failed) every check returns `not_evaluated` and names the status.
- For a cadence or hand-built dict an input counts only when its key is explicitly present and not None.
  `off_scan_hits` and `drift_hz_s` are written by `RadioNarrowbandChannel.measure` (or copied from the stored
  channel measurement by `oddstar vet`); if the measurement has not run they are absent and the check reports it.

Every evaluated result stores the input values it used in its metrics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result

_PRODUCTS = {
    "off_scan_hits": "de-Doppler search of the OFF scans (channel measurement)",
    "drift_hz_s": "measured drift rate of the ON-scan signal (channel measurement)",
    "in_rfi_band": "RFI allocation lookup for the signal frequency",
    "multi_target_match": "session cross-check against unrelated targets",
    "reobserved_confirmed": "independent second-epoch reobservation record",
}
RADIO_QUERY = "breakthrough_listen"
NO_LIVE_DATA = "no ON/OFF radio scans were compared for this star"
NOT_FROM_LIVE_SEARCH = "the live radio search does not produce this input"


def _is_live(target: dict[str, Any]) -> bool:
    return "live_enriched" in target


def _live_status(target: dict[str, Any]) -> str:
    return ((target.get("vetting_archive") or {}).get("query_status") or {}).get(RADIO_QUERY) or "not run"


def _input(target: dict[str, Any], key: str) -> Any:
    """Explicit `target[key]`; for a live star only when `query_status[RADIO_QUERY]` is "ok" (its scans were compared)."""
    if _is_live(target) and _live_status(target) != "ok":
        return None
    return target.get(key)


def _not_evaluated(out_dir: Path, check: str, target: dict[str, Any], key: str, extra: str = "") -> VettingResult:
    note = f"Missing input: {_PRODUCTS[key]} [{key}]"
    if _is_live(target):
        status = _live_status(target)
        note += f" (query_status {RADIO_QUERY}={status}; {NOT_FROM_LIVE_SEARCH if status == 'ok' else NO_LIVE_DATA})"
    elif "scans" not in target:
        note += " (no ON/OFF filterbank cadence for this target)"
    return make_vet_result(
        out_dir, check=check, status="not_evaluated", note=note + extra,
        metrics={"input_values": {key: None}, "missing_inputs": [key], "live_star": _is_live(target), "has_cadence": "scans" in target},
    )


def check_off_scan(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    raw = _input(target, "off_scan_hits")
    if raw is None:
        return _not_evaluated(out_dir, "off_scan", target, "off_scan_hits")
    hits_off = int(raw)
    n_off = sum(1 for sc in target.get("scans") or [] if sc.get("kind") == "OFF")
    return make_vet_result(
        out_dir,
        check="off_scan",
        status="pass" if hits_off == 0 else "fail",
        note=f"Off-target scan detections at matching frequency/drift: {hits_off} (of {n_off} OFF scans in the cadence)",
        metrics={"input_values": {"off_scan_hits": hits_off, "n_off_scans": n_off}, "off_scan_hits": hits_off},
    )


def check_zero_drift(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    raw = _input(target, "drift_hz_s")
    if raw is None:
        return _not_evaluated(out_dir, "zero_drift", target, "drift_hz_s")
    drift = float(raw)
    min_abs = float(cfg.channels["radio_narrowband"]["min_drift_abs_hz_s"])
    return make_vet_result(
        out_dir,
        check="zero_drift",
        status="pass" if abs(drift) >= min_abs else "fail",
        note=f"Measured drift rate={drift:.4f} Hz/s (min |drift|={min_abs} Hz/s)",
        metrics={"input_values": {"drift_hz_s": drift}, "drift_hz_s": drift, "min_abs_drift_hz_s": min_abs},
    )


def check_rfi_band(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    raw = _input(target, "in_rfi_band")
    if raw is None:
        return _not_evaluated(out_dir, "rfi_band", target, "in_rfi_band")
    in_rfi = bool(raw)
    freq = target.get("center_freq_mhz")
    return make_vet_result(
        out_dir,
        check="rfi_band",
        status="fail" if in_rfi else "pass",
        note=(f"Signal at {freq} MHz falls inside known terrestrial/satellite RFI band" if in_rfi
              else f"Center frequency {freq} MHz outside known satellite/terrestrial RFI allocations"),
        metrics={"input_values": {"in_rfi_band": in_rfi, "center_freq_mhz": freq}, "in_rfi_band": in_rfi, "center_freq_mhz": freq},
    )


def check_multi_target(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    raw = _input(target, "multi_target_match")
    if raw is None:
        return _not_evaluated(out_dir, "multi_target", target, "multi_target_match")
    multi = bool(raw)
    return make_vet_result(
        out_dir,
        check="multi_target",
        status="fail" if multi else "pass",
        note="Same frequency/drift signature observed toward unrelated targets in observing session" if multi else "Signal absent toward unrelated targets in session",
        metrics={"input_values": {"multi_target_match": multi}, "multi_target_match": multi},
    )


def check_reobservation(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Pass only on a confirmed independent second epoch; no confirmation is an open hole (caps tier at C), not a fail."""
    raw = _input(target, "reobserved_confirmed")
    if raw is None:
        return _not_evaluated(out_dir, "reobservation", target, "reobserved_confirmed", " (caps tier at C)")
    reobs = bool(raw)
    return make_vet_result(
        out_dir,
        check="reobservation",
        status="pass" if reobs else "not_evaluated",
        note="Confirmed in independent second-epoch radio cadence" if reobs else "No independent second-epoch reobservation confirmed (caps tier at C)",
        metrics={"input_values": {"reobserved_confirmed": reobs}, "reobserved_confirmed": reobs},
    )


def run_radio_vetting(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [
        check_off_scan,
        check_zero_drift,
        check_rfi_band,
        check_multi_target,
        check_reobservation,
    ]
    return [fn(target, cfg, out_dir) for fn in checks]
