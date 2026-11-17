"""3D dust reddening at a star's distance (Bayestar19 via the Argonaut web API) for vetting enrichment.

Only called from vetting enrichment (`live_survey.enrich_star_live_for_vetting`), never from `_row_to_star_dict`,
so adding it never changes a certified channel score. The result is stored at
`star["vetting_archive"]["dust3d"]` and read by `physics.map3d_ebv` / `underluminous_checks.check_extinction`.

Request format (mirrors `dustmaps.bayestar.BayestarWebQuery(version="bayestar2019").query(coords, mode="percentile",
pct=[16, 50, 84])`): HTTP POST of JSON `{"coords": <SkyCoord in dustmaps' JSON encoding>, "mode": "percentile",
"pct": [16, 50, 84], "return_flags": false}` to `<API>/bayestar2019/query`. The response is decoded with the same
tagged-JSON conventions (`{"_type": "np.ndarray", ...}` in `data` base64, `value` or `npy` form).
LIVE-VERIFY: this payload/response format is reproduced from the dustmaps client source and has not been exercised
against the live server from this environment (no network); the parser also accepts plain JSON lists.

Units: Bayestar19 reports reddening in a unit that predicts the same E(g-r)_PS1 as one unit of SFD. Green et al.
(2019, as quoted on argonaut.skymaps.info/usage) give E(B-V) = 0.884 x Bayestar19 (via their Eq. 1, E(g-r)) or
0.996 x Bayestar19 (via Eq. 2, E(r-z)). The factor is the config key `vetting.bayestar19_ebv_per_unit` (default
0.884); the 0.884-0.996 spread is reported as `ebv_conversion_systematic`. LIVE-VERIFY the factor choice.
"""

from __future__ import annotations

import base64
import hashlib
import io
import json
import math
import time
from collections.abc import Callable
from typing import Any

import numpy as np
import requests

from oddstar.io.cache import RawCache
from oddstar.io.storage import equatorial_to_galactic
from oddstar.schemas import ManifestRecord

BAYESTAR19_QUERY_URL = "http://argonaut.skymaps.info/api/v2/bayestar2019/query"
BAYESTAR19_MIN_DEC_DEG = -30.0
BAYESTAR19_ALT_EBV_PER_UNIT = 0.996  # Green et al. 2019 Eq. 2 (E(r-z)) alternative, for the systematic term
PERCENTILES = (16.0, 50.0, 84.0)

Sleeper = Callable[[float], None]


def _ndarray_json(values: list[float]) -> dict[str, Any]:
    return {
        "_type": "np.ndarray",
        "dtype": {"_type": "np.dtype", "descr": "<f8"},
        "shape": [len(values)],
        "value": [float(v) for v in values],
    }


def build_bayestar_payload(l_deg: list[float], b_deg: list[float], d_kpc: list[float]) -> dict[str, Any]:
    """dustmaps-compatible JSON body for a batch of Galactic coordinates with distances."""

    def _q(vals: list[float], unit: str) -> dict[str, Any]:
        return {"_type": "astropy.units.Quantity", "value": _ndarray_json(vals), "unit": unit}

    coords = {
        "_type": "astropy.coordinates.SkyCoord",
        "frame": "galactic",
        "representation": "spherical",
        "lon": _q(l_deg, "deg"),
        "lat": _q(b_deg, "deg"),
        "distance": _q(d_kpc, "kpc"),
    }
    return {"coords": coords, "mode": "percentile", "pct": list(PERCENTILES), "return_flags": False}


def _decode(obj: Any) -> Any:
    """Decode dustmaps tagged JSON (np.ndarray / np.dtype / Quantity) into numpy arrays; pass plain JSON through."""
    if isinstance(obj, list):
        return [_decode(x) for x in obj]
    if not isinstance(obj, dict):
        return obj
    t = obj.get("_type")
    if t == "np.dtype":
        descr = obj.get("descr")
        return np.dtype(descr if isinstance(descr, str) else [tuple(x) for x in descr])
    if t == "np.ndarray":
        dt = _decode(obj.get("dtype", "<f8"))
        if "data" in obj:
            arr = np.frombuffer(base64.b64decode(obj["data"]), dtype=dt)
            return np.reshape(arr, obj.get("shape", arr.shape))
        if "value" in obj:
            return np.asarray(obj["value"], dtype=dt)
        if "npy" in obj:
            return np.load(io.BytesIO(json.loads(obj["npy"]).encode("latin-1")))
        raise ValueError("np.ndarray without data/value/npy")
    if t == "astropy.units.Quantity":
        return _decode(obj.get("value"))
    return {k: _decode(v) for k, v in obj.items()}


def parse_bayestar_percentiles(text: str, n_coords: int) -> np.ndarray:
    """Return an (n_coords, 3) array of the 16/50/84th percentile reddening (Bayestar19 units)."""
    raw = _decode(json.loads(text))
    if isinstance(raw, dict):
        for key in ("result", "data", "reddening", "percentile"):
            if key in raw:
                raw = raw[key]
                break
    arr = np.asarray(raw, dtype=float)
    if arr.ndim == 1 and n_coords == 1 and arr.size == len(PERCENTILES):
        arr = arr.reshape(1, -1)
    if arr.shape != (n_coords, len(PERCENTILES)):
        raise ValueError(f"unexpected Bayestar19 response shape {arr.shape}, expected {(n_coords, len(PERCENTILES))}")
    return arr


def _post_with_retries(url: str, body: bytes, *, timeout_s: float, attempts: int, sleep: Sleeper) -> requests.Response:
    errors: list[str] = []
    for i in range(attempts):
        try:
            resp = requests.post(url, data=body, headers={"content-type": "application/json"}, timeout=timeout_s)
            if resp.status_code == 200:
                return resp
            errors.append(f"HTTP {resp.status_code}")
            if resp.status_code < 500 and resp.status_code != 429:
                break
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
            errors.append(type(exc).__name__)
        if i < attempts - 1:
            sleep(2.0 * (2**i))
    raise RuntimeError(f"POST {url} failed after {len(errors)} attempt(s): {'; '.join(errors)}")


def query_bayestar19_at_distance(
    *,
    ra_deg: float,
    dec_deg: float,
    distance_pc: float | None,
    cache: RawCache,
    ebv_per_unit: float,
    parallax_mas: float | None = None,
    parallax_err_mas: float | None = None,
    distance_err_pc: float | None = None,
    timeout_s: float = 30.0,
    attempts: int = 3,
    sleep: Sleeper = time.sleep,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """E(B-V) from Bayestar19 at the star's distance, with its uncertainty. Never raises.

    sigma_d is `distance_err_pc`, else d * parallax_err / parallax. Queries the map at d - sigma_d, d and d + sigma_d
    in one batch request. ebv_err combines the map's 16-84th
    percentile half-width at d with half the change across d +/- sigma_d. Returns (record, manifests); the record
    always has a `status`: "ok", "no_coverage" (Dec < -30 deg or no map value), or "failed: <reason>".
    """
    if distance_err_pc is None and distance_pc is not None and parallax_mas and parallax_err_mas is not None:
        distance_err_pc = float(distance_pc) * float(parallax_err_mas) / max(1e-9, float(parallax_mas))
    base: dict[str, Any] = {
        "map": "bayestar19",
        "url": BAYESTAR19_QUERY_URL,
        "distance_pc": None if distance_pc is None else float(distance_pc),
        "distance_err_pc": None if distance_err_pc is None else float(distance_err_pc),
        "ebv_per_unit": float(ebv_per_unit),
        "ebv": None,
        "ebv_err": None,
    }
    if dec_deg < BAYESTAR19_MIN_DEC_DEG:
        return {**base, "status": "no_coverage", "note": "Bayestar19 covers Dec > -30 deg only"}, []
    if distance_pc is None or not math.isfinite(float(distance_pc)) or float(distance_pc) <= 0.0:
        return {**base, "status": "failed: no distance for this star"}, []
    try:
        l_deg, b_deg = equatorial_to_galactic(ra_deg, dec_deg)
        d = max(1.0, float(distance_pc))
        sd = max(0.0, float(distance_err_pc or 0.0))
        dists_kpc = [max(0.001, (d - sd) / 1000.0), d / 1000.0, (d + sd) / 1000.0]
        payload = build_bayestar_payload([l_deg] * 3, [b_deg] * 3, dists_kpc)
        body = json.dumps(payload, sort_keys=True).encode("utf-8")
        sig = hashlib.sha256(BAYESTAR19_QUERY_URL.encode("utf-8") + b"|" + body).hexdigest()[:24]
        manifests: list[ManifestRecord] = []
        cached = cache.lookup("dustmaps", "bayestar2019", sig)
        text: str | None = None
        if cached is not None and not cached.offline_fixture:
            try:
                text = cache.read_bytes(cached).decode("utf-8")
                manifests.append(cached)
            except FileNotFoundError:
                text = None
        if text is None:
            resp = _post_with_retries(BAYESTAR19_QUERY_URL, body, timeout_s=timeout_s, attempts=attempts, sleep=sleep)
            text = resp.text
            pct = parse_bayestar_percentiles(text, 3)  # validate before caching
            manifests.append(
                cache.put(
                    archive="dustmaps",
                    release="bayestar2019",
                    request_signature=sig,
                    query_text=f"POST {BAYESTAR19_QUERY_URL} {body.decode('utf-8')}",
                    payload=text.encode("utf-8"),
                    row_count=3,
                    offline_fixture=False,
                )
            )
        else:
            pct = parse_bayestar_percentiles(text, 3)
    except Exception as exc:  # noqa: BLE001 - enrichment must never abort; the status records the failure
        return {**base, "status": f"failed: {type(exc).__name__}: {exc}"}, []

    if not np.all(np.isfinite(pct[1])):
        return {**base, "status": "no_coverage", "note": "map returned no value at this position/distance", "l_deg": l_deg, "b_deg": b_deg}, manifests
    p16, p50, p84 = (float(x) for x in pct[1])
    map_half = 0.5 * (p84 - p16)
    lo_d = float(pct[0][1]) if math.isfinite(pct[0][1]) else p50
    hi_d = float(pct[2][1]) if math.isfinite(pct[2][1]) else p50
    dist_half = 0.5 * abs(hi_d - lo_d)
    err_units = math.hypot(map_half, dist_half)
    return (
        {
            **base,
            "status": "ok",
            "l_deg": l_deg,
            "b_deg": b_deg,
            "bayestar19_units_p16_p50_p84": [p16, p50, p84],
            "bayestar19_units_at_d_minus_plus_sigma": [lo_d, hi_d],
            "ebv": ebv_per_unit * p50,
            "ebv_err": ebv_per_unit * err_units,
            "ebv_conversion_systematic": abs(BAYESTAR19_ALT_EBV_PER_UNIT - ebv_per_unit) * p50,
            "note": (
                f"E(B-V) = {ebv_per_unit} x Bayestar19 (Green et al. 2019 Eq. 1; Eq. 2 gives {BAYESTAR19_ALT_EBV_PER_UNIT}); "
                "LIVE-VERIFY request format and factor"
            ),
        },
        manifests,
    )
