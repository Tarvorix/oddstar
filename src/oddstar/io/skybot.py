"""IMCCE SkyBoT cone search: known Solar System bodies in a field at a given epoch.

Public API (shared with the light-curve vetting checks, keep the signature stable):

    skybot_cone_search(*, ra_deg, dec_deg, radius_arcsec, epoch_jd, observer="500", cache, timeout_s=30.0,
                       sleep=time.sleep) -> (bodies, manifest)

Each body is a dict with at least `name`, `ra_deg`, `dec_deg`, `sep_arcsec` (from the cone centre at `epoch_jd`)
and `v_mag`, plus `number`, `class`, `pos_err_arcsec`, `dra_arcsec_hr` (dRA*cos(Dec)), `ddec_arcsec_hr`,
`dg_au`, `dh_au` and `epoch_jd` when SkyBoT returns them. HTTP failures (after bounded retries) and
unparseable responses raise `ArchiveAccessError`; an empty field returns an empty list.

Endpoint: https://ssp.imcce.fr/webservices/skybot/api/conesearch.php with `-ep` (JD), `-ra`/`-dec` (deg),
`-rd` (deg), `-mime=text`, `-output=basic`, `-loc` (IAU/MPC observatory code, e.g. 500 = geocentre,
C51 = WISE), `-filter=0`, `-objFilter=111` (asteroids, planets, comets) and `-refsys=EQJ2000`.
The text response is a `#`-commented, `|`-separated table preceded by `# Flag: <n>` (1 = bodies found,
0 = none, -1 = error). A JSON body (if the service answers JSON) is also accepted.
"""

from __future__ import annotations

import json
import math
import re
import time
from typing import Any

from oddstar.errors import ArchiveAccessError
from oddstar.io.archive_http import Sleeper, cached_get
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec
from oddstar.schemas import ManifestRecord

SKYBOT_CONESEARCH_URL = "https://ssp.imcce.fr/webservices/skybot/api/conesearch.php"


def _sexa_to_deg(text: str, *, hours: bool) -> float:
    parts = text.replace(":", " ").split()
    if len(parts) == 1:
        val = float(parts[0])
        return val * 15.0 if hours else val
    sign = -1.0 if parts[0].strip().startswith("-") else 1.0
    d = abs(float(parts[0]))
    m = float(parts[1]) if len(parts) > 1 else 0.0
    s = float(parts[2]) if len(parts) > 2 else 0.0
    val = sign * (d + m / 60.0 + s / 3600.0)
    return val * 15.0 if hours else val


def _num(text: Any) -> float | None:
    if text is None:
        return None
    try:
        v = float(str(text).strip())
    except ValueError:
        return None
    return None if math.isnan(v) else v


def _norm(key: str) -> str:
    return re.sub(r"[^a-z]", "", key.lower())


# Normalised header -> output key. Text headers look like "RA(h)", "DE(deg)", "Mv", "Err(arcsec)", "d(arcsec)",
# "dRA(arcsec/h)", "dDEC(arcsec/h)", "Dg(ua)", "Dh(ua)".
_TEXT_COLS = {
    "num": "number",
    "name": "name",
    "rah": "ra",
    "dedeg": "dec",
    "class": "class",
    "mv": "v_mag",
    "errarcsec": "pos_err_arcsec",
    "darcsec": "sep_arcsec",
    "draarcsech": "dra_arcsec_hr",
    "ddecarcsech": "ddec_arcsec_hr",
    "dgua": "dg_au",
    "dhua": "dh_au",
}


def parse_skybot_text(text: str, *, ra_deg: float, dec_deg: float, epoch_jd: float) -> list[dict[str, Any]]:
    """Parse the SkyBoT `-mime=text` response. Raises ArchiveAccessError on an error flag or unknown format."""
    stripped = text.strip()
    if stripped.startswith("[") or stripped.startswith("{"):
        return _parse_skybot_json(stripped, ra_deg=ra_deg, dec_deg=dec_deg, epoch_jd=epoch_jd)
    flag: int | None = None
    header: list[str] | None = None
    rows: list[dict[str, Any]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            m = re.match(r"#\s*Flag\s*:\s*(-?\d+)", line)
            if m:
                flag = int(m.group(1))
                continue
            if "|" in line and "name" in line.lower():
                header = [_norm(c) for c in line.lstrip("#").split("|")]
            continue
        if flag is None:
            break
        if "|" not in line or header is None:
            continue
        vals = [v.strip() for v in line.split("|")]
        raw = {_TEXT_COLS.get(h, h): v for h, v in zip(header, vals)}
        try:
            b_ra = _sexa_to_deg(raw["ra"], hours=True)
            b_dec = _sexa_to_deg(raw["dec"], hours=False)
        except (KeyError, ValueError, IndexError):
            continue
        rows.append(_finish_body(raw, b_ra, b_dec, ra_deg=ra_deg, dec_deg=dec_deg, epoch_jd=epoch_jd))
    if flag is None:
        raise ArchiveAccessError(f"SkyBoT response not understood: {stripped[:160]!r}", retryable=True)
    if flag < 0:
        raise ArchiveAccessError(f"SkyBoT returned error flag {flag}: {stripped[:200]!r}", retryable=False)
    return rows


def _parse_skybot_json(text: str, *, ra_deg: float, dec_deg: float, epoch_jd: float) -> list[dict[str, Any]]:
    try:
        js = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ArchiveAccessError(f"SkyBoT JSON parse failed: {exc}", retryable=True) from exc
    items = js.get("data", js.get("results", [])) if isinstance(js, dict) else js
    if isinstance(js, dict) and int(js.get("flag", 1) or 1) < 0:
        raise ArchiveAccessError(f"SkyBoT returned error: {text[:200]!r}", retryable=False)
    rows: list[dict[str, Any]] = []
    for it in items or []:
        raw = {_TEXT_COLS.get(_norm(k), _norm(k)): v for k, v in dict(it).items()}
        dec_v = raw.get("dec", raw.get("decdeg", raw.get("dedeg")))
        try:
            if _num(raw.get("radeg")) is not None:  # explicit degrees
                b_ra = float(raw["radeg"])
            elif raw.get("ra") is not None:  # SkyBoT "RA(h)": hours, sexagesimal or decimal
                b_ra = _sexa_to_deg(str(raw["ra"]), hours=True)
            else:
                continue
            if dec_v is None:
                continue
            b_dec = _sexa_to_deg(str(dec_v), hours=False)
        except ValueError:
            continue
        rows.append(_finish_body(raw, b_ra, b_dec, ra_deg=ra_deg, dec_deg=dec_deg, epoch_jd=epoch_jd))
    return rows


def _finish_body(raw: dict[str, Any], b_ra: float, b_dec: float, *, ra_deg: float, dec_deg: float, epoch_jd: float) -> dict[str, Any]:
    sep = angular_separation_arcsec(ra_deg, dec_deg, b_ra, b_dec)
    return {
        "number": str(raw.get("number") or "").strip() or None,
        "name": str(raw.get("name") or "").strip(),
        "ra_deg": b_ra,
        "dec_deg": b_dec,
        "sep_arcsec": sep,
        "v_mag": _num(raw.get("v_mag")),
        "class": str(raw.get("class") or "").strip(),
        "pos_err_arcsec": _num(raw.get("pos_err_arcsec")),
        "dra_arcsec_hr": _num(raw.get("dra_arcsec_hr")),
        "ddec_arcsec_hr": _num(raw.get("ddec_arcsec_hr")),
        "dg_au": _num(raw.get("dg_au")),
        "dh_au": _num(raw.get("dh_au")),
        "epoch_jd": float(epoch_jd),
    }


def skybot_cone_search(
    *,
    ra_deg: float,
    dec_deg: float,
    radius_arcsec: float,
    epoch_jd: float,
    observer: str = "500",
    cache: RawCache,
    timeout_s: float = 30.0,
    sleep: Sleeper = time.sleep,
) -> tuple[list[dict[str, Any]], ManifestRecord]:
    """Known Solar System bodies within `radius_arcsec` of (ra, dec) at `epoch_jd` seen from `observer`."""
    params = {
        "-ep": f"{epoch_jd:.6f}",
        "-ra": f"{ra_deg:.6f}",
        "-dec": f"{dec_deg:+.6f}",
        "-rd": f"{radius_arcsec / 3600.0:.6f}",
        "-mime": "text",
        "-output": "basic",
        "-loc": str(observer),
        "-filter": "0",
        "-objFilter": "111",
        "-refsys": "EQJ2000",
    }

    def _validate(payload: bytes) -> None:
        parse_skybot_text(payload.decode("utf-8", errors="replace"), ra_deg=ra_deg, dec_deg=dec_deg, epoch_jd=epoch_jd)

    payload, man, _ = cached_get(
        url=SKYBOT_CONESEARCH_URL,
        params=params,
        archive="ephemeris",
        release="skybot",
        cache=cache,
        timeout_s=timeout_s,
        attempts=3,
        sleep=sleep,
        validate=_validate,
    )
    bodies = parse_skybot_text(payload.decode("utf-8", errors="replace"), ra_deg=ra_deg, dec_deg=dec_deg, epoch_jd=epoch_jd)
    return bodies, man
