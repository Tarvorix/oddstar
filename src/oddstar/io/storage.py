"""HEALPix spherical indexing and partitioned table storage with SQL query engine."""

from __future__ import annotations

import json
import math
import sqlite3
from pathlib import Path
from typing import Any


def ra_dec_to_healpix(ra_deg: float, dec_deg: float, nside: int = 16) -> int:
    """Compute exact HEALPix RING pixel index (Górski et al. 2005) from (RA, Dec) in degrees."""
    ra = math.radians(ra_deg % 360.0)
    dec = max(-90.0, min(90.0, dec_deg))
    z = math.sin(math.radians(dec))
    za = abs(z)
    tt = (ra % (2.0 * math.pi)) / (0.5 * math.pi)  # in [0, 4)
    if za <= 2.0 / 3.0:
        temp1 = nside * (0.5 + tt)
        temp2 = nside * z * 0.75
        jp = int(temp1 - temp2)
        jm = int(temp1 + temp2)
        ir = nside + 1 + jp - jm
        kshift = 1 - (ir & 1)
        ip = (jp + jm - nside + kshift + 1) // 2
        ip = ip % (4 * nside)
        return 2 * nside * (nside - 1) + (ir - 1) * 4 * nside + ip

    tp = tt - int(tt)
    tmp = nside * math.sqrt(3.0 * (1.0 - za))
    jp = int(tp * tmp)
    jm = int((1.0 - tp) * tmp)
    ir = jp + jm + 1
    ip = int(tt * ir) % (4 * ir)
    if z > 0:
        return 2 * ir * (ir - 1) + ip
    return 12 * nside * nside - 2 * ir * (ir + 1) + ip


def healpix_nest_to_radec(pix: int, nside: int = 4096) -> tuple[float, float]:
    """Decode HEALPix NESTED pixel index (used in Gaia source_id >> 35) to (ra_deg, dec_deg)."""
    nside2 = nside * nside
    face = pix // nside2
    ipf = pix % nside2
    ix = 0
    iy = 0
    for b in range(16):
        ix |= ((ipf >> (2 * b)) & 1) << b
        iy |= ((ipf >> (2 * b + 1)) & 1) << b
    jrll = (2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
    jpll = (1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7)
    jr = jrll[face] * nside - ix - iy - 1
    if jr < nside:
        z = 1.0 - (jr * jr) / (3.0 * nside2)
        nr = jr
        kshift = 0
    elif jr > 3 * nside:
        jr_s = 4 * nside - jr
        z = (jr_s * jr_s) / (3.0 * nside2) - 1.0
        nr = jr_s
        kshift = 0
    else:
        z = (2 * nside - jr) * 2.0 / (3.0 * nside)
        nr = nside
        kshift = (jr - nside) & 1
    jp = (jpll[face] * nr + ix - iy + 1 + kshift) // 2
    if jp > 4 * nside:
        jp -= 4 * nside
    if jp < 1:
        jp += 4 * nside
    phi = (jp - (kshift + 1) * 0.5) * (math.pi / (2.0 * nr))
    ra = math.degrees(phi) % 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, z))))
    return ra, dec


def equatorial_to_galactic(ra_deg: float, dec_deg: float) -> tuple[float, float]:
    """Convert ICRS (RA, Dec) in degrees to Galactic (l, b) in degrees."""
    ra = math.radians(ra_deg)
    dec = math.radians(dec_deg)
    # ICRS J2000 North Galactic Pole and ascending node
    ra_ngp = math.radians(192.85948)
    dec_ngp = math.radians(27.12825)
    l_ncp = math.radians(122.93192)
    sin_b = math.sin(dec) * math.sin(dec_ngp) + math.cos(dec) * math.cos(dec_ngp) * math.cos(ra - ra_ngp)
    b = math.asin(max(-1.0, min(1.0, sin_b)))
    y = math.cos(dec) * math.sin(ra - ra_ngp)
    x = math.sin(dec) * math.cos(dec_ngp) - math.cos(dec) * math.sin(dec_ngp) * math.cos(ra - ra_ngp)
    lon = (l_ncp - math.atan2(y, x)) % (2.0 * math.pi)
    return math.degrees(lon), math.degrees(b)


class TableStore:
    """HEALPix-partitioned table store backed by SQLite for fast SQL queries."""

    def __init__(self, tables_dir: Path, nside: int = 16) -> None:
        self.tables_dir = tables_dir
        self.nside = nside
        self.tables_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = self.tables_dir / "catalog.sqlite"

    def write_partitioned(self, table_name: str, rows: list[dict[str, Any]]) -> list[Path]:
        """Write rows partitioned by HEALPix cell and index into SQLite."""
        by_cell: dict[int, list[dict[str, Any]]] = {}
        for row in rows:
            ra = float(row.get("ra_deg", 0.0))
            dec = float(row.get("dec_deg", 0.0))
            cell = ra_dec_to_healpix(ra, dec, nside=self.nside)
            enriched = dict(row)
            enriched["healpix"] = cell
            by_cell.setdefault(cell, []).append(enriched)

        table_dir = self.tables_dir / table_name
        table_dir.mkdir(parents=True, exist_ok=True)
        written: list[Path] = []
        for cell, cell_rows in sorted(by_cell.items()):
            part_dir = table_dir / f"healpix={cell}"
            part_dir.mkdir(parents=True, exist_ok=True)
            part_file = part_dir / "rows.json"
            part_file.write_text(json.dumps(cell_rows, sort_keys=True), encoding="utf-8")
            written.append(part_file)

        with sqlite3.connect(str(self.db_path)) as conn:
            conn.execute(
                f"CREATE TABLE IF NOT EXISTS {table_name} ("
                "source_id TEXT PRIMARY KEY, "
                "ra_deg REAL, "
                "dec_deg REAL, "
                "healpix INTEGER, "
                "payload_json TEXT)"
            )
            conn.execute(f"DELETE FROM {table_name}")
            conn.executemany(
                f"INSERT OR REPLACE INTO {table_name} (source_id, ra_deg, dec_deg, healpix, payload_json) "
                "VALUES (?, ?, ?, ?, ?)",
                [
                    (
                        str(r["source_id"]),
                        float(r.get("ra_deg", 0.0)),
                        float(r.get("dec_deg", 0.0)),
                        int(ra_dec_to_healpix(float(r.get("ra_deg", 0.0)), float(r.get("dec_deg", 0.0)), self.nside)),
                        json.dumps(r, sort_keys=True),
                    )
                    for r in rows
                ],
            )
            conn.commit()
        return written

    def query(self, table_name: str, where_sql: str = "1=1", params: tuple[Any, ...] = ()) -> list[dict[str, Any]]:
        if not self.db_path.exists():
            return []
        with sqlite3.connect(str(self.db_path)) as conn:
            cur = conn.execute(
                f"SELECT payload_json FROM {table_name} WHERE {where_sql}",
                params,
            )
            return [json.loads(row[0]) for row in cur.fetchall()]
