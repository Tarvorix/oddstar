"""Parent sample builder with cut-by-cut accounting and HEALPix table storage."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.io.fixtures_data import generate_field_sample, get_named_object
from oddstar.io.storage import TableStore, ra_dec_to_healpix
from oddstar.schemas import SampleCutStep, SampleSpec


def build_parent_sample(
    channel: str,
    cfg: PipelineConfig,
    config_hash: str,
    root: Path,
    *,
    extra_stars: list[dict[str, Any]] | None = None,
    limit: int | None = None,
    apply_cuts: bool = True,
) -> tuple[SampleSpec, list[dict[str, Any]]]:
    """Construct parent sample for a channel, applying configured cuts and logging every step."""
    n_base = limit if limit is not None else cfg.sample.field_sample_size
    stars = generate_field_sample(
        n_stars=n_base,
        seed=cfg.sample.seed,
        gaia_release=cfg.gaia.release,
    )
    if extra_stars:
        by_id = {s["source_id"]: s for s in stars}
        for es in extra_stars:
            by_id[es["source_id"]] = es
        stars = list(by_id.values())

    cuts: list[SampleCutStep] = []
    cur = list(stars)

    if apply_cuts:
        # Cut 1: Parallax SNR
        before = len(cur)
        cur = [s for s in cur if float(s.get("parallax_snr", 0.0)) > cfg.sample.min_parallax_snr]
        cuts.append(
            SampleCutStep(
                cut_name="parallax_snr",
                description=f"parallax_snr > {cfg.sample.min_parallax_snr}",
                count_before=before,
                count_after=len(cur),
            )
        )

        # Cut 2: Distance limit
        before = len(cur)
        cur = [s for s in cur if float(s.get("distance_pc", 1e9)) <= cfg.sample.max_distance_pc]
        cuts.append(
            SampleCutStep(
                cut_name="distance_pc",
                description=f"distance_pc <= {cfg.sample.max_distance_pc}",
                count_before=before,
                count_after=len(cur),
            )
        )

        # Cut 3: Astrometric RUWE
        before = len(cur)
        cur = [s for s in cur if float(s.get("ruwe", 99.0)) < cfg.sample.max_ruwe]
        cuts.append(
            SampleCutStep(
                cut_name="ruwe",
                description=f"ruwe < {cfg.sample.max_ruwe}",
                count_before=before,
                count_after=len(cur),
            )
        )

        # Cut 4: Galactic plane exclusion
        before = len(cur)
        cur = [s for s in cur if abs(float(s.get("gal_b_deg", 0.0))) > cfg.sample.min_abs_gal_lat_deg]
        cuts.append(
            SampleCutStep(
                cut_name="galactic_latitude",
                description=f"|b| > {cfg.sample.min_abs_gal_lat_deg} deg",
                count_before=before,
                count_after=len(cur),
            )
        )

        # Cut 5: Unique AllWISE cross-match
        if cfg.sample.require_unique_allwise_xmatch:
            before = len(cur)
            cur = [
                s
                for s in cur
                if int(s.get("allwise", {}).get("n_neighbours_gaia_xmatch", 1)) == 1
            ]
            cuts.append(
                SampleCutStep(
                    cut_name="unique_allwise_xmatch",
                    description="Unique Gaia-AllWISE cross-match",
                    count_before=before,
                    count_after=len(cur),
                )
            )

        # Cut 6: Clean 2MASS photometry
        if cfg.sample.require_clean_2mass:
            before = len(cur)
            cur = [s for s in cur if bool(s.get("2mass_clean", True))]
            cuts.append(
                SampleCutStep(
                    cut_name="clean_2mass",
                    description="Clean 2MASS JHKs photometry flags",
                    count_before=before,
                    count_after=len(cur),
                )
            )

    nside = cfg.storage.healpix_nside
    table_dir = root / cfg.storage.tables_dir
    store = TableStore(table_dir, nside=nside)
    table_name = f"sample_{channel}_{cfg.gaia.release}"
    store.write_partitioned(table_name, cur)

    healpix_cells = sorted({ra_dec_to_healpix(float(s["ra_deg"]), float(s["dec_deg"]), nside=nside) for s in cur})
    spec = SampleSpec(
        channel=channel,
        gaia_release=cfg.gaia.release,
        config_hash=config_hash,
        seed=cfg.sample.seed,
        cuts=cuts,
        total_count=len(cur),
        healpix_cells=healpix_cells,
        source_ids=[str(s["source_id"]) for s in cur],
        table_path=str((table_dir / table_name).relative_to(root)),
    )
    return spec, cur


def resolve_star_record(source_id: str, cfg: PipelineConfig) -> dict[str, Any]:
    """Load a named star or generate deterministic archival record for one-star inspection."""
    clean = source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
    star = get_named_object(clean, gaia_release=cfg.gaia.release)
    if star is not None:
        return star
    from oddstar.io.fixtures_data import _make_base_star

    return _make_base_star(clean, gaia_release=cfg.gaia.release)
