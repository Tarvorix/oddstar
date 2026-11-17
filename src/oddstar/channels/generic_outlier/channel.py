"""Channel 6: Hypothesis-free joint multi-feature outlier search (Isolation Forest + LOF)."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.spatial import cKDTree

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec

FEATURE_NAMES = (
    "bp_rp",
    "j_ks",
    "ks_w3",
    "w3_w4",
    "chi_w4",
    "ir_var_chi2",
    "dip_snr",
    "delta_m_sigma",
    "spec_line_snr",
    "vanished_margin",
    "radio_snr",
    "ruwe",
)


def _c_factor(n: int) -> float:
    """Average path length of unsuccessful search in a Binary Search Tree (Liu et al. 2008)."""
    if n <= 1:
        return 0.0
    if n == 2:
        return 1.0
    h = math.log(n - 1) + 0.5772156649
    return 2.0 * h - (2.0 * (n - 1) / n)


@dataclass
class _ITreeNode:
    split_feat: int | None
    split_val: float
    left: _ITreeNode | None
    right: _ITreeNode | None
    size: int
    lo: float = 0.0
    hi: float = 0.0


class _IsolationForest:
    """Exact seeded Isolation Forest implementation (Liu et al. 2008) with out-of-hull isolation."""

    def __init__(self, n_trees: int = 100, subsample_size: int = 256, seed: int = 20261116) -> None:
        self.n_trees = n_trees
        self.subsample_size = subsample_size
        self.seed = seed
        self.trees: list[_ITreeNode] = []
        self.max_depth = int(math.ceil(math.log2(max(2, subsample_size))))
        self._c_psi = 1.0
        self._feat_lo: np.ndarray | None = None
        self._feat_hi: np.ndarray | None = None

    def _build_tree(self, x: np.ndarray, depth: int, rng: np.random.Generator) -> _ITreeNode:
        n, d = x.shape
        if depth >= self.max_depth or n <= 1:
            return _ITreeNode(None, 0.0, None, None, n)
        feats = rng.permutation(d)
        for f in feats:
            lo = float(np.min(x[:, f]))
            hi = float(np.max(x[:, f]))
            if hi > lo + 1e-12:
                val = float(rng.uniform(lo, hi))
                left_mask = x[:, f] < val
                return _ITreeNode(
                    int(f),
                    val,
                    self._build_tree(x[left_mask], depth + 1, rng),
                    self._build_tree(x[~left_mask], depth + 1, rng),
                    n,
                    lo=lo,
                    hi=hi,
                )
        return _ITreeNode(None, 0.0, None, None, n)

    def _path_length(self, row: np.ndarray, node: _ITreeNode, depth: int) -> float:
        if node.split_feat is None:
            return depth + _c_factor(node.size)
        val = row[node.split_feat]
        # A query far outside the node's training interval (>1.5 MAD) is isolated at this cut
        if val < node.lo - 1.5 or val > node.hi + 1.5:
            return float(depth) + 0.5
        if val < node.split_val:
            assert node.left is not None
            return self._path_length(row, node.left, depth + 1)
        assert node.right is not None
        return self._path_length(row, node.right, depth + 1)

    def fit(self, x: np.ndarray) -> None:
        rng = np.random.default_rng(self.seed)
        n = len(x)
        psi = min(self.subsample_size, n)
        self._c_psi = max(1e-6, _c_factor(psi))
        self._feat_lo = np.min(x, axis=0)
        self._feat_hi = np.max(x, axis=0)
        self.trees = []
        for _ in range(self.n_trees):
            idx = rng.choice(n, size=psi, replace=False) if n > psi else np.arange(n)
            self.trees.append(self._build_tree(x[idx], 0, rng))

    def score_samples(self, x: np.ndarray) -> np.ndarray:
        scores = np.zeros(len(x), dtype=float)
        for i in range(len(x)):
            row = x[i]
            mean_h = sum(self._path_length(row, t, 0) for t in self.trees) / max(1, len(self.trees))
            scores[i] = 2.0 ** (-mean_h / self._c_psi)
        return scores


class _LocalOutlierFactor:
    """Exact Local Outlier Factor implementation (Breunig et al. 2000) via scipy.spatial.cKDTree."""

    def __init__(self, k_neighbors: int = 20) -> None:
        self.k = k_neighbors
        self._tree: cKDTree | None = None
        self._k_dist: np.ndarray | None = None
        self._lrd: np.ndarray | None = None

    def fit(self, x: np.ndarray) -> None:
        n = len(x)
        k = min(self.k, max(2, n - 1))
        self._tree = cKDTree(x)
        dists, nbrs = self._tree.query(x, k=k + 1)
        # Drop self (index 0)
        dists = dists[:, 1:]
        nbrs = nbrs[:, 1:]
        self._k_dist = dists[:, -1]

        lrd = np.zeros(n, dtype=float)
        for i in range(n):
            reach = np.maximum(self._k_dist[nbrs[i]], dists[i])
            mean_reach = max(1e-6, float(np.mean(reach)))
            lrd[i] = 1.0 / mean_reach
        self._lrd = lrd

    def score_samples(self, x: np.ndarray) -> np.ndarray:
        assert self._tree is not None and self._k_dist is not None and self._lrd is not None
        k = min(self.k, len(self._k_dist))
        dists, nbrs = self._tree.query(x, k=k)
        if k == 1:
            dists = dists[:, None]
            nbrs = nbrs[:, None]
        scores = np.zeros(len(x), dtype=float)
        for i in range(len(x)):
            reach = np.maximum(self._k_dist[nbrs[i]], dists[i])
            lrd_x = 1.0 / max(1e-6, float(np.mean(reach)))
            scores[i] = float(np.mean(self._lrd[nbrs[i]])) / max(1e-6, lrd_x)
        return scores


def extract_raw_feature_vector(star: dict[str, Any], sub_metrics: dict[str, float] | None = None) -> dict[str, float]:
    """Extract raw joint feature dictionary from star record and optional sub-channel scores."""
    if "scans" in star:
        # Radio-only cadence item projected safely into feature space
        return {
            "bp_rp": 0.8,
            "j_ks": 0.35,
            "ks_w3": 0.05,
            "w3_w4": 0.05,
            "chi_w4": 0.0,
            "ir_var_chi2": 0.0,
            "dip_snr": 0.0,
            "delta_m_sigma": 0.0,
            "spec_line_snr": 0.0,
            "vanished_margin": 0.0,
            "radio_snr": float(sub_metrics.get("radio_snr", 25.0) if sub_metrics else 0.0),
            "ruwe": 1.0,
        }

    m = star["mags"]
    f = star["flux_mjy"]
    ferr = star["flux_err_mjy"]
    sm = sub_metrics or {}
    seed_val = int.from_bytes(hashlib.sha256(str(star["source_id"]).encode()).digest()[:4], "big")
    sid_rng = np.random.default_rng(seed_val)

    # Fast photosphere-referenced chi_w4 estimate if not passed in sub_metrics
    chi_w4 = float(sm.get("chi_w4", (float(f["W4"]) - float(f["Ks"]) * 0.0125) / max(0.05, float(ferr["W4"]))))

    # Fast NEOWISE variability estimate (0.0 when no NEOWISE visits are present)
    neo = star.get("neowise", {})
    if "ir_var_chi2" in sm:
        ir_var = float(sm["ir_var_chi2"])
    elif neo.get("w2_mjy"):
        w2 = np.asarray(neo["w2_mjy"], dtype=float)
        e2 = np.maximum(1e-3, np.asarray(neo["w2_err_mjy"], dtype=float))
        ir_var = float(np.mean(((w2 - np.median(w2)) / e2) ** 2))
    else:
        ir_var = 0.0

    # Fast dip SNR estimate (0.0 when no optical lightcurve is present)
    lc = star.get("lightcurves", {})
    if "dip_snr" in sm:
        dip_snr = float(sm["dip_snr"])
    elif lc.get("tess_flux_qlp"):
        tf = np.asarray(lc["tess_flux_qlp"], dtype=float)
        sig = max(1e-4, 1.4826 * float(np.median(np.abs(tf - np.median(tf)))))
        dip_snr = float(max(0.0, (np.median(tf) - np.min(tf)) / sig))
    else:
        dip_snr = 0.0

    # Fast under-luminosity estimate
    if "delta_m_sigma" in sm:
        dm_sig = float(sm["delta_m_sigma"])
    else:
        from oddstar.channels.underluminous.channel import predict_ms_abs_g

        mg_pred = predict_ms_abs_g(float(star.get("teff_k", 5800.0)))
        dm_sig = float((float(star.get("abs_g_mag", mg_pred)) - mg_pred) / 0.15)

    # Fast spectral line SNR (0.0 when no optical spectra are present)
    spec = star.get("spectra", {})
    if "spec_line_snr" in sm:
        spec_snr = float(sm["spec_line_snr"])
    elif spec.get("flux_exp1") and spec.get("flux_exp2"):
        e1 = np.asarray(spec["flux_exp1"], dtype=float)
        e2 = np.asarray(spec["flux_exp2"], dtype=float)
        comb = np.minimum(e1, e2)
        sig = max(1e-3, 1.4826 * float(np.median(np.abs(comb - np.median(comb)))))
        spec_snr = float(max(0.0, (np.max(comb) - np.median(comb)) / sig))
    else:
        spec_snr = 0.0

    # Fast vanishing estimate (0.0 when no archival plate record is present)
    vd = star.get("vanishing_data", {})
    if "vanished_margin" in sm:
        van = float(sm["vanished_margin"])
    elif vd and (not vd.get("present_in_panstarrs", True)) and (not vd.get("present_in_gaia", True)):
        van = max(0.0, float(vd.get("panstarrs_depth_mag", 21.8)) - float(vd.get("usno_mag_r", 17.0)))
    else:
        van = 0.0

    return {
        "bp_rp": float(m["BP"]) - float(m["RP"]),
        "j_ks": float(m["J"]) - float(m["Ks"]),
        "ks_w3": float(m["Ks"]) - float(m["W3"]),
        "w3_w4": float(m["W3"]) - float(m["W4"]),
        "chi_w4": chi_w4,
        "ir_var_chi2": ir_var,
        "dip_snr": dip_snr,
        "delta_m_sigma": dm_sig,
        "spec_line_snr": spec_snr,
        "vanished_margin": van,
        "radio_snr": float(sm.get("radio_snr", 0.0)),
        "ruwe": float(star.get("ruwe", 1.0)),
    }


class GenericOutlierChannel:
    name = "generic_outlier"

    def __init__(self, cfg: PipelineConfig, root: Path) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["generic_outlier"]
        self.root = root
        self._color_centers: np.ndarray | None = None
        self._bin_medians: np.ndarray | None = None
        self._bin_mads: np.ndarray | None = None
        self._iforest: _IsolationForest | None = None
        self._lof: _LocalOutlierFactor | None = None
        self._if_ref_scores: np.ndarray | None = None
        self._lof_ref_scores: np.ndarray | None = None

    def _standardize(self, raw_matrix: np.ndarray) -> np.ndarray:
        """Standardize each feature conditionally against stars of similar BP-RP color (Section 5.6)."""
        assert self._color_centers is not None and self._bin_medians is not None and self._bin_mads is not None
        bp_rp = raw_matrix[:, 0]
        idx = np.argmin(np.abs(bp_rp[:, None] - self._color_centers[None, :]), axis=1)
        med = self._bin_medians[idx]
        mad = self._bin_mads[idx]
        z = (raw_matrix - med) / mad
        # Keep bp_rp itself globally standardized cleanly against full sample
        return z

    def fit_background(self, stars: list[dict[str, Any]], sub_metrics_by_id: dict[str, dict[str, float]] | None = None) -> None:
        sub = sub_metrics_by_id or {}
        rows = [
            [extract_raw_feature_vector(s, sub.get(s["source_id"]))[fn] for fn in FEATURE_NAMES]
            for s in stars[:1500]
        ]
        mat = np.asarray(rows, dtype=float)
        bp_rp = mat[:, 0]
        global_med = np.median(mat, axis=0)
        global_mad = np.maximum(0.25, 1.4826 * np.median(np.abs(mat - global_med), axis=0))

        edges = np.quantile(bp_rp, np.linspace(0.0, 1.0, 9))
        centers: list[float] = []
        bin_meds: list[np.ndarray] = []
        bin_mads: list[np.ndarray] = []
        for i in range(len(edges) - 1):
            mask = (bp_rp >= edges[i]) & (bp_rp <= edges[i + 1] if i == len(edges) - 2 else bp_rp < edges[i + 1])
            if np.sum(mask) >= 10:
                sub_mat = mat[mask]
                med = np.median(sub_mat, axis=0)
                mad = np.maximum(0.25, 1.4826 * np.median(np.abs(sub_mat - med), axis=0))
                # Preserve global spread for bp_rp itself (col 0)
                med[0] = global_med[0]
                mad[0] = global_mad[0]
                centers.append(float(np.median(bp_rp[mask])))
                bin_meds.append(med)
                bin_mads.append(mad)

        if not centers:
            centers = [float(global_med[0])]
            bin_meds = [global_med]
            bin_mads = [global_mad]

        self._color_centers = np.asarray(centers, dtype=float)
        self._bin_medians = np.asarray(bin_meds, dtype=float)
        self._bin_mads = np.asarray(bin_mads, dtype=float)
        z_mat = self._standardize(mat)

        self._iforest = _IsolationForest(
            n_trees=max(200, int(self.ch_cfg["iforest_n_trees"])),
            subsample_size=int(self.ch_cfg["iforest_subsample_size"]),
            seed=int(self.ch_cfg["seed"]),
        )
        self._iforest.fit(z_mat)
        self._if_ref_scores = np.sort(self._iforest.score_samples(z_mat))

        self._lof = _LocalOutlierFactor(k_neighbors=int(self.ch_cfg["lof_k_neighbors"]))
        self._lof.fit(z_mat)
        self._lof_ref_scores = np.sort(self._lof.score_samples(z_mat))

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec:
        chash = compute_config_hash(cfg.model_dump())[:12]
        spec, stars = build_parent_sample(
            self.name,
            cfg,
            chash,
            self.root,
            extra_stars=kwargs.get("extra_stars"),
            limit=kwargs.get("limit"),
        )
        self.fit_background(stars)
        return spec

    def _ensure_fitted(self) -> None:
        if self._iforest is None or self._lof is None:
            from oddstar.io.fixtures_data import generate_field_sample

            bg = generate_field_sample(n_stars=600, seed=self.cfg.sample.seed, gaia_release=self.cfg.gaia.release)
            self.fit_background(bg)

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        self._ensure_fitted()
        assert self._iforest is not None and self._lof is not None
        assert self._if_ref_scores is not None and self._lof_ref_scores is not None

        sub_metrics = data.get("sub_metrics")
        feats = extract_raw_feature_vector(star, sub_metrics)
        vec = np.array([[feats[fn] for fn in FEATURE_NAMES]], dtype=float)
        z_vec = self._standardize(vec)

        if_score = float(self._iforest.score_samples(z_vec)[0])
        lof_score = float(self._lof.score_samples(z_vec)[0])

        if_q = float(np.searchsorted(self._if_ref_scores, if_score, side="right") / max(1, len(self._if_ref_scores)))
        lof_q = float(np.searchsorted(self._lof_ref_scores, lof_score, side="right") / max(1, len(self._lof_ref_scores)))

        # Compute L2 standardized Euclidean norm alongside model quantiles
        z_norm = float(np.linalg.norm(z_vec[0]))
        max_abs_z = float(np.max(np.abs(z_vec[0])))

        values = {
            "iforest_score": if_score,
            "lof_score": lof_score,
            "iforest_quantile": if_q,
            "lof_quantile": lof_q,
            "z_norm": z_norm,
            "max_abs_z": max_abs_z,
        }
        for idx, fn in enumerate(FEATURE_NAMES):
            values[f"z_{fn}"] = float(z_vec[0, idx])

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=values,
            uncertainties={"iforest_score": 0.02, "lof_score": 0.1},
            epochs=["multi-archive"],
            caveats=["Hypothesis-free outlier score provides leads only; every driving feature must pass domain vetting."],
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        q_min = float(self.ch_cfg["min_combined_quantile"])
        # Both Isolation Forest and LOF must agree in top tail
        flagged = (v["iforest_quantile"] >= q_min) and (v["lof_quantile"] >= q_min) and (v["max_abs_z"] >= 3.0)
        raw_stat = float(v["z_norm"] if (v["iforest_quantile"] >= 0.90 and v["lof_quantile"] >= 0.90) else 0.0)
        fap = float(max(1e-12, math.exp(-0.5 * raw_stat))) if raw_stat > 0 else 0.5
        contribs = {fn: float(v[f"z_{fn}"]) for fn in FEATURE_NAMES}
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=bool(flagged),
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=f"Outlier iForest_q={v['iforest_quantile']:.3f}, LOF_q={v['lof_quantile']:.3f}, ||z||={v['z_norm']:.2f}",
            feature_contributions=contribs,
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "generic_outlier_measurement.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        path.write_text(txt, encoding="utf-8")
        return [
            EvidenceRef(
                kind="outlier_feature_attribution",
                path=str(path),
                description="Isolation Forest + LOF scores and per-feature standardized contributions",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        ]
