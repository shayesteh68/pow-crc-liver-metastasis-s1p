from paths import DATA_DIR, RESULTS_DIR, FIG_DIR, PLOTS_DIR

import argparse
import gc
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib.pyplot as plt
from scipy import sparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--allow-heuristic-doublet-fallback",
    action="store_true",
    help="Allow a labeled gene-count heuristic if Scrublet fails or cannot set a threshold.",
)
parser.add_argument(
    "--thresholds-tsv",
    help="Optional reviewed TSV with patient, tissue, and threshold columns; applies and filters only these explicit thresholds.",
)
args = parser.parse_args()

IN_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_raw.h5ad")
OUT_FP = os.path.join(RESULTS_DIR, "GSE178318_CRC_LM_qc.h5ad")

print(f"[*] Loading raw CRC+LM object: {IN_FP}")
adata = sc.read_h5ad(IN_FP)
required_obs = {"patient", "tissue"}
missing_obs = required_obs - set(adata.obs.columns)
if missing_obs:
    raise KeyError("The AnnData object must contain patient and tissue columns for sample-level doublet detection.")
if adata.obs[["patient", "tissue"]].isna().any().any():
    raise ValueError("Missing patient or tissue labels; refusing to pool distinct libraries for doublet detection.")

if "gene_symbol" in adata.var.columns:
    mito_mask = adata.var["gene_symbol"].astype(str).str.upper().str.startswith("MT-")
else:
    mito_mask = np.array([False] * adata.n_vars)
adata.var["mt"] = mito_mask

print("[*] Calculating basic QC metrics...")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

valid_cells = (
    (adata.obs["n_genes_by_counts"] >= 200)
    & (adata.obs["n_genes_by_counts"] <= 6000)
    & (adata.obs["pct_counts_mt"] <= 20.0)
)
before_qc = adata.n_obs
adata = adata[valid_cells].copy()
print(f"[+] Retained after QC: {adata.n_obs:,} cells (filtered: {before_qc - adata.n_obs:,})")
gc.collect()

doublet_scores = pd.Series(np.nan, index=adata.obs_names, dtype=float)
scrublet_candidates = pd.Series(np.nan, index=adata.obs_names, dtype=float)
predicted_doublets = pd.Series(np.nan, index=adata.obs_names, dtype=float)
doublet_method = pd.Series("scrublet_auto_review_required", index=adata.obs_names, dtype=object)
doublet_thresholds = pd.Series(np.nan, index=adata.obs_names, dtype=float)
doublet_threshold_status = pd.Series("not_evaluated", index=adata.obs_names, dtype=object)
sample_groups = adata.obs["patient"].astype(str) + " | " + adata.obs["tissue"].astype(str)
adata.obs["scrublet_sample_id"] = sample_groups

manual_thresholds = None
if args.thresholds_tsv:
    thresholds_df = pd.read_csv(args.thresholds_tsv, sep="\t", dtype={"patient": str, "tissue": str})
    required_threshold_columns = {"patient", "tissue", "threshold"}
    missing_columns = required_threshold_columns - set(thresholds_df.columns)
    if missing_columns:
        raise ValueError("Threshold TSV is missing required columns: " + ", ".join(sorted(missing_columns)))
    if thresholds_df.duplicated(["patient", "tissue"]).any():
        raise ValueError("Threshold TSV must contain one unique row per patient/tissue library.")
    thresholds_df["threshold"] = pd.to_numeric(thresholds_df["threshold"], errors="coerce")
    if thresholds_df["threshold"].isna().any() or (~np.isfinite(thresholds_df["threshold"])).any():
        raise ValueError("All reviewed thresholds must be finite numeric values.")
    if ((thresholds_df["threshold"] < 0) | (thresholds_df["threshold"] > 1)).any():
        raise ValueError("Reviewed Scrublet thresholds must be between 0 and 1.")
    manual_thresholds = {
        f"{row.patient} | {row.tissue}": float(row.threshold)
        for row in thresholds_df.itertuples(index=False)
    }
    observed_sample_ids = set(sample_groups.unique())
    if set(manual_thresholds) != observed_sample_ids:
        missing = sorted(observed_sample_ids - set(manual_thresholds))
        extra = sorted(set(manual_thresholds) - observed_sample_ids)
        raise ValueError(
            "Threshold TSV must match the observed libraries exactly. Missing: "
            + ", ".join(missing)
            + "; unexpected: "
            + ", ".join(extra)
        )

print("[*] Running doublet detection separately for each patient/tissue library...")
for sample_id in sorted(sample_groups.unique()):
    patient, tissue = sample_id.split(" | ", maxsplit=1)
    sample_mask = sample_groups == sample_id
    print(f"  -> {sample_id}: {int(sample_mask.sum()):,} cells")
    sub_x = adata[sample_mask].X
    sub_x = sub_x.tocsr() if sparse.issparse(sub_x) else sparse.csr_matrix(sub_x)
    scrub = None
    scores = None
    calls = None
    effective_threshold = np.nan
    threshold_status = "scrublet_failed"
    try:
        scrub = scr.Scrublet(sub_x, expected_doublet_rate=0.05, random_state=0)
        scores, calls = scrub.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            verbose=False,
        )
        doublet_scores.loc[sample_mask] = scores
        if manual_thresholds is not None:
            effective_threshold = manual_thresholds[sample_id]
            calls = scrub.call_doublets(threshold=effective_threshold, verbose=False)
            scrublet_candidates.loc[sample_mask] = calls.astype(float)
            predicted_doublets.loc[sample_mask] = calls.astype(float)
            doublet_method.loc[sample_mask] = "scrublet_manual_threshold"
            threshold_status = "reviewed_manual_threshold_applied"
        else:
            if scrub.threshold_ is not None:
                effective_threshold = float(scrub.threshold_)
            if calls is not None:
                scrublet_candidates.loc[sample_mask] = np.asarray(calls, dtype=float)
            simulated_scores = np.asarray(getattr(scrub, "doublet_scores_sim_", []), dtype=float)
            if not np.isfinite(effective_threshold):
                threshold_status = "automatic_threshold_unavailable"
            elif simulated_scores.size == 0 or effective_threshold >= np.nanmax(simulated_scores):
                threshold_status = "automatic_threshold_outside_simulated_range"
            else:
                threshold_status = "automatic_threshold_requires_manual_review"
        doublet_thresholds.loc[sample_mask] = effective_threshold
        doublet_threshold_status.loc[sample_mask] = threshold_status
    except Exception as error:
        if manual_thresholds is not None:
            raise RuntimeError(
                f"Scrublet failed for sample {sample_id!r}; the reviewed threshold cannot be applied without Scrublet scores."
            ) from error
        if not args.allow_heuristic_doublet_fallback:
            raise RuntimeError(
                f"Scrublet failed for sample {sample_id!r}. To permit the labeled gene-count fallback, rerun with --allow-heuristic-doublet-fallback."
            ) from error
        print(f"     [!] Using explicitly enabled heuristic fallback for {sample_id}: {error}")
        high_counts = (
            adata.obs.loc[sample_mask, "n_genes_by_counts"]
            > adata.obs.loc[sample_mask, "n_genes_by_counts"].quantile(0.95)
        )
        scrublet_candidates.loc[sample_mask] = high_counts.to_numpy(dtype=float)
        doublet_method.loc[sample_mask] = "gene_count_top_5_percent"
        doublet_threshold_status.loc[sample_mask] = "heuristic_fallback_requires_review"

    if scores is not None:
        score_quantiles = np.nanquantile(scores, [0.5, 0.9, 0.95, 0.99, 1.0])
        sample_calls = scrublet_candidates.loc[sample_mask]
        sample_candidate_count = int(sample_calls.fillna(0).sum())
        sample_unknown_count = int(sample_calls.isna().sum())
        threshold_text = f"{effective_threshold:.4f}" if np.isfinite(effective_threshold) else "unavailable"
        detected_rate = getattr(scrub, "detected_doublet_rate_", np.nan)
        detectable_fraction = getattr(scrub, "detectable_doublet_fraction_", np.nan)
        estimated_rate = getattr(scrub, "overall_doublet_rate_", np.nan)

        def format_metric(value):
            try:
                return f"{float(value):.4f}"
            except (TypeError, ValueError):
                return "unavailable"

        print(
            f"     threshold={threshold_text} ({threshold_status}); "
            f"candidates={sample_candidate_count}/{len(sample_calls)}; unknown={sample_unknown_count}; "
            f"detected_rate={format_metric(detected_rate)}; "
            f"detectable_fraction={format_metric(detectable_fraction)}; "
            f"estimated_overall_rate={format_metric(estimated_rate)}; "
            "score quantiles (50/90/95/99/100%)="
            + ", ".join(f"{value:.4f}" for value in score_quantiles)
        )

    if scrub is not None and hasattr(scrub, "doublet_scores_sim_"):
        os.makedirs(PLOTS_DIR, exist_ok=True)
        bins = np.linspace(0, 1, 51)
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
        axes[0].hist(scrub.doublet_scores_obs_, bins=bins, color="steelblue")
        axes[1].hist(scrub.doublet_scores_sim_, bins=bins, color="darkorange")
        if np.isfinite(effective_threshold):
            for ax in axes:
                ax.axvline(effective_threshold, color="black", linestyle="--")
        axes[0].set(title="Observed cells", xlabel="Doublet score", ylabel="Cell count")
        axes[1].set(title="Simulated doublets", xlabel="Doublet score", ylabel="Simulation count")
        fig.suptitle(f"Scrublet threshold review: {sample_id} ({threshold_status})")
        fig.tight_layout()
        fig.savefig(os.path.join(PLOTS_DIR, f"scrublet_{patient}_{tissue}.png"), dpi=150)
        plt.close(fig)

    del sub_x
    gc.collect()

adata.obs["doublet_score"] = doublet_scores
adata.obs["scrublet_candidate"] = scrublet_candidates
adata.obs["predicted_doublet"] = predicted_doublets
adata.obs["doublet_method"] = doublet_method
adata.obs["doublet_threshold"] = doublet_thresholds
adata.obs["doublet_threshold_status"] = doublet_threshold_status
print("[+] Doublet methods:")
print(adata.obs["doublet_method"].value_counts(dropna=False))

if manual_thresholds is not None:
    n_doublets = int(predicted_doublets.sum())
    print(f"[+] Doublets under reviewed thresholds: {n_doublets:,} ({n_doublets / max(adata.n_obs, 1) * 100:.2f}%)")
    adata.uns["doublet_filter_status"] = "applied_reviewed_per_library_thresholds"
    keep_cells = ~adata.obs["predicted_doublet"].astype(bool)
    adata = adata[keep_cells].copy()
    print(f"[+] Cells retained after reviewed doublet filtering: {adata.n_obs:,}")
else:
    n_candidates = int(scrublet_candidates.fillna(0).sum())
    n_unknown = int(scrublet_candidates.isna().sum())
    print(
        f"[+] Scrublet candidates: {n_candidates:,}; unscored/unknown: {n_unknown:,}. "
        "No doublet calls were applied; provide a reviewed --thresholds-tsv to filter cells."
    )
    adata.uns["doublet_filter_status"] = "not_applied_auto_threshold_requires_manual_review"

print(f"[*] Saving to: {OUT_FP}")
adata.write_h5ad(OUT_FP, compression="gzip")
print("[+] Step 05 complete.")
