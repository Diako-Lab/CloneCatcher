"""CloneCatcher Module 6 -- ClusterCatcher integration (single_cell mode).

Adds receptor identity to the ClusterCatcher AnnData, producing adata_tcr.h5ad
with clonotype, chain pairing, expansion status, and cluster membership
alongside whatever ClusterCatcher already holds: cell type, mutational
signature exposure, copy number, viral reads.

This is the object the paper's single-cell claims are made from, because it is
the only place a T cell's receptor sits next to the tumour cell state it was
responding to.

The join reads ClusterCatcher at run time rather than copying its columns into
the flat clonotype tables, so a ClusterCatcher rerun does not leave a stale
duplicate behind. The NNLS refit that changes signature weights will change
this object the next time it is built, and will not silently disagree with a
TSV written weeks earlier.

Output: integration/adata_tcr.h5ad
"""

import logging
import sys
from pathlib import Path

import pandas as pd

import barcode_utils as bu

# Columns added to adata.obs. Everything is prefixed so nothing collides with
# ClusterCatcher's own obs, whatever it holds.
PREFIX = "tcr_"


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params

    import anndata as ad

    adata_path = bu.find_clustercatcher_adata(p.clustercatcher_run_dir)
    logging.info("Reading %s", adata_path)
    adata = ad.read_h5ad(adata_path)
    adata.obs_names = adata.obs_names.astype(str)
    logging.info("  %d cells, %d obs columns", adata.n_obs, adata.obs.shape[1])

    bu.assert_canonical(adata.obs_names, source=str(adata_path))

    by_cell = pd.read_csv(snakemake.input.by_cell, sep="\t", dtype=str).fillna("")
    by_sample = pd.read_csv(snakemake.input.by_sample, sep="\t")
    logging.info(
        "  %d TCR cells, %d clonotypes", len(by_cell), by_sample["clonotype_id"].nunique()
    )

    bu.assert_canonical(by_cell["cell_barcode"], source="clonotypes_by_cell.tsv")

    # The tripwire again. These two came from the same run and must overlap.
    overlap, rate = bu.assert_join_rate(
        by_cell["cell_barcode"],
        adata.obs_names,
        min_rate=float(p.min_join_rate),
        label="TCR cells -> AnnData",
    )
    logging.info("TCR cells matched into AnnData: %d (%.1f%%)", overlap, 100 * rate)

    # Per-clonotype attributes, looked up by the shared clonotype_id.
    clono = by_sample.set_index("clonotype_id")
    cell = by_cell.set_index("cell_barcode")

    frame = pd.DataFrame(index=adata.obs_names)
    frame[f"{PREFIX}clonotype_id"] = cell["clonotype_id"].reindex(frame.index)
    frame[f"{PREFIX}cellranger_clonotype_id"] = cell.get(
        "cellranger_clonotype_id", pd.Series(dtype=str)
    ).reindex(frame.index)
    frame[f"{PREFIX}n_chains"] = pd.to_numeric(
        cell["n_chains"].reindex(frame.index), errors="coerce"
    )
    frame[f"{PREFIX}is_paired"] = (
        cell["is_paired"].reindex(frame.index).astype(str).str.lower().isin(("true", "1"))
    )

    ids = frame[f"{PREFIX}clonotype_id"]
    for src, dest in (
        ("frequency", f"{PREFIX}clone_frequency"),
        ("n_observations", f"{PREFIX}clone_size"),
        ("is_expanded", f"{PREFIX}is_expanded"),
        ("cdr3a_aa", f"{PREFIX}cdr3a_aa"),
        ("cdr3b_aa", f"{PREFIX}cdr3b_aa"),
    ):
        if src in clono.columns:
            frame[dest] = ids.map(clono[src])

    # Optional: specificity clusters, when that stage ran.
    clusters_path = Path(snakemake.input.by_sample).parents[1] / "specificity" / "specificity_clusters.tsv"
    if clusters_path.exists():
        clusters = pd.read_csv(clusters_path, sep="\t").set_index("clone_id")
        frame[f"{PREFIX}cluster_id"] = ids.map(clusters["cluster_id"])
        if "cluster_n_samples" in clusters.columns:
            frame[f"{PREFIX}cluster_n_samples"] = ids.map(clusters["cluster_n_samples"])
        logging.info("Added specificity cluster membership from %s", clusters_path)

    # has_tcr is the column most downstream work will filter on, so make it
    # explicit rather than leaving people to test clonotype_id for null.
    frame[f"{PREFIX}has_tcr"] = frame[f"{PREFIX}clonotype_id"].notna()

    for col in frame.columns:
        adata.obs[col] = frame[col].values

    n_tcr = int(adata.obs[f"{PREFIX}has_tcr"].sum())
    logging.info("Cells with a receptor: %d of %d (%.1f%%)",
                 n_tcr, adata.n_obs, 100 * n_tcr / max(adata.n_obs, 1))

    if "final_annotation" in adata.obs.columns:
        tab = (
            adata.obs.loc[adata.obs[f"{PREFIX}has_tcr"], "final_annotation"]
            .value_counts()
            .head(10)
        )
        logging.info("Receptor-bearing cells by annotated type:")
        for label, n in tab.items():
            logging.info("  %-28s %6d", str(label)[:28], n)
        # A sanity read, not an assertion: most receptors should land on cells
        # annotated as T cells. A large fraction elsewhere suggests either
        # ambient contamination or an annotation problem, and is worth a look
        # before anything is concluded from this object.
        t_like = [l for l in tab.index if "T cell" in str(l) or str(l).startswith("T ")]
        if t_like:
            frac = tab[t_like].sum() / max(tab.sum(), 1)
            logging.info("  %.0f%% of the top types are T-cell labelled", 100 * frac)

    out = Path(snakemake.output.adata)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    logging.info("Wrote %s", out)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
