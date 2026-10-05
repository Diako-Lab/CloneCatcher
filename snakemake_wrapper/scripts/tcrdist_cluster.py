"""CloneCatcher Module 5 -- Specificity clustering (tcrdist3).

Groups receptors by shared specificity features and annotates the resulting
clusters against reference databases.

WHAT THIS IS AND IS NOT. tcrdist is a sequence-similarity metric weighted by
the CDR loops that contact peptide-MHC. Receptors close in that space often
share specificity, which makes this a prioritisation: it ranks candidates for
experimental testing. It does not demonstrate binding, and a cluster with no
database hit is uninformative rather than negative, since the databases cover
a small and biased slice of epitope space.

Runs in its own conda environment, because tcrdist3 does not support pandas
2.x and the rest of the pipeline pins 2.2.3 to match ClusterCatcher. Inputs
and outputs are both TSV, so nothing but tabular data crosses that boundary.

Outputs: specificity/specificity_clusters.tsv, specificity/database_hits.tsv,
         specificity/figures/
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

DB_FILES = {"vdjdb": "vdjdb.tsv", "mcpas": "mcpas.tsv"}


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def prepare(df):
    """Shape the shared schema into tcrdist3's expected frame.

    Paired when both chains resolve, beta-only otherwise. Bulk data is always
    beta-only, since the assay cannot pair chains; single-cell is usually
    paired. Paired distances are substantially more informative, so the two
    modes are reported with that difference made explicit rather than silently
    mixed.
    """
    df = df.copy()
    for col in ("cdr3a_aa", "v_a", "j_a", "cdr3b_aa", "v_b", "j_b"):
        df[col] = df.get(col, "").fillna("").astype(str)

    has_b = df["cdr3b_aa"].str.len() > 0
    has_a = df["cdr3a_aa"].str.len() > 0
    df = df[has_b]  # beta is required either way
    if df.empty:
        raise ValueError("No clonotype has a beta chain; tcrdist needs TRB.")

    paired = bool((has_a & has_b).mean() > 0.5)
    logging.info(
        "%d clonotypes with TRB; %.0f%% also have TRA -> running %s",
        len(df),
        100 * has_a[has_b].mean() if has_b.any() else 0,
        "paired alpha-beta" if paired else "beta only",
    )

    # tcrdist3 requires IMGT allele names; CellRanger and TRUST4 both emit
    # gene names without the allele. Appending *01 is the conventional fix.
    def allele(g):
        g = str(g)
        return g if not g or "*" in g else f"{g}*01"

    out = pd.DataFrame(
        {
            "clone_id": df["clonotype_id"].astype(str),
            "subject": df["sample_id"].astype(str),
            "count": df["n_observations"].astype(int),
            "cdr3_b_aa": df["cdr3b_aa"],
            "v_b_gene": df["v_b"].map(allele),
            "j_b_gene": df["j_b"].map(allele),
        }
    )
    if paired:
        out["cdr3_a_aa"] = df["cdr3a_aa"]
        out["v_a_gene"] = df["v_a"].map(allele)
        out["j_a_gene"] = df["j_a"].map(allele)
        out = out[
            (out["cdr3_a_aa"].str.len() > 0)
            & (out["v_a_gene"].str.len() > 0)
            & (out["j_a_gene"].str.len() > 0)
        ]
    out = out[(out["v_b_gene"].str.len() > 0) & (out["j_b_gene"].str.len() > 0)]
    return out.reset_index(drop=True), paired


def cluster(dist, radius):
    """Single-linkage components at a fixed tcrdist radius.

    Deliberately simple. The radius is the parameter that matters and it is
    exposed in the config; adding a clustering algorithm with its own
    hyperparameters on top would obscure rather than improve that.
    """
    n = dist.shape[0]
    adj = dist <= radius
    labels = -np.ones(n, dtype=int)
    current = 0
    for i in range(n):
        if labels[i] != -1:
            continue
        stack, members = [i], []
        while stack:
            j = stack.pop()
            if labels[j] != -1:
                continue
            labels[j] = current
            members.append(j)
            stack.extend(np.flatnonzero(adj[j] & (labels == -1)))
        current += 1
    return labels


def load_db(resources_dir, name):
    path = Path(resources_dir) / DB_FILES[name]
    if not path.exists():
        logging.warning(
            "%s not found at %s; skipping. Fetch with "
            "`CloneCatcher fetch-databases --output-dir %s`",
            name, path, resources_dir,
        )
        return None
    db = pd.read_csv(path, sep="\t", low_memory=False)
    cols = {c.lower(): c for c in db.columns}
    cdr3 = cols.get("cdr3") or cols.get("cdr3.beta.aa") or cols.get("cdr3b")
    epitope = cols.get("epitope") or cols.get("epitope.peptide")
    if not cdr3 or not epitope:
        logging.warning("%s: could not identify CDR3/epitope columns, skipping", name)
        return None
    out = db[[cdr3, epitope]].copy()
    out.columns = ["cdr3b_aa", "epitope"]
    species = cols.get("pathology") or cols.get("antigen.species")
    out["source"] = db[species] if species else name
    out["database"] = name
    return out.dropna(subset=["cdr3b_aa"]).drop_duplicates()


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params
    radius = int(p.radius)

    df = pd.read_csv(snakemake.input.by_sample, sep="\t")
    tcrs, paired = prepare(df)

    from tcrdist.repertoire import TCRrep

    chains = ["alpha", "beta"] if paired else ["beta"]
    logging.info("Computing tcrdist over %d receptors, chains=%s", len(tcrs), chains)

    rep = TCRrep(
        cell_df=tcrs,
        organism="human",
        chains=chains,
        deduplicate=False,
        compute_distances=True,
    )
    dist = rep.pw_alpha + rep.pw_beta if paired else rep.pw_beta

    labels = cluster(np.asarray(dist), radius)
    tcrs["cluster_id"] = [f"c{l:04d}" for l in labels]
    tcrs["chains"] = "alpha_beta" if paired else "beta"

    sizes = tcrs["cluster_id"].value_counts()
    multi = sizes[sizes > 1]
    tcrs["cluster_size"] = tcrs["cluster_id"].map(sizes)
    tcrs["in_multimember_cluster"] = tcrs["cluster_size"] > 1

    # A cluster whose members come from several samples is the interesting
    # case: convergent recombination onto a shared antigen, rather than one
    # patient's clonal lineage.
    spread = tcrs.groupby("cluster_id")["subject"].nunique()
    tcrs["cluster_n_samples"] = tcrs["cluster_id"].map(spread)

    logging.info(
        "%d clusters, %d with >1 member, %d spanning >1 sample (radius %d)",
        len(sizes), len(multi), int((spread > 1).sum()), radius,
    )

    Path(snakemake.output.clusters).parent.mkdir(parents=True, exist_ok=True)
    tcrs.to_csv(snakemake.output.clusters, sep="\t", index=False)

    # --- database annotation ---------------------------------------------
    hits = []
    for name in p.reference_dbs:
        if name not in DB_FILES:
            logging.warning("Unknown database %r, skipping", name)
            continue
        db = load_db(p.resources_dir, name)
        if db is None:
            continue
        merged = tcrs.merge(
            db, left_on="cdr3_b_aa", right_on="cdr3b_aa", how="inner"
        )
        logging.info("%s: %d exact CDR3b matches", name, len(merged))
        if len(merged):
            hits.append(
                merged[
                    ["clone_id", "subject", "cluster_id", "cdr3_b_aa",
                     "epitope", "source", "database"]
                ]
            )

    hit_df = (
        pd.concat(hits, ignore_index=True)
        if hits
        else pd.DataFrame(
            columns=["clone_id", "subject", "cluster_id", "cdr3_b_aa",
                     "epitope", "source", "database"]
        )
    )
    if hit_df.empty:
        logging.info(
            "No database hits. This is common and is not evidence against the "
            "clusters; reference databases cover a narrow slice of epitope space."
        )
    hit_df.to_csv(snakemake.output.hits, sep="\t", index=False)

    # --- figure -----------------------------------------------------------
    figdir = Path(snakemake.output.figdir)
    figdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    counts = sizes.value_counts().sort_index()
    ax.bar(counts.index.astype(int), counts.to_numpy(), color="#4C72B0")
    ax.set_xlabel("Cluster size", fontsize=28)
    ax.set_ylabel("Clusters", fontsize=28)
    ax.set_yscale("log")
    ax.tick_params(labelsize=28)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(figdir / f"cluster_sizes.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
