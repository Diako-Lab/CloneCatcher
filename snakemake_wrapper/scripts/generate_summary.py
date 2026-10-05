"""CloneCatcher Module 7 -- Master summary.

Collects the headline number from every stage into one YAML, so a run can be
checked without opening six tables.

The summary reports what it actually found rather than what it expected. A
stage that did not run is recorded as absent, not as zero: those are different
claims, and conflating them is how the original ClusterCatcher signature bug
reported EXCELLENT on a reconstruction of nothing.

Output: master_summary.yaml
"""

import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import yaml


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def read_tsv(path):
    path = Path(path)
    if not path.exists():
        return None
    try:
        df = pd.read_csv(path, sep="\t")
        return df if len(df) else None
    except Exception as exc:  # noqa: BLE001
        logging.warning("Could not read %s: %s", path, exc)
        return None


def summarize_clonotypes(out_dir):
    df = read_tsv(out_dir / "tcr" / "clonotypes_by_sample.tsv")
    if df is None:
        return {"status": "absent"}

    df["is_expanded"] = df["is_expanded"].astype(str).str.lower().isin(("true", "1"))
    df["is_paired"] = df["is_paired"].astype(str).str.lower().isin(("true", "1"))

    per_sample = {}
    for sample, g in df.groupby("sample_id"):
        per_sample[str(sample)] = {
            "clonotypes": int(len(g)),
            "observations": int(g["n_observations"].sum()),
            "expanded": int(g["is_expanded"].sum()),
            "top_clone_frequency": round(float(g["frequency"].max()), 4),
        }

    return {
        "status": "ok",
        "samples": int(df["sample_id"].nunique()),
        "clonotypes_total": int(len(df)),
        "clonotypes_expanded": int(df["is_expanded"].sum()),
        "paired_fraction": round(float(df["is_paired"].mean()), 4),
        "per_sample": per_sample,
    }


def summarize_cells(out_dir):
    df = read_tsv(out_dir / "tcr" / "clonotypes_by_cell.tsv")
    if df is None:
        return {"status": "absent"}
    df["is_paired"] = df["is_paired"].astype(str).str.lower().isin(("true", "1"))
    return {
        "status": "ok",
        "cells": int(len(df)),
        "clonotypes": int(df["clonotype_id"].nunique()),
        "paired_fraction": round(float(df["is_paired"].mean()), 4),
    }


def summarize_join(out_dir):
    df = read_tsv(out_dir / "tcr" / "barcode_join_report.tsv")
    if df is None:
        return {"status": "absent"}
    overall = df["matched"].sum() / max(df["vdj_cells"].sum(), 1)
    return {
        "status": "ok",
        "overall_join_rate": round(float(overall), 4),
        "worst_sample": str(df.loc[df["join_rate"].idxmin(), "sample_id"]),
        "worst_join_rate": round(float(df["join_rate"].min()), 4),
        "unmapped_total": int(df["unmapped"].sum()),
    }


def summarize_repertoire(out_dir):
    df = read_tsv(out_dir / "repertoire" / "repertoire_stats.tsv")
    if df is None:
        return {"status": "absent"}
    out = {"status": "ok", "samples": int(len(df))}
    if "clonality" in df.columns:
        out["clonality_median"] = round(float(df["clonality"].median()), 4)
        out["clonality_range"] = [
            round(float(df["clonality"].min()), 4),
            round(float(df["clonality"].max()), 4),
        ]
    return out


def summarize_expansion(out_dir):
    df = read_tsv(out_dir / "expansion" / "expansion_results.tsv")
    if df is None:
        return {"status": "absent"}
    test = str(df["test"].iloc[0]) if "test" in df.columns else "unknown"
    out = {"status": "ok", "test": test, "rows": int(len(df))}
    if test == "none":
        out["mode"] = "descriptive (no contrast declared)"
    else:
        if "significant" in df.columns:
            out["significant"] = int(df["significant"].fillna(False).sum())
        if "recurrent" in df.columns:
            out["recurrent_across_blocks"] = int(df["recurrent"].fillna(False).sum())
        if "qvalue" in df.columns and df["qvalue"].notna().any():
            out["min_qvalue"] = round(float(df["qvalue"].min()), 6)
    return out


def summarize_specificity(out_dir):
    clusters = read_tsv(out_dir / "specificity" / "specificity_clusters.tsv")
    if clusters is None:
        return {"status": "absent"}
    sizes = clusters["cluster_id"].value_counts()
    out = {
        "status": "ok",
        "receptors": int(len(clusters)),
        "clusters": int(len(sizes)),
        "multimember_clusters": int((sizes > 1).sum()),
        "largest_cluster": int(sizes.max()),
    }
    if "chains" in clusters.columns:
        out["chains"] = str(clusters["chains"].iloc[0])
    if "cluster_n_samples" in clusters.columns:
        spread = clusters.groupby("cluster_id")["cluster_n_samples"].first()
        out["clusters_spanning_samples"] = int((spread > 1).sum())

    hits = read_tsv(out_dir / "specificity" / "database_hits.tsv")
    out["database_hits"] = 0 if hits is None else int(len(hits))
    if hits is not None and "epitope" in hits.columns:
        out["distinct_epitopes"] = int(hits["epitope"].nunique())
    return out


def summarize_integration(out_dir):
    path = out_dir / "integration" / "adata_tcr.h5ad"
    if not path.exists():
        return {"status": "absent"}
    return {
        "status": "ok",
        "path": str(path),
        "size_mb": round(path.stat().st_size / 1e6, 1),
    }


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params
    out_dir = Path(snakemake.output.summary).parent

    summary = {
        "clonecatcher": {
            "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "input_mode": str(p.mode),
            "samples_requested": sorted(map(str, p.samples)),
            "n_samples_requested": len(p.samples),
        },
        "clonotypes": summarize_clonotypes(out_dir),
        "cells": summarize_cells(out_dir),
        "barcode_join": summarize_join(out_dir),
        "repertoire": summarize_repertoire(out_dir),
        "expansion": summarize_expansion(out_dir),
        "specificity": summarize_specificity(out_dir),
        "integration": summarize_integration(out_dir),
    }

    absent = [k for k, v in summary.items() if isinstance(v, dict) and v.get("status") == "absent"]
    if absent:
        logging.info("Stages with no output: %s", ", ".join(absent))

    for name in ("clonotypes", "repertoire", "expansion", "specificity"):
        logging.info("%-13s %s", name, summary[name])

    Path(snakemake.output.summary).parent.mkdir(parents=True, exist_ok=True)
    with open(snakemake.output.summary, "w") as fh:
        yaml.safe_dump(summary, fh, sort_keys=False, default_flow_style=False)
    logging.info("Wrote %s", snakemake.output.summary)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
