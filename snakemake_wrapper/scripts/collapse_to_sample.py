"""CloneCatcher -- schema convergence.

Both frontends meet here. Everything downstream reads clonotypes_by_sample.tsv
and nothing else, which is what keeps bulk and single-cell analysis from
drifting apart.

This stage owns clonotype identity. It groups on chain sequence and gene usage
and assigns `{sample_id}_ct####`, ranked by abundance within sample. That
definition has to live in one place for both modes, because bulk has no
CellRanger clonotype calls and a definition that differed between arms would
make them incomparable. CellRanger's own call is preserved alongside as
`cellranger_clonotype_id` so a cell can still be traced back to it.

In single_cell mode this stage emits BOTH final tables, so the clonotype_id in
clonotypes_by_cell.tsv is the same one in clonotypes_by_sample.tsv and the
documented join between them actually works.

Outputs: tcr/clonotypes_by_sample.tsv, and in single_cell mode
         tcr/clonotypes_by_cell.tsv
"""

import logging
import pickle
import sys
from pathlib import Path

import pandas as pd

SAMPLE_COLUMNS = [
    "sample_id",
    "clonotype_id",
    "cdr3a_nt",
    "cdr3a_aa",
    "v_a",
    "j_a",
    "cdr3b_nt",
    "cdr3b_aa",
    "v_b",
    "d_b",
    "j_b",
    "n_observations",
    "frequency",
    "is_paired",
    "is_expanded",
    "donor_id",
    "condition",
    "source_mode",
]

CELL_COLUMNS = [
    "cell_barcode",
    "sample_id",
    "clonotype_id",
    "cellranger_clonotype_id",
    "original_barcode",
    "n_chains",
    "is_paired",
]

# Chain identity. Two cells carry the same clonotype when all of these agree.
CHAIN_KEYS = ["cdr3a_aa", "v_a", "j_a", "cdr3b_aa", "v_b", "j_b"]

# Carried through from the first member of each clonotype group.
CARRY = ["cdr3a_nt", "cdr3b_nt", "d_b"]


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def load_sample_annotations(sample_info):
    """donor_id and condition per sample, from the bulk sample pickle.

    Absent in single_cell mode, where a tumor series has no experimental arms
    and the columns stay empty. The expansion stage reads `contrast: null`
    off the config in that case and reports descriptively.
    """
    if not sample_info:
        return {}
    with open(sample_info, "rb") as fh:
        samples = pickle.load(fh)
    out = {}
    for sid, meta in samples.items():
        if isinstance(meta, dict):
            out[str(sid)] = {
                "donor_id": str(meta.get("donor", "") or ""),
                "condition": str(meta.get("condition", "") or ""),
            }
    return out


def read_cells(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if df.empty:
        raise ValueError(f"{path} is empty.")
    df["is_paired"] = df["is_paired"].astype(str).str.lower().isin(("true", "1"))
    for col in CHAIN_KEYS + CARRY:
        if col not in df.columns:
            df[col] = ""
    return df


def collapse_cells(cells):
    """One row per clonotype per sample; n_observations is a cell count."""
    agg = {c: "first" for c in CARRY}
    agg["is_paired"] = "first"
    agg["cell_barcode"] = "count"
    return (
        cells.groupby(["sample_id"] + CHAIN_KEYS, dropna=False)
        .agg(agg)
        .rename(columns={"cell_barcode": "n_observations"})
        .reset_index()
    )


def collapse_bulk(paths):
    """One row per clonotype per sample; n_observations is read support.

    Expects each TRUST4 per-sample table to carry sample_id, the CHAIN_KEYS,
    and n_observations. That contract is the frontend's responsibility.
    """
    frames = []
    for path in map(Path, paths):
        df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
        if df.empty:
            logging.warning("%s is empty, skipping", path)
            continue
        missing = [c for c in ["sample_id", "n_observations"] if c not in df.columns]
        if missing:
            raise ValueError(f"{path} is missing required column(s): {missing}")
        frames.append(df)

    if not frames:
        raise ValueError("No bulk clonotype tables were read.")

    df = pd.concat(frames, ignore_index=True)
    df["n_observations"] = pd.to_numeric(df["n_observations"], errors="coerce").fillna(0)

    for col in CHAIN_KEYS + CARRY:
        if col not in df.columns:
            df[col] = ""

    if "is_paired" in df.columns:
        df["is_paired"] = df["is_paired"].astype(str).str.lower().isin(("true", "1"))
    else:
        df["is_paired"] = (df["cdr3a_aa"] != "") & (df["cdr3b_aa"] != "")

    agg = {c: "first" for c in CARRY}
    agg["is_paired"] = "first"
    agg["n_observations"] = "sum"
    return df.groupby(["sample_id"] + CHAIN_KEYS, dropna=False).agg(agg).reset_index()


def assign_ids_and_frequency(out, threshold):
    """Within-sample frequency, expansion call, and stable clonotype IDs.

    NOTE ON THE DENOMINATOR. `frequency` is a clonotype's share of the
    receptor-bearing observations in its own sample: cells with a recovered
    TCR in single_cell mode, reads assigned to a clonotype in bulk. It is NOT
    a share of all cells in the sample, nor of all T cells.

    That is the right denominator for clonality and for the beta-binomial,
    which compares a clone against the rest of its own repertoire. It is the
    only quantity comparable across libraries of different depth, and it is
    what keeps the two modes on the same footing. Expressing expansion
    relative to a per-tumor T cell fraction instead would make frequencies
    depend on cell type composition, which varies far more between tumors than
    the repertoire does, and across patients with different HLA types the
    comparison would be underpowered besides.
    """
    out["n_observations"] = out["n_observations"].astype(int)

    totals = out.groupby("sample_id")["n_observations"].transform("sum")
    out["frequency"] = (out["n_observations"] / totals.where(totals > 0, 1)).round(6)
    out["is_expanded"] = out["frequency"] >= threshold

    out = out.sort_values(["sample_id", "n_observations"], ascending=[True, False])
    out["clonotype_id"] = (
        out["sample_id"].astype(str)
        + "_ct"
        + (out.groupby("sample_id").cumcount() + 1).astype(str).str.zfill(4)
    )
    return out


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params
    mode = str(p.mode)
    threshold = float(p.expansion_threshold)

    cells = None
    if mode == "single_cell":
        cells = read_cells(snakemake.input[0])
        out = collapse_cells(cells)
    elif mode == "bulk":
        out = collapse_bulk(snakemake.input)
    else:
        raise ValueError(f"Unknown mode {mode!r}")

    out = assign_ids_and_frequency(out, threshold)

    annotations = load_sample_annotations(getattr(p, "sample_info", None))
    out["donor_id"] = out["sample_id"].map(
        lambda s: annotations.get(str(s), {}).get("donor_id", "")
    )
    out["condition"] = out["sample_id"].map(
        lambda s: annotations.get(str(s), {}).get("condition", "")
    )
    out["source_mode"] = mode

    for col in SAMPLE_COLUMNS:
        if col not in out.columns:
            out[col] = ""
    by_sample = out[SAMPLE_COLUMNS]

    logging.info(
        "%s: %d clonotypes across %d samples; %d expanded at freq >= %.3f",
        mode,
        len(by_sample),
        by_sample["sample_id"].nunique(),
        int(by_sample["is_expanded"].sum()),
        threshold,
    )
    for sample, g in by_sample.groupby("sample_id"):
        logging.info(
            "  %-8s %5d clonotypes  %5d expanded  top freq %.3f",
            sample,
            len(g),
            int(g["is_expanded"].sum()),
            g["frequency"].max(),
        )

    Path(snakemake.output.by_sample).parent.mkdir(parents=True, exist_ok=True)
    by_sample.to_csv(snakemake.output.by_sample, sep="\t", index=False)

    # --- cell table, with the SAME clonotype_id ---------------------------
    if cells is None:
        return

    merged = cells.merge(
        out[["sample_id"] + CHAIN_KEYS + ["clonotype_id"]],
        on=["sample_id"] + CHAIN_KEYS,
        how="left",
        validate="many_to_one",
    )

    unassigned = merged["clonotype_id"].isna().sum()
    if unassigned:
        # Every cell was grouped from this same frame, so a miss here means
        # the merge keys disagree with the grouping keys.
        raise ValueError(
            f"{unassigned} cells did not receive a clonotype_id. The merge "
            "keys and the grouping keys have diverged."
        )

    for col in CELL_COLUMNS:
        if col not in merged.columns:
            merged[col] = ""
    by_cell = merged[CELL_COLUMNS].sort_values(["sample_id", "cell_barcode"])

    logging.info(
        "cell table: %d cells mapped onto %d clonotypes",
        len(by_cell),
        by_cell["clonotype_id"].nunique(),
    )

    Path(snakemake.output.by_cell).parent.mkdir(parents=True, exist_ok=True)
    by_cell.to_csv(snakemake.output.by_cell, sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
