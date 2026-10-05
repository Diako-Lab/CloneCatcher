"""CloneCatcher Module 1b -- CellRanger VDJ parse (single_cell mode).

Reads `filtered_contig_annotations.csv` for one sample and collapses it to one
row per cell with TRA and TRB resolved. CellRanger has already assembled the
contigs and called clonotypes, so this stage filters and pivots rather than
reconstructing.

Barcodes stay BARE here. Promotion into ClusterCatcher's namespace happens in
Module 2, which is the only stage with the AnnData in hand to verify the join.

Output: tcr/vdj/{sample}/contigs.tsv
"""

import logging
import sys
from pathlib import Path

import pandas as pd

# Published CellRanger output paths only. Copies also exist deep inside
# SC_MULTI_CS pipeline internals, but those paths carry content hashes, are
# not a supported interface, and can be removed by cleanup.
VDJ_PATTERNS = (
    "multi_runs/{s}/outs/per_sample_outs/{s}/vdj_t/filtered_contig_annotations.csv",
    "multi_runs/{s}/outs/multi/vdj_t/filtered_contig_annotations.csv",
    "{s}/outs/per_sample_outs/{s}/vdj_t/filtered_contig_annotations.csv",
    "{s}/outs/filtered_contig_annotations.csv",
)

CHAIN_COLUMNS = {
    "TRA": ("cdr3a_nt", "cdr3a_aa", "v_a", "j_a", "umis_a"),
    "TRB": ("cdr3b_nt", "cdr3b_aa", "v_b", "j_b", "umis_b"),
}


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def find_contigs(cellranger_dir, sample):
    cellranger_dir = Path(cellranger_dir)
    for pattern in VDJ_PATTERNS:
        path = cellranger_dir / pattern.format(s=sample)
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No filtered_contig_annotations.csv for sample {sample!r} under "
        f"{cellranger_dir}. Checked the published outs/ paths only. If this "
        "sample's multi config declares no VDJ library it should have been "
        "filtered out upstream; run `CloneCatcher inspect-barcodes` to check."
    )


def pick_best_chain(group, chain):
    """Highest-UMI productive contig for one chain in one cell.

    A cell can carry more than one productive contig per chain, from allelic
    inclusion or ambient contamination. Taking the top-UMI contig is the
    standard resolution and matches what CellRanger uses for its own
    clonotype definition.
    """
    rows = group[group["chain"] == chain]
    if rows.empty:
        return None
    return rows.sort_values("umis", ascending=False).iloc[0]


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params

    sample = str(p.sample)
    path = find_contigs(p.cellranger_dir, sample)
    logging.info("Sample %s: reading %s", sample, path)

    df = pd.read_csv(path)
    logging.info("  %d contig rows", len(df))

    # --- filter ---------------------------------------------------------
    before = len(df)
    if "is_cell" in df.columns:
        df = df[df["is_cell"].astype(str).str.lower().isin(("true", "1"))]
    if "high_confidence" in df.columns:
        df = df[df["high_confidence"].astype(str).str.lower().isin(("true", "1"))]
    if p.require_productive and "productive" in df.columns:
        df = df[df["productive"].astype(str).str.lower().isin(("true", "1"))]
    df = df[df["chain"].isin(CHAIN_COLUMNS)]
    df = df[df["umis"] >= int(p.min_umis)]
    logging.info("  %d rows after filtering (dropped %d)", len(df), before - len(df))

    if df.empty:
        raise ValueError(
            f"Sample {sample}: no contigs survived filtering. Check min_umis "
            f"({p.min_umis}) and require_productive ({p.require_productive})."
        )

    # --- collapse to one row per cell ------------------------------------
    records = []
    for barcode, group in df.groupby("barcode", sort=True):
        rec = {"original_barcode": str(barcode), "sample_id": sample}
        n_chains = 0
        for chain, (nt, aa, v, j, umi) in CHAIN_COLUMNS.items():
            best = pick_best_chain(group, chain)
            if best is None:
                rec.update({nt: "", aa: "", v: "", j: "", umi: 0})
                continue
            n_chains += 1
            rec[nt] = str(best.get("cdr3_nt", "") or "")
            rec[aa] = str(best.get("cdr3", "") or "")
            rec[v] = str(best.get("v_gene", "") or "")
            rec[j] = str(best.get("j_gene", "") or "")
            rec[umi] = int(best["umis"])

        # D gene is TRB only.
        best_b = pick_best_chain(group, "TRB")
        rec["d_b"] = str(best_b.get("d_gene", "") or "") if best_b is not None else ""

        # CellRanger clonotype IDs are 'clonotype1', 'clonotype2', ... and are
        # unique only WITHIN a sample. Namespacing them the same way barcodes
        # are namespaced keeps them from colliding once samples are pooled.
        raw = group["raw_clonotype_id"].dropna()
        raw = raw[raw.astype(str).str.startswith("clonotype")]
        rec["clonotype_id"] = f"{sample}_{raw.iloc[0]}" if len(raw) else ""

        rec["n_chains"] = n_chains
        rec["is_paired"] = n_chains == 2
        records.append(rec)

    out = pd.DataFrame.from_records(records)

    if p.require_paired:
        before = len(out)
        out = out[out["is_paired"]]
        logging.info("  require_paired: kept %d of %d cells", len(out), before)

    unassigned = (out["clonotype_id"] == "").sum()
    if unassigned:
        logging.warning(
            "  %d cells have no CellRanger clonotype assignment", unassigned
        )

    logging.info(
        "Sample %s: %d cells, %d paired (%.1f%%), %d clonotypes",
        sample,
        len(out),
        int(out["is_paired"].sum()),
        100 * out["is_paired"].mean() if len(out) else 0.0,
        out.loc[out["clonotype_id"] != "", "clonotype_id"].nunique(),
    )

    Path(snakemake.output.contigs).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(snakemake.output.contigs, sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
