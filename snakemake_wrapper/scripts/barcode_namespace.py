"""CloneCatcher Module 2 -- Barcode namespace (single_cell mode).

Promotes bare CellRanger VDJ barcodes into ClusterCatcher's canonical
namespace and asserts the join against the GEX matrix.

This is the stage where the original ClusterCatcher namespace bug would have
been caught on day one. It does not write an empty table on a failed join; it
raises, with an example barcode from each side.

A VDJ cell absent from the GEX matrix is legitimate, since VDJ libraries
recover T cells that GEX QC drops. The expected rate is below 100%, which is
why min_join_rate is configured per run rather than hardcoded.

Writes an INTERMEDIATE, not a final table. Clonotype identity is assigned in
collapse_to_sample.py, which emits clonotypes_by_cell.tsv with IDs that
actually join to clonotypes_by_sample.tsv.

Outputs: tcr/cells_namespaced.tsv, tcr/barcode_join_report.tsv
"""

import logging
import sys
from pathlib import Path

import pandas as pd

import barcode_utils as bu

CELL_COLUMNS = [
    "cell_barcode",
    "sample_id",
    "original_barcode",
    "cellranger_clonotype_id",
    "n_chains",
    "is_paired",
    "cdr3a_nt",
    "cdr3a_aa",
    "v_a",
    "j_a",
    "umis_a",
    "cdr3b_nt",
    "cdr3b_aa",
    "v_b",
    "d_b",
    "j_b",
    "umis_b",
]


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params
    min_rate = float(p.min_join_rate)
    on_unmapped = str(p.on_unmapped).lower()

    if on_unmapped not in ("drop", "raise"):
        raise ValueError(
            f"barcode.on_unmapped must be 'drop' or 'raise', got {on_unmapped!r}"
        )

    # --- the GEX side ----------------------------------------------------
    obs, adata_path = bu.load_clustercatcher_obs(
        p.clustercatcher_run_dir, columns=["sample_id", "original_barcode"]
    )
    logging.info("ClusterCatcher AnnData: %s", adata_path)
    logging.info("  %d cells, %d samples", len(obs), obs["sample_id"].nunique())

    # Fails loudly on a pre-fix run rather than producing an empty join.
    bu.assert_canonical(obs.index, source=str(adata_path))

    gex_by_sample = {
        s: set(g.index.astype(str))
        for s, g in obs.groupby(obs["sample_id"].astype(str))
    }

    # --- the VDJ side ----------------------------------------------------
    frames, report = [], []

    for path in map(Path, snakemake.input.contigs):
        df = pd.read_csv(path, sep="\t", dtype={"original_barcode": str})
        if df.empty:
            logging.warning("%s is empty, skipping", path)
            continue

        sample = str(df["sample_id"].iloc[0])
        gex = gex_by_sample.get(sample)
        if gex is None:
            raise bu.NamespaceError(
                f"Sample {sample!r} has VDJ output but does not appear in the "
                f"ClusterCatcher AnnData, which holds: "
                f"{sorted(gex_by_sample)[:10]}. The two runs cover different "
                "samples, or sample_id differs between them."
            )

        df["cell_barcode"] = bu.promote(df["original_barcode"], sample)

        n_vdj = len(df)
        matched = df["cell_barcode"].isin(gex)
        overlap, rate = int(matched.sum()), float(matched.mean())

        logging.info(
            "%-8s VDJ %6d  GEX %6d  matched %6d  (%.1f%%)",
            sample,
            n_vdj,
            len(gex),
            overlap,
            100 * rate,
        )

        report.append(
            {
                "sample_id": sample,
                "vdj_cells": n_vdj,
                "gex_cells": len(gex),
                "matched": overlap,
                "join_rate": round(rate, 4),
                "unmapped": n_vdj - overlap,
            }
        )

        # The tripwire. Raises on zero overlap or below the configured floor.
        bu.assert_join_rate(
            df["cell_barcode"], gex, min_rate=min_rate, label=f"VDJ->GEX [{sample}]"
        )

        if not matched.all():
            n_unmapped = int((~matched).sum())
            if on_unmapped == "raise":
                example = df.loc[~matched, "cell_barcode"].iloc[0]
                raise bu.NamespaceError(
                    f"{sample}: {n_unmapped} VDJ cells are absent from the GEX "
                    f"matrix (e.g. {example!r}) and on_unmapped is 'raise'."
                )
            logging.info(
                "  dropping %d unmapped VDJ cells (QC-filtered on the GEX side)",
                n_unmapped,
            )
            df = df[matched]

        frames.append(df)

    if not frames:
        raise ValueError("No VDJ contigs were read. Check the upstream stage.")

    cells = pd.concat(frames, ignore_index=True)

    for col in CELL_COLUMNS:
        if col not in cells.columns:
            cells[col] = ""
    cells = cells[CELL_COLUMNS].sort_values(["sample_id", "cell_barcode"])

    # Barcodes are unique per sample, so a duplicate here means a sample was
    # read twice or two samples share an id.
    dupes = cells["cell_barcode"].duplicated().sum()
    if dupes:
        raise bu.NamespaceError(
            f"{dupes} duplicate cell barcodes after pooling. A sample was "
            "processed twice, or two samples share a sample_id."
        )

    rep = pd.DataFrame(report).sort_values("sample_id")
    overall = rep["matched"].sum() / max(rep["vdj_cells"].sum(), 1)
    logging.info(
        "Pooled: %d cells across %d samples, overall join rate %.1f%%",
        len(cells),
        len(rep),
        100 * overall,
    )

    for out in (snakemake.output.cells, snakemake.output.report):
        Path(out).parent.mkdir(parents=True, exist_ok=True)

    cells.to_csv(snakemake.output.cells, sep="\t", index=False)
    rep.to_csv(snakemake.output.report, sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
