#!/usr/bin/env python
"""
Inspect the barcode namespace of a ClusterCatcher run and report the join rate
against CellRanger VDJ output.

Run this before configuring CloneCatcher in single_cell mode. It answers three
questions that every downstream join depends on:

  1. What is `sample_id` in this run -- SRA accession or multi-run name?
  2. Are the AnnData barcodes canonical, or from before the namespace fix?
  3. What fraction of VDJ cells will actually join to the GEX matrix?

Usage:
    python inspect_barcodes.py --run-dir /master/jlehle/WORKING/SC/fastq/OPSSC/GSE226620
    python inspect_barcodes.py --run-dir <dir> --sample HPV01
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

CANONICAL = re.compile(r"^[^_]+_[ACGT]+-\d+$")   # SRR13177101_AAACCT...-1
LEGACY = re.compile(r"^[ACGT]+-\d+-\S+$")        # AAACCT...-1-SRR13177101  (pre-fix)
BARE = re.compile(r"^[ACGT]+(-\d+)?$")           # AAACCT...-1

# Searched in order. The signature outputs are absent during a signature rerun,
# so fall back to the most recent AnnData that exists.
ADATA_CANDIDATES = [
    "signatures/adata_final.h5ad",
    "adata_final.h5ad",
    "dysregulation/adata_cancer_detected.h5ad",
    "viral_integration/adata_viral_integrated.h5ad",
    "viral_integration/adata_with_virus.h5ad",
    "annotation/adata_annotated.h5ad",
]


def classify(barcode):
    if CANONICAL.match(barcode):
        return "canonical"
    if LEGACY.match(barcode):
        return "legacy (pre-fix)"
    if BARE.match(barcode):
        return "bare"
    return "unrecognized"


def find_adata(run_dir):
    for rel in ADATA_CANDIDATES:
        path = run_dir / rel
        if path.exists():
            return path
    return None


def find_vdj(run_dir, sample):
    """Locate filtered_contig_annotations.csv for a sample.

    Only the published `outs/` paths are searched. CellRanger also leaves
    copies deep inside SC_MULTI_CS pipeline internals, and those are not a
    supported interface -- their paths carry content hashes and they can be
    removed by cleanup without warning.
    """
    patterns = [
        f"cellranger/multi_runs/{sample}/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv",
        f"cellranger/multi_runs/{sample}/outs/multi/vdj_t/filtered_contig_annotations.csv",
        f"cellranger/{sample}/outs/per_sample_outs/{sample}/vdj_t/filtered_contig_annotations.csv",
        f"cellranger/{sample}/outs/filtered_contig_annotations.csv",
    ]
    for rel in patterns:
        path = run_dir / rel
        if path.exists():
            return path
    return None


def read_multi_libraries(run_dir, sample):
    """Return the [libraries] rows from a sample's multi config, as
    (fastq_id, feature_type) pairs. Empty list if no config is found.

    A sample whose config declares only Gene Expression produced no VDJ
    output by design, which is a different situation from a missing file.
    """
    candidates = [
        run_dir / "cellranger" / "multi_runs" / f"{sample}_multi_config.csv",
        run_dir / "cellranger" / "multi_runs" / sample / "outs" / "config.csv",
    ]
    for path in candidates:
        if not path.exists():
            continue
        rows, in_block = [], False
        for line in path.read_text().splitlines():
            line = line.strip()
            if line.startswith("["):
                in_block = line.lower().startswith("[libraries]")
                continue
            if not in_block or not line or line.startswith("fastq_id"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 3:
                rows.append((parts[0], parts[2]))
        if rows:
            return rows
    return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--sample", default=None,
                    help="Sample to test the VDJ join against. Defaults to the first found.")
    args = ap.parse_args()

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        sys.exit(f"Not a directory: {run_dir}")

    import scanpy as sc

    # ---- 1. AnnData barcode namespace -------------------------------------
    adata_path = find_adata(run_dir)
    if adata_path is None:
        sys.exit(f"No AnnData found under {run_dir}. Tried:\n  " +
                 "\n  ".join(ADATA_CANDIDATES))

    print(f"AnnData:  {adata_path.relative_to(run_dir)}")
    adata = sc.read_h5ad(adata_path)
    obs_names = adata.obs_names.astype(str)

    print(f"Cells:    {len(obs_names):,}")
    print("\nFirst 5 barcodes:")
    for bc in obs_names[:5]:
        print(f"  {bc!r}   [{classify(bc)}]")

    forms = pd.Series([classify(b) for b in obs_names]).value_counts()
    print("\nBarcode forms:")
    for form, n in forms.items():
        print(f"  {form:<18} {n:>8,}")

    if "legacy (pre-fix)" in forms:
        print("\n  *** This run predates the barcode namespace fix. ***")
        print("  *** Run repair_barcodes.py before using it as CloneCatcher input. ***")

    # ---- 2. What sample_id actually is -------------------------------------
    print("\nobs columns relevant to identity:")
    for col in ("sample_id", "original_barcode", "batch", "donor", "condition"):
        if col in adata.obs.columns:
            vals = adata.obs[col].astype(str)
            uniq = vals.unique()
            shown = ", ".join(map(repr, uniq[:6]))
            more = f" ... (+{len(uniq) - 6})" if len(uniq) > 6 else ""
            print(f"  {col:<18} {len(uniq):>4} unique: {shown}{more}")
        else:
            print(f"  {col:<18} ABSENT")

    if "sample_id" in adata.obs.columns:
        sample_ids = set(adata.obs["sample_id"].astype(str))
        prefixes = {b.split("_", 1)[0] for b in obs_names if "_" in b}
        if prefixes and prefixes != sample_ids:
            print("\n  *** sample_id does not match the barcode prefix. ***")
            print(f"  sample_id only:      {sorted(sample_ids - prefixes)[:5]}")
            print(f"  barcode prefix only: {sorted(prefixes - sample_ids)[:5]}")
        elif prefixes:
            print("\n  sample_id matches the barcode prefix.")

        # Does sample_id look like an SRA accession or a multi-run name?
        sra_like = sum(bool(re.match(r"^[SED]RR\d+$", s)) for s in sample_ids)
        print(f"  {sra_like}/{len(sample_ids)} sample_ids look like SRA run accessions.")

    # Compare against the cellranger directory names on disk.
    cr_dir = run_dir / "cellranger"
    if cr_dir.is_dir():
        dirs = sorted(p.name for p in cr_dir.iterdir()
                      if p.is_dir() and p.name != "multi_runs")
        print(f"\ncellranger/ sample directories ({len(dirs)}): {dirs[:8]}")
    multi_dir = cr_dir / "multi_runs"
    if multi_dir.is_dir():
        dirs = sorted(p.name for p in multi_dir.iterdir() if p.is_dir())
        print(f"cellranger/multi_runs/ ({len(dirs)}): {dirs[:8]}")

    # ---- 3. VDJ join rate --------------------------------------------------
    sample = args.sample
    if sample is None and "sample_id" in adata.obs.columns:
        sample = str(adata.obs["sample_id"].astype(str).iloc[0])
    if sample is None:
        print("\nNo sample to test the VDJ join against; pass --sample.")
        return

    print(f"\nVDJ join test for sample {sample!r}")

    libraries = read_multi_libraries(run_dir, sample)
    if libraries:
        print("  multi config declares:")
        for fastq_id, feature_type in libraries:
            print(f"    {fastq_id:<16} {feature_type}")
        has_vdj = any("vdj" in ft.lower() for _, ft in libraries)
        if not has_vdj:
            print("\n  *** This sample has no VDJ library. ***")
            print("  cellranger multi ran in gene-expression-only mode, so there")
            print("  is no vdj_t output to join. The sample cannot be analyzed in")
            print("  single_cell mode. Re-run multi with a VDJ-T library, or")
            print("  exclude it via `samples` in the CloneCatcher config.")
            return

    vdj_path = find_vdj(run_dir, sample)
    if vdj_path is None:
        print("  filtered_contig_annotations.csv not found under outs/.")
        if not libraries:
            print("  No multi config found either; check the cellranger layout.")
        return
    print(f"  {vdj_path.relative_to(run_dir)}")

    contigs = pd.read_csv(vdj_path)
    vdj_bare = contigs["barcode"].astype(str).unique()
    print(f"  VDJ cells: {len(vdj_bare):,}")
    print(f"  VDJ barcode example: {vdj_bare[0]!r}   [{classify(vdj_bare[0])}]")

    promoted = {f"{sample}_{b}" for b in vdj_bare}
    gex = set(obs_names)
    overlap = len(promoted & gex)
    rate = overlap / max(len(promoted), 1)

    print(f"  Promoted example:    {sorted(promoted)[0]!r}")
    print(f"  Join rate: {overlap:,}/{len(promoted):,} ({rate:.1%})")

    if rate == 0:
        print("\n  *** Zero overlap. This is a namespace mismatch, not missing data. ***")
        print(f"  GEX example:      {sorted(gex)[0]!r}")
        print(f"  Promoted example: {sorted(promoted)[0]!r}")
    elif rate < 0.30:
        print("\n  Low join rate. Expected below 100% because VDJ recovers T cells")
        print("  that GEX QC drops, but under 30% suggests a key problem.")
    else:
        print("\n  Join looks sane. Set barcode.min_join_rate just below this.")


if __name__ == "__main__":
    main()
