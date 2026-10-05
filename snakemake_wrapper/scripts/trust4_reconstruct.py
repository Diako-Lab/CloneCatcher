"""CloneCatcher Module 1a -- TRUST4 reconstruction (bulk mode).

Assembles TCR contigs from one sample's 5' RACE libraries and emits the shared
clonotype schema.

A LIMIT WORTH KNOWING. Bulk TCR sequencing cannot pair alpha with beta. The
chains come from different mRNA molecules in a pooled lysate, and nothing in
the library preserves which cell they shared. TRUST4 reports each chain
independently, so in bulk mode `is_paired` is always False and a clonotype is
defined on TRB alone -- which is the conventional bulk definition, beta being
the more diverse chain and the one carrying most of the specificity signal.

The consequence downstream is that specificity clustering in bulk mode runs on
single-chain receptors, which is weaker evidence than paired. That is a
property of the assay, not of this pipeline, and it is why the single-cell
arm is the stronger place to establish a neoantigen-TCR link.

Output: tcr/trust4/{sample}/clonotypes.tsv
"""

import logging
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

# TRUST4's report.tsv columns, in order.
REPORT_COLUMNS = [
    "count",
    "frequency",
    "CDR3nt",
    "CDR3aa",
    "V",
    "D",
    "J",
    "C",
    "cid",
    "cid_full_length",
]

CHAIN_PREFIX = {"TRA": "a", "TRB": "b"}


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def find_trust4():
    """Locate run-trust4, preferring the vendored build."""
    vendored = Path("external/TRUST4/run-trust4")
    if vendored.exists():
        return str(vendored.resolve())
    found = shutil.which("run-trust4")
    if found:
        return found
    raise FileNotFoundError(
        "run-trust4 not found. Build the vendored copy with "
        "`cd external/TRUST4 && make`, or install the bioconda package."
    )


def sample_fastqs(sample_info, sample):
    with open(sample_info, "rb") as fh:
        samples = pickle.load(fh)
    if sample not in samples:
        raise KeyError(f"Sample {sample!r} not in {sample_info}")
    meta = samples[sample]
    r1 = str(meta["fastq_r1"])
    r2 = str(meta.get("fastq_r2", "") or "")
    return r1, r2


def chain_of(v_gene, j_gene, c_gene):
    """Assign a row to TRA or TRB from whichever gene call is present."""
    for gene in (v_gene, j_gene, c_gene):
        g = str(gene or "")
        if g.startswith("TRA"):
            return "TRA"
        if g.startswith("TRB"):
            return "TRB"
    return None


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params
    sample = str(snakemake.wildcards.sample)
    outdir = Path(p.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    r1, r2 = sample_fastqs(p.sample_info, sample)
    trust4 = find_trust4()

    cmd = [trust4, "-f", str(p.imgt), "-o", sample, "--od", str(outdir),
           "-t", str(snakemake.threads)]
    cmd += ["-1", r1, "-2", r2] if r2 else ["-u", r1]

    logging.info("Sample %s: %s", sample, " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.stdout:
        logging.info("TRUST4 stdout:\n%s", proc.stdout.strip())
    if proc.returncode != 0:
        logging.error("TRUST4 stderr:\n%s", proc.stderr.strip())
        raise RuntimeError(f"TRUST4 failed for {sample} (exit {proc.returncode})")

    report = outdir / f"{sample}_report.tsv"
    if not report.exists():
        raise FileNotFoundError(
            f"TRUST4 exited 0 but wrote no report for {sample}: {report}"
        )

    df = pd.read_csv(report, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    missing = [c for c in ("count", "CDR3aa", "V", "J") if c not in df.columns]
    if missing:
        raise ValueError(
            f"{report} is missing expected column(s) {missing}. TRUST4's "
            f"report format may have changed; columns present: {list(df.columns)}"
        )
    logging.info("  %d raw rows", len(df))

    df["chain"] = [
        chain_of(v, j, c)
        for v, j, c in zip(df["V"], df["J"], df.get("C", [""] * len(df)))
    ]
    before = len(df)
    df = df[df["chain"].isin(CHAIN_PREFIX)]
    df = df[df["count"] >= int(p.min_abundance)]
    df = df[df["CDR3aa"].astype(str).str.len() > 0]
    # TRUST4 marks unresolvable CDR3s with a bare underscore or a stop codon.
    df = df[~df["CDR3aa"].astype(str).str.contains(r"[\*_]", regex=True)]
    logging.info("  %d rows after filtering (dropped %d)", len(df), before - len(df))

    keep = [c for c in str(p.chains).split(",") if c]
    df = df[df["chain"].isin(keep)]

    if df.empty:
        raise ValueError(
            f"Sample {sample}: no clonotypes survived filtering. Check "
            f"min_abundance ({p.min_abundance}) and that the libraries are 5' "
            "chemistry -- TRUST4 cannot reconstruct receptors from 3' data."
        )

    # Emit the shared schema. Alpha and beta land in their own column sets and
    # are never joined, because bulk cannot pair them.
    records = []
    for _, row in df.iterrows():
        pre = CHAIN_PREFIX[row["chain"]]
        rec = {
            "sample_id": sample,
            "cdr3a_nt": "", "cdr3a_aa": "", "v_a": "", "j_a": "",
            "cdr3b_nt": "", "cdr3b_aa": "", "v_b": "", "d_b": "", "j_b": "",
            "n_observations": int(row["count"]),
            "is_paired": False,
        }
        rec[f"cdr3{pre}_nt"] = str(row.get("CDR3nt", "") or "")
        rec[f"cdr3{pre}_aa"] = str(row["CDR3aa"])
        rec[f"v_{pre}"] = str(row["V"] or "")
        rec[f"j_{pre}"] = str(row["J"] or "")
        if pre == "b":
            rec["d_b"] = str(row.get("D", "") or "")
        records.append(rec)

    out = pd.DataFrame.from_records(records)
    counts = out.apply(
        lambda r: "TRB" if r["cdr3b_aa"] else "TRA", axis=1
    ).value_counts()
    logging.info(
        "Sample %s: %d clonotypes (%s)",
        sample,
        len(out),
        ", ".join(f"{k} {v}" for k, v in counts.items()),
    )

    Path(snakemake.output.clonotypes).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(snakemake.output.clonotypes, sep="\t", index=False)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
