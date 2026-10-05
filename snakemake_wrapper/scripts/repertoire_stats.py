"""CloneCatcher Module 3 -- Repertoire statistics (Immunarch).

Clonality, diversity, and V/J usage per sample.

Python wrapping an R call, rather than a .R script, so scripts/ stays uniform
with ClusterCatcher and so the input reshaping and output validation live in
the same language as the rest of the pipeline. The same pattern ClusterCatcher
uses for inferCNV.

Outputs: repertoire/repertoire_stats.tsv, repertoire/figures/
"""

import logging
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

# Immunarch's expected input column names.
IMMUNARCH_COLUMNS = {
    "n_observations": "Clones",
    "frequency": "Proportion",
    "cdr3_nt": "CDR3.nt",
    "cdr3_aa": "CDR3.aa",
    "v_name": "V.name",
    "d_name": "D.name",
    "j_name": "J.name",
}

R_SCRIPT = r"""
suppressPackageStartupMessages({
  library(immunarch); library(dplyr); library(readr); library(ggplot2)
})

args     <- commandArgs(trailingOnly = TRUE)
in_dir   <- args[1]
out_tsv  <- args[2]
fig_dir  <- args[3]
methods  <- strsplit(args[4], ",")[[1]]
top_n    <- as.integer(args[5])

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

files   <- list.files(in_dir, pattern = "\\.tsv$", full.names = TRUE)
samples <- tools::file_path_sans_ext(basename(files))
data    <- lapply(files, function(f) readr::read_tsv(f, show_col_types = FALSE))
names(data) <- samples
immdata <- list(data = data, meta = data.frame(Sample = samples))

rows <- list()

# Clonality: 1 - normalised Shannon. 0 is a flat repertoire, 1 is monoclonal.
inv <- repDiversity(immdata$data, "inv.simp")
div <- repDiversity(immdata$data, "div")

for (s in samples) {
  d <- immdata$data[[s]]
  p <- d$Proportion / sum(d$Proportion)
  shannon <- -sum(p * log(p), na.rm = TRUE)
  rows[[s]] <- data.frame(
    sample_id     = s,
    n_clonotypes  = nrow(d),
    n_observations= sum(d$Clones),
    shannon       = shannon,
    clonality     = if (nrow(d) > 1) 1 - shannon / log(nrow(d)) else NA_real_,
    top_clone_freq= max(p, na.rm = TRUE),
    top10_freq    = sum(sort(p, decreasing = TRUE)[1:min(10, length(p))])
  )
}
stats <- dplyr::bind_rows(rows)

if ("chao1" %in% methods) {
  ch <- repDiversity(immdata$data, "chao1")
  stats$chao1 <- as.numeric(ch[match(stats$sample_id, rownames(ch)), "Estimator"])
}
if ("gini" %in% methods) {
  stats$gini <- sapply(samples, function(s) {
    x <- sort(immdata$data[[s]]$Clones); n <- length(x)
    if (n < 2) return(NA_real_)
    sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
  })
}

readr::write_tsv(stats, out_tsv)

save_both <- function(plot, name) {
  for (ext in c("pdf", "png")) {
    ggplot2::ggsave(file.path(fig_dir, paste0(name, ".", ext)),
                    plot, width = 10, height = 7, dpi = 300)
  }
}

try({ save_both(vis(repExplore(immdata$data, "volume")), "clonotype_counts") }, silent = TRUE)
try({ save_both(vis(repClonality(immdata$data, "top", .head = c(10, 100, 1000))),
                "clonal_proportion") }, silent = TRUE)
try({ save_both(vis(geneUsage(immdata$data, "hs.trbv", .norm = TRUE)),
                "trbv_usage") }, silent = TRUE)

cat("immunarch: wrote", nrow(stats), "sample rows\n")
"""


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def to_immunarch(df, chain="b"):
    """Reshape the shared schema into Immunarch's expected columns.

    Keyed on the beta chain by default. Beta is the conventional clonotype
    definition, it is the only chain bulk mode resolves, and using it in both
    modes keeps the statistics comparable between arms.
    """
    out = pd.DataFrame(
        {
            "Clones": df["n_observations"].astype(int),
            "Proportion": df["frequency"].astype(float),
            "CDR3.nt": df[f"cdr3{chain}_nt"].fillna(""),
            "CDR3.aa": df[f"cdr3{chain}_aa"].fillna(""),
            "V.name": df[f"v_{chain}"].fillna(""),
            "D.name": df.get("d_b", pd.Series([""] * len(df))).fillna(""),
            "J.name": df[f"j_{chain}"].fillna(""),
        }
    )
    return out[out["CDR3.aa"].astype(str).str.len() > 0]


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params

    df = pd.read_csv(snakemake.input.by_sample, sep="\t")
    logging.info("%d clonotypes across %d samples", len(df), df["sample_id"].nunique())

    fig_dir = Path(snakemake.output.figdir)
    fig_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        in_dir = tmp / "immunarch_in"
        in_dir.mkdir()

        written = 0
        for sample, g in df.groupby("sample_id"):
            conv = to_immunarch(g)
            if conv.empty:
                logging.warning("%s has no beta-chain clonotypes, skipping", sample)
                continue
            conv.to_csv(in_dir / f"{sample}.tsv", sep="\t", index=False)
            written += 1

        if not written:
            raise ValueError(
                "No sample had a usable beta chain. Immunarch statistics are "
                "keyed on TRB; check the upstream frontend."
            )

        script = tmp / "repertoire_stats.R"
        script.write_text(R_SCRIPT)

        cmd = [
            "Rscript", str(script), str(in_dir), str(snakemake.output.stats),
            str(fig_dir), ",".join(map(str, p.diversity_methods)),
            str(p.top_n_clones),
        ]
        logging.info("Running: %s", " ".join(cmd))
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.stdout:
            logging.info("R stdout:\n%s", proc.stdout.strip())
        if proc.returncode != 0:
            logging.error("R stderr:\n%s", proc.stderr.strip())
            raise RuntimeError(
                f"Immunarch failed (exit {proc.returncode}). If the error names a "
                "missing package, the environment did not solve; build it with "
                "`snakemake --use-conda --conda-create-envs-only`."
            )

    stats = pd.read_csv(snakemake.output.stats, sep="\t")
    logging.info("Wrote %d sample rows", len(stats))
    for _, r in stats.iterrows():
        logging.info(
            "  %-8s %6d clonotypes  clonality %.3f  top clone %.3f",
            r["sample_id"], r["n_clonotypes"],
            r.get("clonality", float("nan")), r.get("top_clone_freq", float("nan")),
        )


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
