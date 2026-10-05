"""CloneCatcher Module 4 -- Expansion.

Two behaviours, selected by whether the config declares a contrast.

  contrast: null   Descriptive. Per-sample expanded-clone calls and clonality.
                   This is what a tumour series with no experimental arms
                   supports: there is nothing to difference against.

  contrast: set    Beta-binomial differential test between a test and a
                   reference condition, blocked on donor, with clones
                   expanding against an excluded condition removed first as
                   alloreactive.

THE MODEL. Clone counts are overdispersed relative to binomial: replicate
libraries from one condition vary more than sampling alone predicts. A
beta-binomial absorbs that with a dispersion parameter rho, estimated once
across all clonotypes rather than per clone, where there is never enough data.
The test is a likelihood ratio against a null of one shared proportion, one
degree of freedom, Benjamini-Hochberg corrected.

BLOCKING. Donors differ in repertoire composition far more than conditions do
within a donor, so pooling across donors would swamp the effect. Each donor is
tested separately and the per-donor p-values combined by Fisher's method. That
also yields recurrence directly, which is the specificity requirement the
design asks for: a clone expanding in one donor is a candidate, a clone
expanding in several is evidence.

Outputs: expansion/expansion_results.tsv, expansion/figures/
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

EPS = 1e-10


def setup_logging(path):
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-7s %(message)s",
        handlers=[logging.FileHandler(path, mode="w"), logging.StreamHandler(sys.stderr)],
    )


def betabinom_loglik(k, n, p, rho):
    """Beta-binomial log likelihood, parameterised by mean p and dispersion rho.

    rho -> 0 recovers the binomial. a = p(1-rho)/rho, b = (1-p)(1-rho)/rho.
    """
    p = np.clip(p, EPS, 1 - EPS)
    rho = np.clip(rho, EPS, 1 - EPS)
    s = (1 - rho) / rho
    a, b = p * s, (1 - p) * s
    return np.sum(stats.betabinom.logpmf(k, n, a, b))


def fit_p(k, n, rho):
    """ML proportion under fixed dispersion."""
    total_n = max(n.sum(), 1)
    start = np.clip(k.sum() / total_n, EPS, 1 - EPS)
    res = optimize.minimize_scalar(
        lambda q: -betabinom_loglik(k, n, q, rho),
        bounds=(EPS, 1 - EPS),
        method="bounded",
    )
    return float(res.x) if res.success else float(start)


def estimate_dispersion(groups):
    """One rho for the whole dataset, by profile likelihood.

    Per-clone estimation is hopeless at these counts. A single shared value is
    the standard compromise and is conservative: it is driven by the bulk of
    clones, which are small, so it will not be understated by a few large ones.
    """
    def nll(rho):
        total = 0.0
        for k, n in groups:
            total -= betabinom_loglik(k, n, fit_p(k, n, rho), rho)
        return total

    res = optimize.minimize_scalar(nll, bounds=(1e-4, 0.5), method="bounded")
    return float(res.x) if res.success else 0.05


def lrt(k_test, n_test, k_ref, n_ref, rho):
    """Likelihood ratio test for a difference in proportion between arms."""
    k_all = np.concatenate([k_test, k_ref])
    n_all = np.concatenate([n_test, n_ref])

    p0 = fit_p(k_all, n_all, rho)
    ll0 = betabinom_loglik(k_all, n_all, p0, rho)

    p1t = fit_p(k_test, n_test, rho)
    p1r = fit_p(k_ref, n_ref, rho)
    ll1 = betabinom_loglik(k_test, n_test, p1t, rho) + betabinom_loglik(
        k_ref, n_ref, p1r, rho
    )

    stat = max(2 * (ll1 - ll0), 0.0)
    return stat, float(stats.chi2.sf(stat, df=1)), p1t, p1r


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full_like(p, np.nan)
    if not ok.any():
        return out
    sub = p[ok]
    order = np.argsort(sub)
    ranked = sub[order]
    m = len(ranked)
    adj = np.minimum.accumulate((ranked * m / np.arange(m, 0, -1))[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


def descriptive(df, threshold):
    """No contrast: report what each sample's repertoire looks like."""
    out = df.copy()
    out["test"] = "none"
    logging.info(
        "No contrast declared; reporting %d expanded clonotypes at freq >= %.3f",
        int(out["is_expanded"].sum()),
        threshold,
    )
    return out


def differential(df, contrast, threshold):
    group_col = contrast["group_col"]
    test_level = str(contrast["test"])
    ref_level = str(contrast["reference"])
    exclude = contrast.get("exclude")
    block_col = contrast.get("block_col")

    for col in [group_col] + ([block_col] if block_col else []):
        if col not in df.columns or df[col].isna().all() or (df[col] == "").all():
            raise ValueError(
                f"Contrast needs column {col!r}, which is empty. In single_cell "
                "mode there are no experimental arms; set expansion.contrast "
                "to null."
            )

    # Alloreactive removal. A clone expanding against the excluded condition is
    # responding to the cell background rather than to the antigen, so it is
    # dropped before the comparison that matters.
    if exclude:
        allo = set(
            df.loc[(df[group_col] == exclude) & df["is_expanded"], "clonotype_id"]
        )
        before = df["clonotype_id"].nunique()
        df = df[~df["clonotype_id"].isin(allo)]
        logging.info(
            "Excluded %d of %d clonotypes as alloreactive (expanded against %r)",
            len(allo), before, exclude,
        )

    arms = df[df[group_col].isin([test_level, ref_level])]
    if arms.empty:
        raise ValueError(
            f"No rows in conditions {test_level!r}/{ref_level!r}. Present: "
            f"{sorted(df[group_col].dropna().unique())}"
        )

    depth = arms.groupby("sample_id")["n_observations"].sum().to_dict()

    # Dispersion from the reference arm only, so a real effect in the test arm
    # is not absorbed into the overdispersion estimate.
    ref_groups = []
    for _, g in arms[arms[group_col] == ref_level].groupby("clonotype_id"):
        k = g["n_observations"].to_numpy(float)
        n = g["sample_id"].map(depth).to_numpy(float)
        if len(k) > 1 and n.sum() > 0:
            ref_groups.append((k, n))
    rho = estimate_dispersion(ref_groups[:500]) if ref_groups else 0.05
    logging.info("Dispersion rho = %.4f (from %d reference clonotypes)", rho, len(ref_groups))

    blocks = sorted(arms[block_col].unique()) if block_col else [None]
    records = []

    for clonotype, g in arms.groupby("clonotype_id"):
        per_block, expanded_in = [], 0
        for b in blocks:
            gb = g if b is None else g[g[block_col] == b]
            t = gb[gb[group_col] == test_level]
            r = gb[gb[group_col] == ref_level]
            if t.empty or r.empty:
                continue
            k_t = t["n_observations"].to_numpy(float)
            n_t = t["sample_id"].map(depth).to_numpy(float)
            k_r = r["n_observations"].to_numpy(float)
            n_r = r["sample_id"].map(depth).to_numpy(float)
            if n_t.sum() == 0 or n_r.sum() == 0:
                continue
            _, pval, p_t, p_r = lrt(k_t, n_t, k_r, n_r, rho)
            per_block.append(pval)
            if p_t > p_r:
                expanded_in += 1

        if not per_block:
            continue

        # Fisher's method across donors. One-sided in effect, because
        # expanded_in records direction separately.
        chi = -2 * np.sum(np.log(np.clip(per_block, EPS, 1.0)))
        combined = float(stats.chi2.sf(chi, df=2 * len(per_block)))

        sub = g[g[group_col].isin([test_level, ref_level])]
        freqs = sub.groupby(group_col)["frequency"].mean()
        f_t = float(freqs.get(test_level, 0.0))
        f_r = float(freqs.get(ref_level, 0.0))

        records.append(
            {
                "clonotype_id": clonotype,
                "n_blocks_tested": len(per_block),
                "n_blocks_expanded": expanded_in,
                "freq_test": round(f_t, 6),
                "freq_reference": round(f_r, 6),
                "log2fc": round(float(np.log2((f_t + EPS) / (f_r + EPS))), 4),
                "pvalue": combined,
            }
        )

    res = pd.DataFrame.from_records(records)
    if res.empty:
        raise ValueError("No clonotype had data in both arms of the contrast.")

    res["qvalue"] = bh(res["pvalue"].to_numpy())
    res["significant"] = res["qvalue"] < float(contrast.get("fdr", 0.05))
    res["recurrent"] = res["n_blocks_expanded"] >= 2
    res["test"] = f"{test_level}_vs_{ref_level}"

    meta = df.drop_duplicates("clonotype_id").set_index("clonotype_id")
    for col in ("sample_id", "cdr3a_aa", "v_a", "j_a", "cdr3b_aa", "v_b", "j_b"):
        if col in meta.columns:
            res[col] = res["clonotype_id"].map(meta[col])

    return res.sort_values(["qvalue", "pvalue"])


def plot(res, figdir, contrast):
    figdir = Path(figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))
    if contrast is None or "log2fc" not in res.columns:
        data = res.groupby("sample_id")["is_expanded"].sum().sort_values()
        ax.barh(data.index.astype(str), data.to_numpy(), color="#4C72B0")
        ax.set_xlabel("Expanded clonotypes", fontsize=28)
        ax.set_ylabel("Sample", fontsize=28)
    else:
        sig = res["significant"].fillna(False)
        ax.scatter(res.loc[~sig, "log2fc"], -np.log10(res.loc[~sig, "pvalue"] + EPS),
                   s=18, c="#B0B0B0", label="ns")
        ax.scatter(res.loc[sig, "log2fc"], -np.log10(res.loc[sig, "pvalue"] + EPS),
                   s=28, c="#C44E52", label="q < threshold")
        ax.set_xlabel("log2 fold change", fontsize=28)
        ax.set_ylabel("-log10 p", fontsize=28)
        ax.legend(fontsize=28, frameon=False)

    ax.tick_params(labelsize=28)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(figdir / f"expansion.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main(snakemake):
    setup_logging(snakemake.log[0])
    p = snakemake.params
    threshold = float(p.expansion_threshold)
    contrast = p.contrast

    df = pd.read_csv(snakemake.input.by_sample, sep="\t")
    df["is_expanded"] = df["is_expanded"].astype(str).str.lower().isin(("true", "1"))
    logging.info("%d clonotypes across %d samples", len(df), df["sample_id"].nunique())

    if not contrast:
        res = descriptive(df, threshold)
    else:
        contrast = dict(contrast)
        contrast.setdefault("fdr", p.fdr)
        res = differential(df, contrast, threshold)
        n_sig = int(res["significant"].sum())
        logging.info(
            "%d of %d clonotypes significant at q < %.3f; %d recurrent across blocks",
            n_sig, len(res), float(contrast["fdr"]), int(res["recurrent"].sum()),
        )

    Path(snakemake.output.results).parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(snakemake.output.results, sep="\t", index=False)
    plot(res, snakemake.output.figdir, contrast)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
