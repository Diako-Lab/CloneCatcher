"""Barcode namespace helpers shared across CloneCatcher stages.

ClusterCatcher's canonical cell name is ``{sample_id}_{cellranger_barcode}``.
CellRanger VDJ output carries the bare ``{cellranger_barcode}`` with no sample
identity. Everything that crosses that gap goes through this module.

Two rules are enforced here rather than left to each caller:

  * A canonical name is **prefixed**, never rebuilt. The CellRanger suffix is
    part of the barcode; we neither synthesize nor strip one.
  * A failed join is *unknown* data, never *absent* data. Zero overlap raises.

This module holds no Snakemake rule body and is never run directly. It is
imported by barcode_namespace.py and join_clustercatcher.py.
"""

import re
from pathlib import Path

import pandas as pd

# {sample_id}_{barcode}-{suffix}. sample_id may itself contain underscores,
# which is why every parse below uses rsplit rather than split: the barcode
# half is ACGT plus a numeric suffix and can never contain an underscore.
CANONICAL = re.compile(r"^.+_[ACGT]+-\d+$")

# Pre-fix form from ClusterCatcher v1.0.2 through v1.4.2. Joins to nothing.
LEGACY = re.compile(r"^[ACGT]+-\d+-\S+$")

# SComatic internals and raw CellRanger output. Carries no sample identity.
BARE = re.compile(r"^[ACGT]+(-\d+)?$")

# Searched in order. The signature outputs are absent during a signature
# rerun, so fall back to the most recent AnnData that exists. Barcodes are
# set at annotation and never change downstream, so any of these answers the
# namespace question identically.
ADATA_CANDIDATES = (
    "signatures/adata_final.h5ad",
    "adata_final.h5ad",
    "dysregulation/adata_cancer_detected.h5ad",
    "viral_integration/adata_viral_integrated.h5ad",
    "viral_integration/adata_with_virus.h5ad",
    "annotation/adata_annotated.h5ad",
)


class NamespaceError(ValueError):
    """Raised when barcodes are not in the namespace a stage requires."""


def classify(barcode):
    """Name the namespace a single barcode belongs to."""
    barcode = str(barcode)
    if CANONICAL.match(barcode):
        return "canonical"
    if LEGACY.match(barcode):
        return "legacy"
    if BARE.match(barcode):
        return "bare"
    return "unrecognized"


def split_canonical(name):
    """Split a canonical name into (sample_id, barcode).

    Uses rsplit so a sample_id containing underscores still parses. The
    barcode half is ACGT plus a numeric suffix and never contains one.
    """
    name = str(name)
    if "_" not in name:
        raise NamespaceError(f"Not a canonical cell name: {name!r}")
    return tuple(name.rsplit("_", 1))


def assert_canonical(barcodes, source="input", probe_n=1000):
    """Raise unless barcodes are in ClusterCatcher canonical form.

    Checks a probe rather than the full set; a namespace error is uniform
    across a file in every case we have seen, and the full scan is wasteful
    on a 170k-row frame.
    """
    s = pd.Index(barcodes).astype(str).to_series()
    probe = s.head(probe_n)
    if probe.empty:
        raise NamespaceError(f"{source}: no barcodes to check.")

    legacy = probe[probe.str.match(LEGACY)]
    if len(legacy):
        raise NamespaceError(
            f"{source}: barcodes are in the pre-fix "
            f"'{{barcode}}-{{suffix}}-{{sample_id}}' form "
            f"(e.g. {legacy.iloc[0]!r}). This ClusterCatcher run predates the "
            "namespace fix. Run repair_barcodes.py on it first."
        )

    if probe.str.match(BARE).all():
        raise NamespaceError(
            f"{source}: barcodes are bare and carry no sample identity "
            f"(e.g. {probe.iloc[0]!r}). This is a per-sample intermediate, "
            "not a collected output."
        )

    bad = probe[~probe.str.match(CANONICAL)]
    if len(bad):
        raise NamespaceError(
            f"{source}: {len(bad)}/{len(probe)} barcodes are not canonical "
            f"(e.g. {bad.iloc[0]!r}); expected "
            "'{sample_id}_{barcode}-{suffix}'."
        )
    return True


def assert_join_rate(left, right, min_rate=0.40, label="join"):
    """Raise when two barcode sets barely overlap.

    This is the tripwire. A low rate on a correct join is a finding worth
    inspecting; a zero rate is a bug by definition, because two tables drawn
    from the same run always share cells.

    Returns (overlap, rate).
    """
    l = set(map(str, left))
    r = set(map(str, right))
    overlap = len(l & r)
    rate = overlap / max(len(l), 1)

    if rate < min_rate:
        l_ex = sorted(l)[0] if l else "(empty)"
        r_ex = sorted(r)[0] if r else "(empty)"
        detail = (
            "Zero overlap is a namespace mismatch, not missing data."
            if overlap == 0
            else f"Below the configured floor of {min_rate:.0%}."
        )
        raise NamespaceError(
            f"{label}: only {overlap:,}/{len(l):,} ({rate:.1%}) of barcodes "
            f"matched. {detail}\n"
            f"  left  example: {l_ex!r}\n"
            f"  right example: {r_ex!r}"
        )
    return overlap, rate


def promote(bare_barcodes, sample_id):
    """Lift bare CellRanger barcodes into the canonical namespace.

    Prefixes sample_id onto the barcode CellRanger wrote. Does not touch the
    suffix: appending a hardcoded '-1' is precisely how the original
    ClusterCatcher namespace bug was introduced.
    """
    sample_id = str(sample_id)
    if "_" in sample_id:
        # Not fatal -- rsplit still parses it -- but worth surfacing, since
        # any consumer using the documented split('_', 1) will mis-parse.
        import warnings

        warnings.warn(
            f"sample_id {sample_id!r} contains an underscore. Parse canonical "
            "names with rsplit('_', 1), not split('_', 1).",
            stacklevel=2,
        )
    return [f"{sample_id}_{bc}" for bc in map(str, bare_barcodes)]


def build_barcode_map(adata, sample_col="sample_id"):
    """Return {sample_id: {bare_barcode: canonical_name}} from an AnnData.

    Preferred over promote() wherever the cell exists in the GEX matrix,
    because it looks the canonical name up rather than reconstructing it.
    Requires the 'original_barcode' bridge column that ClusterCatcher writes
    on the same line that builds the canonical name.
    """
    if "original_barcode" not in adata.obs.columns:
        raise NamespaceError(
            "AnnData has no 'original_barcode' column, so bare barcodes "
            "cannot be looked up. This run predates the bridge column."
        )
    out = {}
    for s, o, idx in zip(
        adata.obs[sample_col].astype(str),
        adata.obs["original_barcode"].astype(str),
        adata.obs_names.astype(str),
    ):
        out.setdefault(s, {})[o] = idx
    return out


def find_clustercatcher_adata(run_dir):
    """Locate the most recent usable ClusterCatcher AnnData."""
    run_dir = Path(run_dir)
    for rel in ADATA_CANDIDATES:
        path = run_dir / rel
        if path.exists():
            return path
    raise FileNotFoundError(
        f"No ClusterCatcher AnnData under {run_dir}. Tried: "
        + ", ".join(ADATA_CANDIDATES)
    )


def load_clustercatcher_obs(run_dir, columns=None):
    """Read just the obs frame from a ClusterCatcher run.

    Reads backed so a large h5ad does not have to come into memory for what
    is only ever a barcode and cell-type lookup.
    """
    import anndata as ad

    path = find_clustercatcher_adata(run_dir)
    adata = ad.read_h5ad(path, backed="r")
    obs = adata.obs.copy()
    obs.index = adata.obs_names.astype(str)
    if columns is not None:
        keep = [c for c in columns if c in obs.columns]
        obs = obs[keep]
    return obs, path
