"""Command-line interface.

Exposes the stages that are ported and runnable today. Each command works
against the existing R pipeline's output files, so you can reproduce the
corrections without refetching anything.

    resko rank results/A12B_compound_protein_evidence_corrected.csv
    resko audit-cache results/A17D_api_cache
    resko coverage eef1a_master_seed_list.csv SIDER/drug_names.tsv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer

from . import __version__
from .chembl import audit_cache_dir
from .ranking import ScoringMode, rank_candidates
from .sider import coverage_report, load_drug_names, match_drug_names

app = typer.Typer(add_completion=False, help=__doc__)


@app.command()
def version() -> None:
    """Print the package version."""
    typer.echo(__version__)


@app.command()
def rank(
    evidence_csv: Path = typer.Argument(..., help="Compound-protein evidence table"),
    out: Path = typer.Option(None, "--out", "-o", help="Write ranking CSV here"),
    legacy: bool = typer.Option(
        False, "--legacy", help="Reproduce the R original (relation ignored)"
    ),
    compare: bool = typer.Option(
        False, "--compare", help="Show both rankings side by side"
    ),
) -> None:
    """Rank candidates from an evidence table.

    By default applies the relation-aware correction (CODE_REVIEW.md S1-1).
    Use --compare to see what the correction changes.
    """
    evidence = pd.read_csv(evidence_csv, low_memory=False)
    mode = ScoringMode.legacy() if legacy else ScoringMode()
    ranked = rank_candidates(evidence, mode)

    if compare:
        other = rank_candidates(evidence, ScoringMode.legacy() if not legacy else ScoringMode())
        merged = ranked[["molecule_chembl_id", "display_name", "rank"]].merge(
            other[["molecule_chembl_id", "rank"]],
            on="molecule_chembl_id",
            suffixes=("_this", "_other"),
        )
        typer.echo(merged.to_string(index=False))
    else:
        columns = [
            "rank", "molecule_chembl_id", "display_name", "best_protein",
            "best_record_score", "n_proteins", "n_censored",
            "n_proteomics", "all_evidence_censored",
        ]
        typer.echo(ranked[columns].to_string(index=False))

    if out:
        ranked.to_csv(out, index=False)
        typer.echo(f"\nwrote {out}", err=True)


@app.command("audit-cache")
def audit_cache(
    cache_dir: Path = typer.Argument(..., help="ChEMBL response cache directory"),
    show_all: bool = typer.Option(False, "--all", help="List every file, not just strays"),
) -> None:
    """Report cache files that do not match the canonical key convention.

    Detects the orphaned-file / schema-divergence problem in
    CODE_REVIEW.md S2-1.
    """
    audit = audit_cache_dir(cache_dir)
    paged = audit[audit["paged"]].copy()
    paged["canonical"] = paged["canonical_name"].astype(bool)

    typer.echo(
        f"{len(audit)} files: {int(paged['canonical'].sum())} canonical, "
        f"{int((~paged['canonical']).sum())} orphaned, "
        f"{int((~audit['paged']).sum())} un-paged"
    )
    if len(paged):
        typer.echo("\nkeys per record, by convention:")
        typer.echo(
            paged.groupby("canonical")[["n_fields_first", "n_fields_union"]]
            .agg(["min", "max"])
            .to_string()
        )
    strays = paged[~paged["canonical"]]
    if len(strays):
        typer.secho(
            f"\n{len(strays)} orphaned file(s) no current code path reads:",
            fg=typer.colors.YELLOW,
        )
        typer.echo(strays[["file", "n_records", "n_fields_union"]].to_string(index=False))
    if show_all:
        typer.echo("\n" + audit.to_string(index=False))


@app.command()
def coverage(
    seed_csv: Path = typer.Argument(..., help="Seed list CSV (one compound per row)"),
    drug_names_tsv: Path = typer.Argument(..., help="SIDER drug_names.tsv"),
    seed_column: str = typer.Option("Drug", "--column", help="Seed name column"),
) -> None:
    """Report SIDER coverage of a seed list with boundary-anchored matching.

    Category labels ("Didemnin analogues") are excluded from the denominator
    since they can never match a database row (CODE_REVIEW.md S3-1).
    """
    seeds_frame = pd.read_csv(seed_csv)
    if seed_column not in seeds_frame.columns:
        raise typer.BadParameter(
            f"column {seed_column!r} not in {list(seeds_frame.columns)}"
        )
    seeds = seeds_frame[seed_column].dropna().astype(str).tolist()
    matches = match_drug_names(seeds, load_drug_names(drug_names_tsv))
    report = coverage_report(matches)

    typer.echo(f"seeds supplied:        {report['n_seeds_supplied']}")
    typer.echo(f"category labels:      -{report['n_category_labels_excluded']}")
    typer.echo(f"compounds searched:    {report['n_compounds_searched']}")
    typer.echo(f"matched in SIDER:      {report['n_matched']}")
    typer.echo(f"coverage:              {report['coverage_fraction']:.1%}")
    if report["matched_names"]:
        typer.echo("matched: " + ", ".join(report["matched_names"]))


if __name__ == "__main__":
    app()
