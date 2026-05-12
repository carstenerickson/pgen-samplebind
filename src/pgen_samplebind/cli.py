"""Click entry point. Per LLD §3.16.

main() runs click in non-standalone mode so the LLD owns exit-code routing
explicitly rather than inheriting click's defaults (click.UsageError defaults
to exit 2, which collides with HLD's exit-2 reservation for I/O failure).
"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

import click

from . import __version__
from .commands.afs_cmd import run_afs
from .commands.hash_cmd import run_hash
from .commands.inspect_cmd import run_inspect
from .commands.merge_cmd import run_merge
from .commands.validate_cmd import run_validate
from .errors import PgenSamplebindError
from .types import ExitCode, MergePolicy


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(version=__version__, prog_name="pgen-samplebind")
def cli() -> None:
    """pgen-samplebind: bind PFILE/BFILE/EIGENSTRAT genotype datasets sharing a variant set.

    The missing plink2 --pmerge case for ancient-DNA / population-genetics workflows.
    """


@cli.command("hash")
@click.argument("input_path", type=click.Path(exists=False, path_type=Path))
@click.option(
    "--emit-canonical",
    is_flag=True,
    default=False,
    help="Print the canonicalized bytestream that would be hashed (for diagnosis).",
)
def hash_command(input_path: Path, emit_canonical: bool) -> None:
    """Emit canonical variant-set hash for cross-format panel-identity verification."""
    run_hash(input_path, emit_canonical)


@cli.command("inspect")
@click.argument("input_path", type=click.Path(exists=False, path_type=Path))
@click.option(
    "--json",
    "json_output",
    is_flag=True,
    default=False,
    help="Emit machine-readable JSON instead of TSV-style text.",
)
def inspect_command(input_path: Path, json_output: bool) -> None:
    """Structured summary of one input: format, samples, variants, populations, sex distribution."""
    run_inspect(input_path, json_output)


@cli.command("merge")
@click.argument("inputs", nargs=-1, required=True, type=click.Path(path_type=Path))
@click.option(
    "--target",
    "target_paths",
    type=click.Path(path_type=Path),
    multiple=True,
    help="Mark one or more inputs as targets (single-sample / small-cohort mode); "
    "activates asymmetric strand-check + per-target call-rate gate "
    "(--target-min-call-rate). Repeat --target to append multiple targets. "
    "Targets are appended after the positional inputs; canonical remains input[0]. "
    "Under --on-collision suffix, a single target uses `_target`; multiple targets "
    "use `_target_<input_idx>` to disambiguate.",
)
@click.option(
    "-o",
    "--out",
    "output_prefix",
    type=click.Path(path_type=Path),
    required=True,
    help="Output PFILE prefix.",
)
@click.option("--variant-key", type=click.Choice(["chr_pos", "id"]), default="chr_pos")
@click.option("--on-mismatch", type=click.Choice(["drop", "error"]), default="drop")
@click.option(
    "--on-missing",
    type=click.Choice(["fill_missing", "drop_variant", "error"]),
    default="fill_missing",
)
@click.option("--on-extra", type=click.Choice(["warn", "drop", "error"]), default="warn")
@click.option("--on-strand", type=click.Choice(["drop", "flip", "error"]), default="flip")
@click.option("--trust-strand", is_flag=True, default=False)
@click.option(
    "--on-collision",
    type=click.Choice(["error", "first", "suffix"]),
    default="error",
)
@click.option("--id-column", default="IID", help=".psam column for identity ops.")
@click.option(
    "--population-column",
    default=None,
    help=".psam column holding population labels (default: POP/PHENO/PHENO1 fallback).",
)
@click.option("--target-min-call-rate", type=float, default=0.40)
@click.option(
    "--validate-strand-fail-pct",
    type=float,
    default=10.0,
    help="Exit 1 if ambiguous-strand drops exceed N%% of intersection (default 10).",
)
@click.option(
    "--relabel-from",
    "relabel_from",
    type=click.Path(path_type=Path),
    default=None,
    help="TSV-driven relabel of the population (POP) column. Two-col form "
    "(no header) maps POP→POP across inputs; N-col form (with --relabel-input-col "
    "and --relabel-output-col) joins per-sample on the input column.",
)
@click.option(
    "--relabel-input-col",
    default=None,
    help="For N-col --relabel-from TSVs: which column matches each sample's "
    "id_column value (default IID; override with --id-column).",
)
@click.option(
    "--relabel-output-col",
    default=None,
    help="For N-col --relabel-from TSVs: which column to map to (becomes the new POP value).",
)
@click.option(
    "--report",
    "report_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Per-variant action TSV (streamed; constant memory).",
)
@click.option(
    "--report-json",
    "report_json_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Run-level summary JSON (~few KB; rows excluded by default).",
)
@click.option(
    "--report-json-include-rows",
    is_flag=True,
    default=False,
    help="Include per-variant rows in --report-json output (buffered; warns at "
    ">100 MB predicted size; prefer --report TSV for streaming at scale).",
)
@click.option("--quiet", is_flag=True, default=False, help="Suppress progress to stdout.")
@click.option(
    "--block-size",
    type=int,
    default=2048,
    help="Variants per pgenlib read block (default: 2048).",
)
def merge_command(
    inputs: tuple[Path, ...],
    target_paths: tuple[Path, ...],
    output_prefix: Path,
    variant_key: str,
    on_mismatch: str,
    on_missing: str,
    on_extra: str,
    on_strand: str,
    trust_strand: bool,
    on_collision: str,
    id_column: str,
    population_column: str | None,
    target_min_call_rate: float,
    validate_strand_fail_pct: float,
    relabel_from: Path | None,
    relabel_input_col: str | None,
    relabel_output_col: str | None,
    report_path: Path | None,
    report_json_path: Path | None,
    report_json_include_rows: bool,
    quiet: bool,
    block_size: int,
) -> None:
    """Bind inputs into one output PFILE."""
    policy = MergePolicy(
        on_mismatch=on_mismatch,  # type: ignore[arg-type]
        on_missing=on_missing,  # type: ignore[arg-type]
        on_extra=on_extra,  # type: ignore[arg-type]
        on_strand=on_strand,  # type: ignore[arg-type]
        on_collision=on_collision,  # type: ignore[arg-type]
        trust_strand=trust_strand,
        variant_key=variant_key,  # type: ignore[arg-type]
        target_min_call_rate=target_min_call_rate,
        validate_strand_fail_pct=validate_strand_fail_pct,
        population_column=population_column,
        id_column=id_column,
        block_size=block_size,
        report_json_include_rows=report_json_include_rows,
    )
    run_merge(
        input_paths=inputs,
        target_paths=target_paths,
        output_prefix=output_prefix,
        policy=policy,
        report_path=report_path,
        report_json_path=report_json_path,
        quiet=quiet,
        relabel_from=relabel_from,
        relabel_input_col=relabel_input_col,
        relabel_output_col=relabel_output_col,
    )


@cli.command("validate")
@click.argument("inputs", nargs=-1, required=True, type=click.Path(path_type=Path))
@click.option("--variant-key", type=click.Choice(["chr_pos", "id"]), default="chr_pos")
@click.option("--on-mismatch", type=click.Choice(["drop", "error"]), default="drop")
@click.option(
    "--on-missing",
    type=click.Choice(["fill_missing", "drop_variant", "error"]),
    default="fill_missing",
)
@click.option("--on-extra", type=click.Choice(["warn", "drop", "error"]), default="warn")
@click.option("--on-strand", type=click.Choice(["drop", "flip", "error"]), default="flip")
@click.option("--trust-strand", is_flag=True, default=False)
@click.option(
    "--on-collision",
    type=click.Choice(["error", "first", "suffix"]),
    default="error",
)
@click.option("--id-column", default="IID")
@click.option("--population-column", default=None)
@click.option(
    "--validate-strand-fail-pct",
    type=float,
    default=10.0,
    help="Exit 1 if ambiguous-strand drops exceed N%% of intersection (default 10).",
)
@click.option(
    "--relabel-from",
    "relabel_from",
    type=click.Path(path_type=Path),
    default=None,
    help="TSV-driven POP relabel (per HLD §Relabeling). 2-col header-less or "
    "N-col with --relabel-input-col / --relabel-output-col.",
)
@click.option("--relabel-input-col", default=None)
@click.option("--relabel-output-col", default=None)
@click.option(
    "--report",
    "report_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Per-variant action TSV.",
)
@click.option(
    "--report-json",
    "report_json_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Run-level summary JSON.",
)
@click.option(
    "--report-json-include-rows",
    is_flag=True,
    default=False,
    help="Include per-variant rows in --report-json output.",
)
@click.option("--quiet", is_flag=True, default=False)
def validate_command(
    inputs: tuple[Path, ...],
    variant_key: str,
    on_mismatch: str,
    on_missing: str,
    on_extra: str,
    on_strand: str,
    trust_strand: bool,
    on_collision: str,
    id_column: str,
    population_column: str | None,
    validate_strand_fail_pct: float,
    relabel_from: Path | None,
    relabel_input_col: str | None,
    relabel_output_col: str | None,
    report_path: Path | None,
    report_json_path: Path | None,
    report_json_include_rows: bool,
    quiet: bool,
) -> None:
    """Check alignment of inputs without writing output. Exits 0 if alignment
    OK, 1 if any Exit-1 gate fires (HLD §Exit-1 validation gates), 3 on
    invariant violation (multi-allelic input, duplicate canonical keys,
    --on-collision error)."""
    policy = MergePolicy(
        on_mismatch=on_mismatch,  # type: ignore[arg-type]
        on_missing=on_missing,  # type: ignore[arg-type]
        on_extra=on_extra,  # type: ignore[arg-type]
        on_strand=on_strand,  # type: ignore[arg-type]
        on_collision=on_collision,  # type: ignore[arg-type]
        trust_strand=trust_strand,
        variant_key=variant_key,  # type: ignore[arg-type]
        validate_strand_fail_pct=validate_strand_fail_pct,
        population_column=population_column,
        id_column=id_column,
        report_json_include_rows=report_json_include_rows,
    )
    run_validate(
        input_paths=inputs,
        policy=policy,
        report_path=report_path,
        report_json_path=report_json_path,
        quiet=quiet,
        relabel_from=relabel_from,
        relabel_input_col=relabel_input_col,
        relabel_output_col=relabel_output_col,
    )


@cli.command("afs")
@click.argument("input_path", type=click.Path(exists=False, path_type=Path))
@click.option(
    "-o",
    "--out",
    "output_dir",
    type=click.Path(path_type=Path),
    required=True,
    help="Output directory for the three AFS TSVs + manifest.json.",
)
@click.option(
    "--population-column",
    default="POP",
    help="`.psam` column to aggregate by (default: POP — set by pgen-samplebind merge).",
)
@click.option(
    "--populations",
    "populations_filter",
    multiple=True,
    help="Restrict to a subset of populations (repeat flag for each). Default: all.",
)
@click.option(
    "--no-pseudohaploid-adjust",
    "adjust_pseudohaploid",
    is_flag=True,
    default=True,
    flag_value=False,
    help="Skip pseudohaploid called-allele adjustment (treat all samples as diploid). "
    "Default: adjust when input has a PSEUDOHAPLOID column.",
)
@click.option(
    "--include-sex-chrom",
    is_flag=True,
    default=False,
    help="Include sex chromosomes (chr 23/24/25/26). Default: autosomes only.",
)
@click.option("--block-size", type=int, default=2048)
@click.option("--quiet", is_flag=True, default=False)
def afs_command(
    input_path: Path,
    output_dir: Path,
    population_column: str,
    populations_filter: tuple[str, ...],
    adjust_pseudohaploid: bool,
    include_sex_chrom: bool,
    block_size: int,
    quiet: bool,
) -> None:
    """Per-population allele-frequency-spectrum TSVs from one PFILE/BFILE/EIGENSTRAT input.

    Emits three TSVs + manifest matching AT2's `*_to_afs()` shape so a small R
    loader can feed them into downstream f-statistic work without a BED
    intermediate. Bridge until `pfile_to_afs()` lands in AT2 upstream.
    """
    run_afs(
        input_path=input_path,
        output_dir=output_dir,
        population_column=population_column,
        populations_filter=populations_filter if populations_filter else None,
        adjust_pseudohaploid=adjust_pseudohaploid,
        include_sex_chrom=include_sex_chrom,
        block_size=block_size,
        quiet=quiet,
    )


def main() -> None:
    """Console-script entry point.

    Runs click in non-standalone mode and routes every exit through the
    HLD-pinned ExitCode enum.
    """
    try:
        cli(standalone_mode=False)
    except click.UsageError as e:
        click.echo(f"usage error: {e.format_message()}", err=True)
        sys.exit(ExitCode.USAGE_ERROR)
    except click.exceptions.Abort:
        sys.exit(ExitCode.USAGE_ERROR)
    except PgenSamplebindError as e:
        click.echo(f"error: {e}", err=True)
        sys.exit(e.exit_code)
    except NotImplementedError as e:
        # Day 1 stubs for merge/validate.
        click.echo(f"not implemented: {e}", err=True)
        sys.exit(ExitCode.USAGE_ERROR)
    except Exception:
        traceback.print_exc()
        sys.exit(ExitCode.INVARIANT_VIOLATION)
