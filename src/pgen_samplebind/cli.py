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
from .commands.hash_cmd import run_hash
from .commands.inspect_cmd import run_inspect
from .errors import PgenSamplebindError
from .types import ExitCode


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
def merge_command() -> None:
    """Bind inputs into one output PFILE. (Deferred to project Day 3 per HLD.)"""
    raise NotImplementedError("merge subcommand deferred to project Day 3 per HLD §Project plan")


@cli.command("validate")
def validate_command() -> None:
    """Check alignment, no output written. (Deferred to project Day 4 per HLD.)"""
    raise NotImplementedError("validate subcommand deferred to project Day 4 per HLD §Project plan")


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
