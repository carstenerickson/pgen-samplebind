"""Smoke test for `hash` subcommand. Full HLD test 10 (cross-format invariance)
deferred until BFILE/EIGENSTRAT inputs are wired (project Day 6).
"""

from __future__ import annotations

import re

from click.testing import CliRunner

from pgen_samplebind.cli import cli
from pgen_samplebind.types import InputDescriptor

_SHA256_LINE = re.compile(r"^sha256:[0-9a-f]{64}$")


def test_hash_emits_sha256_prefix(synth_panel_tiny: InputDescriptor) -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["hash", str(synth_panel_tiny.path)])
    assert result.exit_code == 0, result.output
    line = result.output.strip()
    assert _SHA256_LINE.match(line), f"unexpected output: {line!r}"


def test_hash_deterministic_for_same_input(synth_panel_tiny: InputDescriptor) -> None:
    runner = CliRunner()
    r1 = runner.invoke(cli, ["hash", str(synth_panel_tiny.path)])
    r2 = runner.invoke(cli, ["hash", str(synth_panel_tiny.path)])
    assert r1.exit_code == 0
    assert r2.exit_code == 0
    assert r1.output == r2.output


def test_emit_canonical_bytestream_is_tsv(synth_panel_tiny: InputDescriptor) -> None:
    runner = CliRunner()
    result = runner.invoke(cli, ["hash", "--emit-canonical", str(synth_panel_tiny.path)])
    assert result.exit_code == 0, result.output
    # First line should be 4 tab-separated fields: chrom\tpos\tref\talt
    first_line = result.output.split("\n", 1)[0]
    fields = first_line.split("\t")
    assert len(fields) == 4
    assert fields[0].isdigit()  # chrom is int
    assert fields[1].isdigit()  # pos is int
    assert fields[2] in {"A", "C", "G", "T"}
    assert fields[3] in {"A", "C", "G", "T"}
