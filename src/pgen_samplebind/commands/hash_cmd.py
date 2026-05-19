"""`hash` subcommand orchestrator. Sequence in"""

from __future__ import annotations

import sys
from pathlib import Path

from ..formats import prepared_input
from ..hashing import hash_input_with_canonical
from ..pvar import check_max_alleles


def run_hash(input_path: Path, emit_canonical: bool) -> None:
    """Read input, canonicalize variants, SHA-256, emit.

    --emit-canonical writes the canonicalized bytestream to stdout (no hash
    line) for diagnosis. Default writes `sha256:<hex>` followed by newline.
    """
    with prepared_input(input_path) as desc:
        check_max_alleles(desc.pgen_path)
        h, canonical = hash_input_with_canonical(desc)

    if emit_canonical:
        sys.stdout.buffer.write(canonical)
        sys.stdout.flush()
    else:
        print(h.render())
