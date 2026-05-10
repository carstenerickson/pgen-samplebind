"""Exception hierarchy. cli.main() catches PgenSamplebindError and exits with `e.exit_code`."""

from __future__ import annotations

from .types import ExitCode


class PgenSamplebindError(Exception):
    """Base for all tool-internal errors."""

    exit_code: ExitCode = ExitCode.INVARIANT_VIOLATION


class ValidationError(PgenSamplebindError):
    """Alignment failure, call-rate gate breach, threshold violation."""

    exit_code = ExitCode.VALIDATION_FAILURE


class IOFailure(PgenSamplebindError):
    """Cannot read input, cannot write output, advisory lock held, plink2 subprocess failed."""

    exit_code = ExitCode.IO_FAILURE


class InvariantViolation(PgenSamplebindError):
    """Encoding mismatch beyond resolution; multi-allelic input; duplicate canonical keys;
    --on-* error policy triggered."""

    exit_code = ExitCode.INVARIANT_VIOLATION


class UsageError(PgenSamplebindError):
    """Bad CLI argument combination not catchable by click's own validation."""

    exit_code = ExitCode.USAGE_ERROR
