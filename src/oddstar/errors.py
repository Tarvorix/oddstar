"""Structured error definitions and stable exit codes for oddstar CLI & MCP server."""

from __future__ import annotations


EXIT_OK = 0
EXIT_INTERNAL = 1
EXIT_CONFIG = 2
EXIT_UNCERTIFIED = 3
EXIT_ARCHIVE = 4
EXIT_VALIDATION = 5
EXIT_NOT_FOUND = 6
EXIT_POLICY = 7
EXIT_ISOLATION = 8


class OddstarError(Exception):
    """Base error carrying structured CLI contract metadata."""

    code: str = "INTERNAL_ERROR"
    exit_code: int = EXIT_INTERNAL
    retryable: bool = False

    def __init__(
        self,
        message: str,
        *,
        hint: str = "",
        code: str | None = None,
        exit_code: int | None = None,
        retryable: bool | None = None,
    ) -> None:
        super().__init__(message)
        self.message = message
        self.hint = hint
        if code is not None:
            self.code = code
        if exit_code is not None:
            self.exit_code = exit_code
        if retryable is not None:
            self.retryable = retryable

    def to_dict(self) -> dict[str, object]:
        return {
            "ok": False,
            "code": self.code,
            "message": self.message,
            "hint": self.hint,
            "retryable": self.retryable,
        }


class ConfigError(OddstarError):
    code = "CONFIG_ERROR"
    exit_code = EXIT_CONFIG


class UncertifiedSearchError(OddstarError):
    code = "UNCERTIFIED_CHANNEL"
    exit_code = EXIT_UNCERTIFIED


class ArchiveAccessError(OddstarError):
    code = "ARCHIVE_ERROR"
    exit_code = EXIT_ARCHIVE
    retryable = True


class ValidationGateError(OddstarError):
    code = "VALIDATION_FAILED"
    exit_code = EXIT_VALIDATION


class ArtifactNotFoundError(OddstarError):
    code = "NOT_FOUND"
    exit_code = EXIT_NOT_FOUND


class PolicyViolationError(OddstarError):
    code = "POLICY_VIOLATION"
    exit_code = EXIT_POLICY


class IsolationError(OddstarError):
    code = "ISOLATION_VIOLATION"
    exit_code = EXIT_ISOLATION
