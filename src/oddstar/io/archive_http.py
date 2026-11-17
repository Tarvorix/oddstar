"""Shared HTTP helpers for vetting query modules: bounded retries and content-addressed RawCache provenance.

Same house style as `oddstar.io.optical_timeseries` (retries with exponential backoff on timeouts, connection
errors, HTTP 429 and 5xx; one cache blob per request signature), kept separate so the infrared vetting queries
do not depend on private helpers of the light-curve module.
"""

from __future__ import annotations

import hashlib
import json
import time
from collections.abc import Callable
from typing import Any

import requests

from oddstar.errors import ArchiveAccessError
from oddstar.io.cache import RawCache
from oddstar.schemas import ManifestRecord

Sleeper = Callable[[float], None]


def response_bytes(resp: Any) -> bytes:
    """Return the body of a `requests.Response` (or a test double exposing only `.text`)."""
    content = getattr(resp, "content", None)
    if isinstance(content, (bytes, bytearray)):
        return bytes(content)
    return str(getattr(resp, "text", "") or "").encode("utf-8")


def get_with_retries(
    url: str,
    params: dict[str, Any] | None,
    *,
    timeout_s: float,
    attempts: int = 3,
    backoff_s: float = 2.0,
    sleep: Sleeper = time.sleep,
) -> Any:
    """GET with bounded retries; raises ArchiveAccessError after `attempts` or on a non-retryable HTTP status."""
    errors: list[str] = []
    for i in range(attempts):
        try:
            resp = requests.get(url, params=params, timeout=timeout_s)
            if resp.status_code == 200:
                return resp
            errors.append(f"HTTP {resp.status_code}")
            if resp.status_code < 500 and resp.status_code != 429:
                break
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
            errors.append(type(exc).__name__)
        if i < attempts - 1:
            sleep(backoff_s * (2**i))
    raise ArchiveAccessError(f"GET {url} failed after {len(errors)} attempt(s): {'; '.join(errors)}", retryable=True)


def request_signature(url: str, params: dict[str, Any]) -> tuple[str, str]:
    desc = f"GET {url} {json.dumps(params, sort_keys=True)}"
    return desc, hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]


def cached_get(
    *,
    url: str,
    params: dict[str, Any],
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: float,
    attempts: int = 3,
    sleep: Sleeper = time.sleep,
    validate: Callable[[bytes], None] | None = None,
) -> tuple[bytes, ManifestRecord, bool]:
    """Fetch raw bytes once into the content-addressed cache.

    `validate` is called on fresh payloads before caching; it raises to reject an error page so that a bad
    response is never cached as data.
    """
    desc, sig = request_signature(url, params)
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return cache.read_bytes(cached), cached, True
        except FileNotFoundError:
            pass  # evicted blob: refetch
    resp = get_with_retries(url, params, timeout_s=timeout_s, attempts=attempts, sleep=sleep)
    payload = response_bytes(resp)
    if validate is not None:
        validate(payload)
    man = cache.put(
        archive=archive,
        release=release,
        request_signature=sig,
        query_text=desc,
        payload=payload,
        row_count=0,
        offline_fixture=False,
    )
    return payload, man, False
