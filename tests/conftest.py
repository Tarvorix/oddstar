"""Shared test setup: clear process-wide archive state so one test's simulated outage cannot leak into the next."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from oddstar.io import dust3d, live_survey


@pytest.fixture(autouse=True)
def _reset_archive_state() -> Iterator[None]:
    dust3d.reset_unreachable()
    live_survey._vizier_mirror_serves.clear()
    yield
    dust3d.reset_unreachable()
    live_survey._vizier_mirror_serves.clear()
