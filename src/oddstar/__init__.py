"""oddstar: Reproducible multi-channel anomalous star search pipeline."""

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/oddstar-mplconfig")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

__version__ = "0.1.0"
