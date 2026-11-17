# oddstar — Reproducible Anomalous Star Search Pipeline

`oddstar` is a modular, offline-reproducible astrophysical search pipeline across eight observational channels spanning infrared excess, mid-infrared time-domain variability, optical dimming, narrowband drifting radio signals, optical under-luminosity, multi-feature outlier analysis, unresolved spectral emission lines, and historical vanishing/appearing sources.

## Architecture

1. **Ingest (`oddstar.io`)**: Stage 1 is the sole network-capable stage. Queries go through rate-limited archive fetchers into a content-addressed raw cache (`data/raw/<archive>/<sha256>`) with an enforced `10.0 GB` budget and SHA-256 manifest logging.
2. **Sample (`oddstar.sample`)**: Builds parent samples with explicit cut accounting and HEALPix-partitioned tabular storage (`nside=16`).
3. **Channels (`oddstar.channels`)**: Eight independent implementations of the unified `Channel` protocol (`sample`, `measure`, `score`, `evidence`).
4. **Validation Gate (`oddstar.validate`)**: Four mandatory tests per channel (blind benchmark recovery, impostor rejection, synthetic injection completeness grid, empirical false-alarm calibration). `oddstar search` is locked out without a valid certificate matching the exact git commit and config SHA-256 hash.
5. **Vetting (`oddstar.vet`)**: Tri-state (`pass`, `fail`, `not_evaluated`) artifact and astrophysical impostor checks, each backed by an on-disk evidence file.
6. **Fusion & Dossier (`oddstar.fuse`, `oddstar.dossier`)**: False-alarm calibration, rule-based Tier `A`/`B`/`C`/`Explained`/`Rejected` assignment, 12-row hypothesis matrix, cross-channel pattern matching, look-elsewhere accounting, and plain-language dossiers.
7. **Adversarial Review & Lessons (`oddstar.review`, `review_kit/`)**: Isolated raw-archive re-verification across three rounds and CI-enforced lesson capture.
8. **CLI & MCP Server (`oddstar.cli`, `oddstar.server`)**: JSON-first CLI with deterministic exit codes and an MCP stdio tool server wrapping every command.

## Quick Start

```bash
.venv/bin/oddstar doctor --json
.venv/bin/pytest -v
.venv/bin/oddstar inspect 4582145718700356096 --json
```
