# `review_kit/` — Isolated Reviewer Environment (Ground Rule R6)

This standalone workspace is reserved strictly for independent candidate review.

## Isolation Guarantees
- **No pipeline code**: `verify_target.py` enforces a runtime import hook blocking `oddstar.channels`, `oddstar.vet`, `oddstar.fuse`, and `oddstar.validate`.
- **No pipeline cache/runs**: Path checks block access to `data/`, `runs/`, `dossiers/`, and `validation/`.
- **Minimal dependencies only**: `numpy`, `scipy`, `astropy`, `astroquery`, `pyvo`, `matplotlib`, and standard library only.
- **Schema compliance**: Emits JSON records conforming strictly to `review/review_record.schema.json`.
