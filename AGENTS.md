# AGENTS.md — Coding & Operational Rules for `oddstar`

`oddstar` is a reproducible, agent-operable astrophysical pipeline that searches public sky archives across eight independent observational channels for stars exhibiting behaviour unexplained by known astrophysics or instrumental artifacts. Every candidate is systematically vetted against known impostors, fused into calibrated false-alarm tiers with explicit open holes, and subjected to isolated multi-round adversarial review from raw archive data before entering an append-only registry. Null results are published alongside completeness grids measured via synthetic signal injection and empirical false-alarm calibration.

## Ground Rules (R1–R10)

| # | Rule | Enforcement |
| --- | --- | --- |
| R1 | **Recover known results blind.** Each channel must find a frozen set of known objects without being told which they are. | Validation harness (`oddstar validate`) hides benchmark IDs inside a field-star sample; `search` refuses without a valid certificate. |
| R2 | **Reject known impostors and instrument glitches.** Each channel must turn away a frozen set of known false positives. | Same harness negative suite; any impostor promoted above Tier C fails the validation gate. |
| R3 | **Measure sensitivity by injection.** Inject realistic synthetic signals into real data and measure recovery; calibrate false alarms separately. | `validate` writes a completeness grid and empirical false-alarm calibration with confidence intervals per channel. |
| R4 | **Freeze thresholds before looking.** Cuts, thresholds and predictions are fixed in a versioned config before a search run. | Config SHA-256 hash is stamped on every output; changes require a `configs/CHANGELOG.md` entry and are CI-checked. |
| R5 | **Never trust your own result.** Every candidate goes to an independent agent with a non-leading prompt, repeatedly. | Three-round review protocol; no candidate is finalized in the registry without stored review records. |
| R6 | **Reviewers re-verify from raw data.** They re-fetch from the archive and re-measure, not read pipeline intermediates. | `review_kit/` runtime isolation barrier blocks importing `oddstar.channels` / `oddstar.vet` or reading `data/` / `runs/`. |
| R7 | **Turn every finding into a lesson.** Each reviewer finding becomes a generic check, test or skill. | CI (`oddstar.ci_checks`) fails if any finding is closed without a linked test/check or written justification. |
| R8 | **Everything is reproducible.** Any candidate can be regenerated with one command from recorded inputs. | Content-addressed raw cache manifests with SHA-256 hashes, pinned lock file, git commit, and deterministic RNG seed on every run. |
| R9 | **Say what is and is not ruled out.** Publish null results with sensitivity and coverage. | Mandatory 12-row hypothesis table and tri-state checks (`pass`, `fail`, `not_evaluated`); `not_evaluated` never counts as `pass`. |
| R10 | **Agents share mistakes.** A second agent can repeat the first one's error. | Round 3 review requires an independent model family where available plus a human checklist over `not_evaluated` rows. The round-3 ledger event records `independent_family`; `oddstar status` reports rounds, families and checklist counts, and CI fails tracked docs that call a review independent while `independent_family` is false. |

## CLI Command Reference

| Command | Purpose |
| --- | --- |
| `oddstar doctor` | Check Python environment, lock file, cache cap, and archive definitions |
| `oddstar sources list` | List all registered public archives, releases, and rate limits |
| `oddstar fetch <archive> --ids/--cone/--adql` | Fetch raw data into the content-addressed cache and write manifest |
| `oddstar sample build <channel>` | Build parent star sample and record counts after each cut |
| `oddstar validate <channel>` | Run all 4 validation tests and issue a frozen certificate |
| `oddstar search <channel>` | Score the sample (blocked unless a matching certificate exists) |
| `oddstar vet <run>` | Run all applicable artifact/impostor checks on flagged targets |
| `oddstar radio followup --from-run <run>` | Fetch and search the Breakthrough Listen recordings of chosen stars (Tier C or better by default, `--ids`, or `--all-covered`), one file at a time; `search radio_narrowband` itself downloads nothing |
| `oddstar fuse <runs...>` | Join channel scores on `source_id`, calibrate FAP, assign tiers and rank |
| `oddstar dossier <source_id>` | Generate `dossier.json`, `dossier.md`, and evidence plots |
| `oddstar inspect <source_id>` | One-star mode: run all channels, vetting, fusion, and dossier for one ID |
| `oddstar inspect <source_id> --register` | Also register the one-star result in the ledger and export its committed evidence |
| `oddstar evidence export <source_id> [--run-id]` | Write the compact committed evidence directory `registry/evidence/<source_id>/<run_id>/` (no ledger write) |
| `oddstar evidence verify [<source_id>]` | Verify ledger evidence references against directories and SHA-256 hashes; count `legacy_without_evidence` |
| `oddstar runs export <run_id>...` | Write the compact committed record of search or fusion runs to `registry/runs/` (sample ids to `registry/samples/`; no ledger write) |
| `oddstar runs verify` | Check every committed run record against its sample files and its own counts |
| `oddstar report build` / `oddstar report check` | Rebuild `docs/candidates.md` and `reviewed_candidate.md` from committed outputs; fail if they differ from a fresh build |
| `oddstar review packet <ids>` | Build blinded reviewer packet mixed with controls and known objects |
| `oddstar review ingest <file>` | Validate and record reviewer JSON, flagging discrepancies |
| `oddstar review ingest --all` | Ingest every stored record under `review/records/` (idempotent; stored records are never rewritten) |
| `oddstar lessons add` | Create a lesson file in `lessons/` and link automated test/check |
| `oddstar lessons register <file>` | Register an existing lesson file in the ledger (closes its finding when `closes_finding`) |
| `oddstar lessons amend <id> --set FIELD=VALUE --reason` | Correct lesson metadata; registered lessons get a `lesson_updated` ledger event |
| `oddstar lessons record-sha` | Set each lesson's `commit_sha` to the commit that introduced its files |
| `oddstar null-report <run>` | Emit sensitivity completeness report and upper limits for a channel run |
| `oddstar status` | Report certified channels, active runs, open review findings per source, review status, and evidence counts |
| `oddstar show <path> --rows 20` | Inspect tabular/JSON/JSONL artifacts safely in bounded pages |
| `oddstar schema [<type>]` | Print or regenerate JSON Schemas in `schemas/` |

## Definition of Done

A change in this repository is complete only when:
1. All offline unit, integration, vetting fixture, and CI governance tests pass (`pytest`).
2. JSON Schemas in `schemas/` and `review/review_record.schema.json` are regenerated cleanly with zero drift.
3. Affected channel documentation in `docs/channels/<name>.md` and `configs/CHANGELOG.md` (if thresholds changed) are updated.
4. No forbidden nomenclature or leaked answer keys exist in code, tests, docs, or outputs.

## Never-Do List

- Never change a threshold, benchmark list or impostor list in the same commit as search results.
- Never mark a check `pass` without an evidence file.
- Never describe an object as a technosignature, megastructure or alien anything in code, logs, commit messages or outputs.
- Never read reviewer packets' answer keys, and never give a reviewer pipeline outputs.
- Never commit raw data.
