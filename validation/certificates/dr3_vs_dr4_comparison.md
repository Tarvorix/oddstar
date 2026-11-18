# Gaia DR3 vs DR4 Re-Validation Comparison Report

- Config hash: DR3 `fe9c2aa57d3c`, DR4 `d70e74340b82`
- Validation runs on the offline fixture set (synthetic stars and literature analogues), not on live archive data.
- Open case 4582145718700356096 parallax SNR: DR3=119.2 -> DR4=168.1 (1.41x)

| Channel | DR3 Passed | DR4 Passed | DR3 Blind | DR4 Blind | DR3 Inj DP | DR4 Inj DP | DR3 FAR | DR4 FAR |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `ir_excess` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0000 |
| `ir_timedomain` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0000 |
| `optical_dips` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0000 |
| `radio_narrowband` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0000 |
| `underluminous` | True | True | 1.00 | 1.00 | 0.95 | 0.95 | 0.0000 | 0.0000 |
| `generic_outlier` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0067 |
| `spectral_lines` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0000 |
| `vanishing` | True | True | 1.00 | 1.00 | 1.00 | 1.00 | 0.0000 | 0.0000 |
