# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.
