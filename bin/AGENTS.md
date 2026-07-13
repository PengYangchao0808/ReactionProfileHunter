# bin/AGENTS.md

## OVERVIEW
Thin CLI wrappers that insert the repo root into `sys.path` and call `rph_core.v4_orchestrator.main()`. Enables running the V4 pipeline without `pip install`.

## WHERE TO LOOK
- `rph_run`: primary V4 entry point — calls `rph_core.v4_orchestrator.main()`
- `rph_watch`: WSL live status viewer — calls `rph_core.v4_watch.main()`

## V4 CLI FLAGS
```
bin/rph_run --csv <trusted.csv> --rx-id <id> --output <dir> [--config <yaml>] [--stop-after {s0,s1,s2,s3,s4}] [--log-level <level>]
bin/rph_watch --output <run> [--watch]
```

## CONVENTIONS
- Keep `bin/` wrappers minimal — all logic lives in `rph_core/`.
- Same behavior as `python -m rph_core`.

## ANTI-PATTERNS
- Parsing or modifying pipeline config in bin wrappers — `rph_core/v4_orchestrator.py` owns config resolution.
- Adding step logic or side effects here.
- Referencing the deleted `rph_core.orchestrator` module — V4 uses `v4_orchestrator`.
