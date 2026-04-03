# bin/AGENTS.md

## OVERVIEW
Thin CLI wrappers that insert the repo root into `sys.path` and call `rph_core.orchestrator.main()`. Enables running the pipeline without `pip install`.

## WHERE TO LOOK
- `rph_run` and `rph`: identical behavior — add `PROJECT_ROOT` to `sys.path`, then `main()`

## CONVENTIONS
- Keep `bin/` wrappers minimal — all logic lives in `rph_core/`.
- Same behavior as `python -m rph_core` and installed `rph_run` console script.

## ANTI-PATTERNS
- Parsing or modifying pipeline config in bin wrappers — `rph_core/orchestrator.py` owns config resolution.
- Adding step logic or side effects here.
