# V4 cleanup log

This repository contains only the V4 S0--S4 runtime.  The V3 archive is
maintained outside this worktree and is deliberately not part of the V4
source, test, configuration, or execution surface.

## Removed runtime surfaces

- V3 scheduling, progress display, and legacy checkpoint management.
- `QCTaskRunner`, Berny/QST2/IRC rescue, and the Gau+xTB integration.
- S2 preoptimization, neutral-precursor, and xTB path-search entrypoints.
- Legacy S1 anchor/funnel implementations and V3 S3 `step3_opt`.
- Lewis-acid preoptimization and condition/DR feature aggregation, which are
  outside the V4 DFT-pipeline responsibility.
- Deprecated V3 tests, rescue scripts, root-level experiment runners, and
  Gau+xTB helper sources.

## Retained V4 contract

```text
S0 mechanism -> S1 CENSO-LITE -> S2 PEB -> S3 low-level -> S4 high-level
```

- `irc_valid` remains a nullable result compatibility field; V4 does not
  submit or validate IRC jobs.
- S2 is implemented by `PEBScanner` and its private `PEBScanEngine`.
- Generated QC outputs and the existing V3 output backup are ignored by Git
  and are not deleted by this source cleanup.

## Guard

`tests/test_v4_no_legacy_runtime.py` prevents reintroduction of retired V3
modules and configuration keys into the V4 runtime.
