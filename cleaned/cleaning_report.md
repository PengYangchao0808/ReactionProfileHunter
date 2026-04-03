# Reaxys Data Cleaning Report
============================================================

## Version Snapshot

- **Git commit**: uncommitted
- **ENABLE_ATOM_MAP**: True
- **ENABLE_CORE_EXTRACTION**: True
- **ENABLE_MECH_CLASSIFIER**: True
- **ENABLE_COND_MATRIX**: False
- **Input file**: all_reactions.json
- **Valid records**: 74
- **Report generated at**: 2026-03-15T14:34:25

## Overview

- **Total records**: 74
- **GOOD**: 25 (33.8%)
- **WARNING**: 47 (63.5%)
- **ERROR**: 2 (2.7%)

## Field Completeness

| Field | Non-null | Completeness |
|-------|----------|--------------|
| precursor_smiles | 74 | 100.0% |
| product_smiles_main | 72 | 97.3% |
| yield | 73 | 98.6% |
| de | 33 | 44.6% |
| ee | 28 | 37.8% |
| temp_celsius | 69 | 93.2% |
| solvent | 0 | 0.0% |
| base | 0 | 0.0% |

## Yield Statistics

- **Records with yield**: 73/74 (98.6%)
- **Average yield**: 60.8%
- **Median yield**: 65.0%
- **Min yield**: 0.0%
- **Max yield**: 91.0%

## Stereochemistry

- **Records with de**: 33 (44.6%)
- **Records with ee**: 28 (37.8%)
- **Records with stereo markers**: 18 (24.3%)
- **Stereo consistency OK**: 16
- **Stereo consistency MISMATCH**: 46
- **Stereo consistency NA (no claim)**: 12

## Mass Balance Validation

- **Passed**: 0/74 (0.0%)
- **Not Evaluable**: 70/74 (94.6%)  (missing oxidant / known reaction pathway)
- **Failed**: 2/74 (2.7%)

## Leaving Group Distribution

- **UNKNOWN**: 2 (2.7%)

## Solvent Distribution

| Solvent | Count | Percentage |
|---------|-------|------------|

## Temperature Statistics

- **Records with temperature**: 69/74
- **Average temperature**: -51.0°C
- **Min temperature**: -78°C
- **Max temperature**: 25°C

## Requires Manual Review

### ERROR Records

| rx_id | Issue |
|-------|-------|
| N/A | SMILES error: EMPTY_PRODUCT; Mass balance failed;  |
| N/A | SMILES error: EMPTY_PRODUCT; Mass balance failed;  |

### WARNING Records

_(Deduplicated: 47 rows → 22 unique reactions)_

| rx_id | Issue |
|-------|-------|
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Unknown leaving group: C5H6O2; Missing yield data; |
| N/A | Stereochemistry claimed but no markers |
| N/A | Unknown leaving group: C9H6O2; Stereochemistry cla |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |
| N/A | Stereochemistry claimed but no markers |

## T4a: Reaction Type Distribution

| reaction_type | count | percentage | avg_confidence |
|---------------|-------|------------|----------------|
| 4+3 | 70 | 94.6% | N/A |
| derivatization | 2 | 2.7% | N/A |
| unknown | 2 | 2.7% | N/A |

## T4b: Atom Mapping Success Rate

| precursor_type | total | success_count | success_rate |
|----------------|-------|---------------|--------------|
| allenamide | 71 | 70 | 98.6% |
| unknown | 3 | 2 | 66.7% |

## T4c: Core Extraction Success Rate

| reaction_type | total | success_count | success_rate |
|---------------|-------|---------------|--------------|
| 4+3 | 70 | 64 | 91.4% |
| derivatization | 2 | 0 | 0.0% |
| unknown | 2 | 0 | 0.0% |

## T4d: Condition Field Parsing Coverage

| field | parsed_count | total | coverage_rate |
|-------|--------------|-------|---------------|
| temp_celsius | 69 | 74 | 93.2% |
| solvent | 0 | 74 | 0.0% |
| base | 0 | 74 | 0.0% |
| time_total_h | 33 | 74 | 44.6% |

## T4e: Yield Distribution Histogram

Bins: 0-20% | 20-40% | 40-60% | 60-80% | 80-100%
Complete             : #### (4) | ###### (6) | ################### (19) | ########################### (27) | ################# (17)
Censored/Upper-bound :  (0) |  (0) |  (0) |  (0) |  (0)

## T4f: Leaving Group Distribution (only [5+2])

| leaving_group | count | percentage |
|---------------|-------|------------|
| N/A | 0 | 0.0% |

## T4g: Solvent/Reagent Vocabulary TOP K + Hash

### Solvent TOP 10

| solvent | count |
|---------|-------|

- **Solvent unique values**: 0
- **Solvent SHA256**: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`

### Base/Reagent TOP 10

| base | count |
|------|-------|

- **Base unique values**: 0
- **Base SHA256**: `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`

## T4h: Version Info + Runtime

- **RDKit version**: 2025.09.3
- **RXNMapper version**: 0.4.3
- **Runtime**: N/A
- **Report generated at**: 2026-03-15T14:34:25

## T4i: Mapping Availability Statistics

- **Usable mappings (OK/LOW_CONFIDENCE)**: 72 (97.3%)
- **Unusable mappings**: 2 (2.7%)
- **Average confidence**: 0.567
- **Confidence distribution**: <0.5: 33, 0.5-0.8: 26, >=0.8: 13

## T4j: Ring Size Distribution by Precursor Type

| precursor_type | 5 | 6 | 7 | 8 | 9 | other/None |
|----------------|---|---|---|---|---|------------|
| allenamide | 0 | 0 | 70 | 0 | 0 | 1 |
| unknown | 0 | 0 | 0 | 0 | 0 | 3 |

## T4k: Precursor Subtype Distribution

| precursor_type | precursor_subtype | count | percentage |
|----------------|-------------------|-------|------------|
| allenamide | general | 71 | 95.9% |
| unknown | unknown | 3 | 4.1% |

## T4l: Core Extraction Fail Distribution

| Reason | Count | Percentage |
|--------|-------|------------|
| AMBIGUOUS_AFTER_EVAL | 8 | 80.0% |
| OTHER | 2 | 20.0% |

## T4m: Mass Balance Oxidation Diagnosis

| Category | Count | Percentage |
|----------|-------|------------|
| ALLENAMIDE_OXIDATION_EXPECTED | 0 | 0.0% |
| OXIDANT_NOT_RECORDED | 70 | 94.6% |
| OXIDANT_POSSIBLE | 0 | 0.0% |
| MISSING_COREACTANT_O | 0 | 0.0% |
| MISSING_COREACTANT_OTHER | 0 | 0.0% |
| PASSED | 0 | 0.0% |

