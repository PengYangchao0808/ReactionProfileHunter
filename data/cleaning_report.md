# Reaxys Data Cleaning Report
============================================================

## Version Snapshot

- **Git commit**: uncommitted
- **ENABLE_ATOM_MAP**: True
- **ENABLE_CORE_EXTRACTION**: True
- **ENABLE_MECH_CLASSIFIER**: True
- **ENABLE_COND_MATRIX**: False
- **Input file**: xiong_2003_for_cleaner.xlsx
- **Valid records**: 20
- **Report generated at**: 2026-04-27T23:10:53

## Overview

- **Total records**: 20
- **GOOD**: 20 (100.0%)
- **WARNING**: 0 (0.0%)
- **ERROR**: 0 (0.0%)

## Field Completeness

| Field | Non-null | Completeness |
|-------|----------|--------------|
| precursor_smiles | 20 | 100.0% |
| product_smiles_main | 20 | 100.0% |
| yield | 20 | 100.0% |
| de | 16 | 80.0% |
| ee | 0 | 0.0% |
| temp_celsius | 20 | 100.0% |
| solvent | 20 | 100.0% |
| base | 20 | 100.0% |

## Yield Statistics

- **Records with yield**: 20/20 (100.0%)
- **Average yield**: 58.0%
- **Median yield**: 72.5%
- **Min yield**: 5.0%
- **Max yield**: 90.0%

## Stereochemistry

- **Records with de**: 16 (80.0%)
- **Records with ee**: 0 (0.0%)
- **Records with stereo markers**: 20 (100.0%)
- **Stereo consistency OK**: 16
- **Stereo consistency MISMATCH**: 0
- **Stereo consistency NA (no claim)**: 4

## Mass Balance Validation

- **Passed**: 0/20 (0.0%)
- **Not Evaluable**: 20/20 (100.0%)  (missing oxidant / known reaction pathway)
- **Failed**: 0/20 (0.0%)

## Leaving Group Distribution


## Solvent Distribution

| Solvent | Count | Percentage |
|---------|-------|------------|
| CH2Cl2 | 20 | 100.0% |

## Temperature Statistics

- **Records with temperature**: 20/20
- **Average temperature**: -25.6°C
- **Min temperature**: -78°C
- **Max temperature**: 25°C

## T4a: Reaction Type Distribution

| reaction_type | count | percentage | avg_confidence |
|---------------|-------|------------|----------------|
| 4+3 | 20 | 100.0% | N/A |

## T4b: Mapping Trust Rate (chemically validated)

| precursor_type | total | accepted | trusted_mechanistic |
|----------------|-------|----------|---------------------|
| allenamide | 20 | 20 | 20 |

## T4c: Core Extraction Status

| reaction_type | total | OK | AMBIGUOUS | ERROR |
|---------------|-------|----|-----------|-------|
| 4+3 | 20 | 20 | 0 | 0 |

## T4d: Condition Field Parsing Coverage

| field | parsed_count | total | coverage_rate |
|-------|--------------|-------|---------------|
| temp_celsius | 20 | 20 | 100.0% |
| solvent | 20 | 20 | 100.0% |
| base | 20 | 20 | 100.0% |
| time_total_h | 0 | 20 | 0.0% |

## T4e: Yield Distribution Histogram

Bins: 0-20% | 20-40% | 40-60% | 60-80% | 80-100%
Complete             :  (0) | ### (3) | ### (3) | ######## (8) | #### (4)
Censored/Upper-bound : ## (2) |  (0) |  (0) |  (0) |  (0)

## T4f: Leaving Group Distribution (only [5+2])

| leaving_group | count | percentage |
|---------------|-------|------------|
| N/A | 0 | 0.0% |

## T4g: Solvent/Reagent Vocabulary TOP K + Hash

### Solvent TOP 10

| solvent | count |
|---------|-------|
| CH2Cl2 | 20 |

- **Solvent unique values**: 1
- **Solvent SHA256**: `fd363c86931063b6f0409b88df70a72bc3149ff2c985f60070020dfbf55c7296`

### Base/Reagent TOP 10

| base | count |
|------|-------|
| DMDO | 20 |

- **Base unique values**: 1
- **Base SHA256**: `73a09822a5776a63fa14ea085a268850a733c333801c3ee67dada554c006c901`

## T4h: Version Info + Runtime

- **RDKit version**: 2025.09.3
- **RXNMapper version**: 0.4.3
- **Runtime**: N/A
- **Report generated at**: 2026-04-27T23:10:53

## T4i: Mapping Availability Statistics

- **Usable mappings (OK/LOW_CONFIDENCE)**: 20 (100.0%)
- **Unusable mappings**: 0 (0.0%)
- **Average confidence**: 0.824
- **Confidence distribution**: <0.5: 0, 0.5-0.8: 5, >=0.8: 15

## T4j: Ring Size Distribution by Precursor Type

| precursor_type | 5 | 6 | 7 | 8 | 9 | other/None |
|----------------|---|---|---|---|---|------------|
| allenamide | 0 | 0 | 17 | 0 | 0 | 3 |

## T4k: Precursor Subtype Distribution

| precursor_type | precursor_subtype | count | percentage |
|----------------|-------------------|-------|------------|
| allenamide | general | 20 | 100.0% |

## T4l: Core Extraction Fail Distribution

No extraction failures

## T4m: Atom Mapping Trust Summary

- **Mapping generated**: 20 (total records)
- **Mapping accepted**: 20 (100.0%)
- **Mapping needs review**: 0 (0.0%)
- **Mapping failed**: 0 (0.0%)
- **Identity validation PASS**: 20
- **Identity validation REVIEW**: 0
- **Identity validation VETO**: 0
- **Trusted for graph delta**: 20 (100.0%)
- **Trusted for mechanistic interpretation**: 20 (100.0%)

## T4n: Resolver Status Distribution

| resolver_status | count | percentage |
|-----------------|-------|------------|
| RESOLVED | 20 | 100.0% |

## T4o: Resolver Name Distribution

| resolver_name | count | percentage |
|---------------|-------|------------|
| dearomatization_cycloaddition | 20 | 100.0% |

## T4m: Mass Balance Oxidation Diagnosis

| Category | Count | Percentage |
|----------|-------|------------|
| ALLENAMIDE_OXIDATION_EXPECTED | 20 | 100.0% |
| OXIDANT_NOT_RECORDED | 0 | 0.0% |
| OXIDANT_POSSIBLE | 0 | 0.0% |
| MISSING_COREACTANT_O | 0 | 0.0% |
| MISSING_COREACTANT_OTHER | 0 | 0.0% |
| PASSED | 0 | 0.0% |

