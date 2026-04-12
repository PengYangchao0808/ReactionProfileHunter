# mechanism_classifier/AGENTS.md

## OVERVIEW
S0 mechanism classification: builds directed acyclic graphs (DAGs) from SMILES representing reaction mechanisms. Identifies cycloaddition mode ([4+3], [5+2], [4+2], [3+2]), topology (inter/intra), forming/breaking bonds, and precursor/product relationships. Optional — when S0 fails, S2 falls back to SMARTS auto-detection.

## WHERE TO LOOK
| File | Role |
|------|------|
| `classifier.py` | `MechanismClassifier` — main entry; `classify_single()`, `classify_from_csv()`, `export_to_json()`, `export_summary()` |
| `graph_builder.py` | `GraphBuilder` — constructs `MechanismGraph` from reaction records; nodes, edges, pathways |
| `models.py` | Dataclasses: `MechanismGraph`, `GraphNode`, `GraphEdge`, `PathwayInfo`; enums: `CycloMode`, `TopologyType`, `NodeState` |
| `clean_adapter.py` | Adapter for "cleaner" data format (SMILES + bond change annotations) |

## KEY DATA MODEL
```
MechanismGraph (models.py:137)
├── reaction_id, reaction_type, cyclo_mode (CycloMode), topology (TopologyType)
├── nodes: List[GraphNode]  (reactant/intermediate/product/complex)
├── edges: List[GraphEdge]  (forming_bonds, breaking_bonds, ts_type)
├── pathways: List[PathwayInfo]
└── to_json() / from_json() serialization

CycloMode: C4_PLUS_3, C5_PLUS_2, C4_PLUS_2, C3_PLUS_2, UNKNOWN
TopologyType: INTER, INTRA_TYPE_I, INTRA_TYPE_II, INTRA_UNKNOWN
```

## INTEGRATION WITH PIPELINE
- Orchestrator calls `_run_s0()` which uses `MechanismClassifier.classify_single()`
- Output: `S0_Mechanism/mechanism_graph.json` + `mechanism_summary.json`
- S0 summary provides `forming_bonds` to S2 — authoritative source when available
- When S0 fails (degrade policy): S2 uses `SMARTSMatcher` fallback

## CONVENTIONS
- Graph nodes use atom-mapped SMILES for identity tracking
- `forming_bonds` and `breaking_bonds` are stored as `List[Tuple[int, int]]` (0-based atom indices)
- Cycloaddition mode is normalized from free-text (e.g., "[4+3]" → `CycloMode.C4_PLUS_3`)

## ANTI-PATTERNS
- Reimplementing bond identification logic that `GraphBuilder._create_edges()` already handles
- Bypassing S0 when cleaner_data is available — S0 provides more reliable forming bonds than SMARTS
