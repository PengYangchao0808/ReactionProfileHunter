"""
Step 4: Feature Miner
=======================

V6.2: 特征提取模块 - 插件管线（extract-only，无QC）

Author: QCcalc Team
Date: 2026-01-09
Updated: 2026-02-02 (V6.2: Step1/Step2 机理感知特征)
"""

import logging
import json
from pathlib import Path
from typing import Optional, Tuple, List, Any
from datetime import datetime

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.ui import get_progress_manager
from rph_core.utils.forming_bonds_resolver import resolve_forming_bonds
from rph_core.steps.step4_features.mech_packager import resolve_mechanism_context
from .status import FeatureStatus, aggregate_plugin_status
from .schema import FeatureSchema, DEPLOYABLE_COLUMNS_V1, COMPUTATIONAL_FEATURES, QA_METADATA_FEATURES, validate_feature_layer

logger = logging.getLogger(__name__)


class FeatureMiner(LoggerMixin):
    """
    V6.2: 特征提取器 (Step 4) - 插件管线版本（extract-only）

    V6.2 变更:
    - 新增 Step1 activation 特征提取器 (s1_* 前缀)
    - 新增 Step2 cyclization 特征提取器 (s2_* 前缀)
    - 集成 CDFT、GEDT、虚频校验
    - 12个 extractors (原10个 + 2个新)

    输入:
    - ts_final: TS 结构 (来自 Step 3)
    - reactant: 底物 (来自 Step 2)
    - product: 产物 (来自 Step 1)
    - forming_bonds: 形成键的原子索引 ((i, j), (k, l))
    - ts_fchk, ts_orca_out: TS fchk/ORCA输出（用于TS质量解析/交互分析，纯解析）
    - reactant_fchk, reactant_orca_out: 底物fchk/ORCA输出（可选，纯解析）
    - product_fchk, product_orca_out: 产物fchk/ORCA输出（可选，纯解析）

    输出:
    - features_raw.csv: 完整特征表（所有可用特征）
    - features_mlr.csv: MLR就绪特征表（列由 config.step4.mlr.columns 决定；未配置则使用默认去重列集）
    - feature_meta.json: 元数据（schema, trace, status）
    """

    def __init__(self, config: dict[str, Any]):
        """
        初始化特征提取器

        Args:
            config: 配置字典
        """
        self.config = config
        self.logger.info("FeatureMiner 初始化完成 (V6.2 extract-only plugin pipeline)")

    def run(
        self,
        ts_final: Optional[Path],
        reactant: Optional[Path],
        product: Optional[Path],
        output_dir: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], ...]] = None,
        fragment_indices: Optional[Tuple[List[int], List[int]]] = None,
        sp_matrix_report: Optional[object] = None,  # 保留参数兼容性
        ts_fchk: Optional[Path] = None,  # V6.1: fchk input
        ts_orca_out: Optional[Path] = None,  # V6.1: ORCA output
        reactant_fchk: Optional[Path] = None,  # V6.1: reactant fchk
        reactant_orca_out: Optional[Path] = None,  # V6.1: reactant ORCA
        product_fchk: Optional[Path] = None,  # V6.1: product fchk
        product_orca_out: Optional[Path] = None,  # V6.1: product ORCA
        # V6.2: Log file inputs for TS quality and frequency analysis
        ts_log: Optional[Path] = None,
        reactant_log: Optional[Path] = None,
        product_log: Optional[Path] = None,
        # V6.2: S1 inputs for Step1 activation features
        s1_dir: Optional[Path] = None,
        s1_shermo_summary_file: Optional[Path] = None,
        s1_hoac_thermo_file: Optional[Path] = None,
        s1_precursor_xyz: Optional[Path] = None,
        s1_conformer_energies_file: Optional[Path] = None,
    ) -> Path:
        """
        V6.1: 执行特征提取（extract-only插件管线）

        Args:
            ts_final: TS 结构 (XYZ 或 .log)
            reactant: 底物结构 (XYZ 或 .log)
            product: 产物结构 (XYZ 或 .log)
            output_dir: 输出目录
            forming_bonds: 形成键的原子索引 ((i, j), (k, l))
            fragment_indices: 双片段索引 (fragment_A_atoms, fragment_B_atoms)
            sp_matrix_report: S3.5 SP 矩阵报告（保留兼容性）
            ts_fchk: TS fchk 文件（用于TS质量解析/交互分析）
            ts_orca_out: TS ORCA 输出（用于TS质量解析）
            reactant_fchk: 底物fchk 文件
            reactant_orca_out: 底物ORCA 输出
            product_fchk: 产物fchk 文件
            product_orca_out: 产物ORCA 输出

        Returns:
            features_raw.csv: 特征表路径（主输出文件）
        """
        from .context import FeatureContext, FeatureResult
        from .extractors import list_extractors
        from .schema import (
            write_features_raw_csv,
            write_features_mlr_csv,
            get_schema_signature
        )

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Step 4 开始: 特征提取 (V6.1 extract-only plugin pipeline)")
        pm = get_progress_manager()

        work_dir = output_dir.parent
        s3_dir = (work_dir / "S3_TransitionAnalysis" if (work_dir / "S3_TransitionAnalysis").exists()
                  else work_dir / "S3_TS" if (work_dir / "S3_TS").exists()
                  else None)

        artifacts_index = None
        if s3_dir is not None:
            index_path = s3_dir / "artifacts_index.json"
            if index_path.exists():
                try:
                    artifacts_index = json.loads(index_path.read_text())
                except Exception as e:
                    self.logger.warning(f"Failed to read artifacts_index.json: {e}")

        if ts_final is None or reactant is None or product is None:
            mech_cfg = self.config.get('step4', {}).get('mechanism_packaging', {})
            mech_context = resolve_mechanism_context(output_dir, mech_cfg)
            if ts_final is None and mech_context.s3_ts_final:
                ts_final = mech_context.s3_ts_final
            if product is None and mech_context.s1_product:
                product = mech_context.s1_product
            if reactant is None:
                reactant = mech_context.s3_reactant_sp or mech_context.s2_intermediate

        if forming_bonds is None:
            forming_cfg = self.config.get('step4', {}).get('forming_bonds', {}) or {}
            resolved = resolve_forming_bonds(
                product_xyz=product,
                ts_xyz=ts_final,
                s3_dir=s3_dir,
                s4_dir=output_dir,
                config=forming_cfg,
                write_meta=forming_cfg.get('write_meta', True)
            )
            forming_bonds = resolved.forming_bonds
            if resolved.warnings:
                for w in resolved.warnings:
                    self.logger.warning(f"Forming bonds resolver: {w}")

        context = FeatureContext(
            ts_xyz=ts_final,
            reactant_xyz=reactant,
            product_xyz=product,
            ts_fchk=ts_fchk,
            reactant_fchk=reactant_fchk,
            product_fchk=product_fchk,
            ts_orca_out=ts_orca_out,
            s3_dir=s3_dir,
            artifacts_index=artifacts_index,
            forming_bonds=forming_bonds,
            fragment_indices=fragment_indices,
            sp_report=sp_matrix_report,
            job_run_policy="disallow",
            # V6.2 P0: Inject S3 path handles with s3_ prefix for extractors
            s3_ts_fchk=ts_fchk,
            s3_reactant_fchk=reactant_fchk,
            s3_ts_log=ts_log,
            s3_reactant_g_kcal=(getattr(sp_matrix_report, "g_reactant", None) if sp_matrix_report else None),
            # V6.2 P0: Inject log file handles
            ts_log=ts_log,
            reactant_log=reactant_log,
            product_log=product_log,
            # V6.2 P0: Inject S1 path handles for Step1 activation features
            s1_dir=s1_dir,
            s1_shermo_summary_file=s1_shermo_summary_file,
            s1_hoac_thermo_file=s1_hoac_thermo_file,
            s1_precursor_xyz=s1_precursor_xyz,
            s1_conformer_energies_file=s1_conformer_energies_file,
        )

        enabled_plugins = self.config.get('step4', {}).get('enabled_plugins', None)
        if enabled_plugins is None:
            self.logger.info("No enabled_plugins specified, using all registered extractors")
            extractors = list_extractors()
        else:
            extractors = list_extractors()
            self.logger.info(f"Enabled plugins: {enabled_plugins}")
            extractors = [e for e in extractors if e.get_plugin_name() in enabled_plugins]

        all_features = {}
        plugin_traces = {}
        all_warnings = []  # List of warning dicts: {'code': str, 'plugin': str, 'severity': str, 'detail': str}

        for extractor in extractors:
            plugin_name = extractor.get_plugin_name()
            self.logger.debug(f"Running plugin: {plugin_name}")
            pm.update_step("s4", description=f"S4: Running plugin '{plugin_name}'...")

            trace = extractor.run(context)

            # Aggregate features
            if trace._extracted_features:
                all_features.update(trace._extracted_features)

            plugin_traces[plugin_name] = trace

            # Collect warnings from plugin trace
            if hasattr(trace, 'warnings') and trace.warnings:
                for w in trace.warnings:
                    # Handle both string warnings and dict warnings
                    if isinstance(w, dict):
                        all_warnings.append(w)
                    else:
                        # Convert string warning code to structured format
                        all_warnings.append({
                            'code': w,
                            'plugin': plugin_name,
                            'severity': 'warn',
                            'detail': str(w)
                        })

        # V6.2: Config snapshot and provenance
        config_snapshot = self.config.get('step4', {}).get('step4_features', {})
        
        # Build provenance
        provenance = {
            "extract_mode": "extract-only",
            "plugin_pipeline_version": "6.2",
            "multiwfn_status": "unknown"
        }
        
        # Check Multiwfn status
        if 'multiwfn_features' in plugin_traces:
            mw_trace = plugin_traces['multiwfn_features']
            mw_features = mw_trace._extracted_features
            provenance["multiwfn_status"] = mw_features.get('mw_status', 'not_run')
            provenance["multiwfn_cache_hit"] = mw_features.get('mw_cache_hit', False)

        aggregated_status = self._aggregate_plugin_status(plugin_traces)
        enabled_plugin_names = [e.get_plugin_name() for e in extractors]
        schema_sig = get_schema_signature(features=all_features, enabled_plugins=enabled_plugin_names)
        
        sample_id = output_dir.parent.name
        
        all_features["schema_version"] = "6.2"
        all_features["schema_signature"] = schema_sig
        all_features["feature_status"] = aggregated_status.value
        all_features["sample_id"] = sample_id
        all_features["qc.warnings_count"] = len(all_warnings)

        schema_warnings = FeatureSchema().validate_row(all_features)
        for warning_msg in schema_warnings:
            all_warnings.append({
                'code': 'W_SCHEMA_VALIDATION',
                'plugin': 'schema',
                'severity': 'warn',
                'detail': warning_msg,
            })

        # Deduplicate warnings by (code, plugin)
        seen_warnings = set()
        deduped_warnings = []
        for w in all_warnings:
            key = (w.get('code'), w.get('plugin'))
            if key not in seen_warnings:
                seen_warnings.add(key)
                deduped_warnings.append(w)
        deduped_warnings_list = deduped_warnings  # For FeatureResult.warnings field
        all_features["qc.warnings_count"] = len(deduped_warnings_list)

        # Set qc.sample_weight based on policy
        # Policy: 1.0 only when feature_status == OK AND TS validity is ok
        #        0.0 for INVALID_INPUTS / FAILED / TS invalid
        #        0.5 for PARTIAL (optional)
        ts_valid = True  # Default to valid
        if aggregated_status == FeatureStatus.OK:
            sample_weight = 1.0
        elif aggregated_status in [FeatureStatus.INVALID_INPUTS, FeatureStatus.FAILED]:
            sample_weight = 0.0
        else:  # PARTIAL or other
            sample_weight = 0.0  # Degraded samples get 0 weight
        
        all_features["qc.sample_weight"] = sample_weight

        # Build artifact_presence tracking
        artifact_presence = {
            "ts_xyz": ts_final is not None and ts_final.exists(),
            "ts_fchk": ts_fchk is not None and ts_fchk.exists(),
            "ts_log": ts_log is not None and ts_log.exists(),
            "reactant_fchk": reactant_fchk is not None and reactant_fchk.exists(),
            "product_fchk": product_fchk is not None and product_fchk.exists(),
            "s1_shermo_summary_file": s1_shermo_summary_file is not None and s1_shermo_summary_file.exists(),
            "s1_hoac_thermo_file": s1_hoac_thermo_file is not None and s1_hoac_thermo_file.exists(),
        }

        feature_result = FeatureResult(
            features=all_features,
            schema_version="6.2",
            schema_signature=schema_sig,
            feature_status=aggregated_status,
            method="extract-only",
            solvent="",
            temperature_K=config_snapshot.get('temperature_K', 298.15),
            enabled_plugins=enabled_plugin_names,
            warnings=deduped_warnings_list,
            plugin_traces=plugin_traces,
            config_snapshot=config_snapshot,
            provenance=provenance,
            artifact_presence=artifact_presence
        )

        features_raw_csv = output_dir / "features_raw.csv"
        features_mlr_csv = output_dir / "features_mlr.csv"
        feature_meta_json = output_dir / "feature_meta.json"

        self.logger.info("Writing 3-file output contract...")

        write_features_raw_csv(features_raw_csv, all_features)

        step4_cfg = self.config.get('step4', {}) or {}
        mlr_cfg = step4_cfg.get('mlr', {}) or {}
        mlr_columns = mlr_cfg.get('columns')
        if mlr_columns is not None and not isinstance(mlr_columns, list):
            self.logger.warning(
                f"step4.mlr.columns should be a list, got {type(mlr_columns).__name__}; using default"
            )
            mlr_columns = None

        write_features_mlr_csv(features_mlr_csv, all_features, mlr_columns)

        self.logger.info("Writing V6.4 three-layer output...")

        deployable_csv = output_dir / "deployable_features.csv"
        qa_metadata_csv = output_dir / "qa_metadata.csv"

        import pandas as pd

        df_full = pd.DataFrame([all_features])

        deployable_cols = [c for c in DEPLOYABLE_COLUMNS_V1 if c in df_full.columns]
        df_deployable = df_full[deployable_cols] if deployable_cols else pd.DataFrame()
        if not df_deployable.empty:
            df_deployable.to_csv(deployable_csv, index=False)
            self.logger.info(f"✓ deployable_features.csv: {len(deployable_cols)} columns (Layer 1 + Layer 2)")

        qa_cols = [c for c in QA_METADATA_FEATURES if c in df_full.columns]
        df_qa = df_full[qa_cols] if qa_cols else pd.DataFrame()
        if not df_qa.empty:
            df_qa.to_csv(qa_metadata_csv, index=False)
            self.logger.info(f"✓ qa_metadata.csv: {len(qa_cols)} columns (Layer 3)")

        feature_result.to_json(feature_meta_json)

        if all_warnings:
            self.logger.warning(f"Total warnings from plugins: {len(all_warnings)}")
            warn_by_plugin = {}
            for w in all_warnings:
                p = w.get('plugin', 'unknown')
                warn_by_plugin.setdefault(p, []).append(w.get('detail', w.get('code')))
            
            for p, details in warn_by_plugin.items():
                self.logger.warning(f"  Plugin '{p}': {len(details)} warnings")
                for d in details[:3]:
                    self.logger.warning(f"    - {d}")
                if len(details) > 3:
                    self.logger.warning(f"    - ... and {len(details)-3} more")

        if feature_result.feature_status == FeatureStatus.OK:
            overall_status = "ok"
        elif feature_result.feature_status == FeatureStatus.PARTIAL:
            overall_status = "partial"
        elif feature_result.feature_status == FeatureStatus.FAILED:
            overall_status = "failed"
        else:
            overall_status = "invalid_inputs"

        try:
            status_file = output_dir / ".rph_step_status.json"
            step_status = {
                "step": "s4",
                "status": overall_status,
                "timestamp": datetime.now().isoformat(),
                "plugins_count": len(extractors),
                "warnings_count": len(all_warnings)
            }
            status_file.write_text(json.dumps(step_status, indent=2))
        except Exception as e:
            self.logger.debug(f"Failed to write .rph_step_status.json: {e}")

        self.logger.info(f"✓ 特征提取完成: {features_raw_csv} (status: {overall_status})")
        self.logger.info(f"✓ features_mlr.csv: {features_mlr_csv}")
        self.logger.info(f"✓ feature_meta.json: {feature_meta_json}")

        return features_raw_csv

    def _aggregate_plugin_status(self, traces: dict[str, Any]) -> FeatureStatus:
        if not traces:
            return FeatureStatus.MISSING_INPUTS
            
        plugin_statuses = {name: trace.status for name, trace in traces.items()}
        return aggregate_plugin_status(plugin_statuses)
