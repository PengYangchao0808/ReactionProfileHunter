"""
Step 4: Feature Miner
========================

V6.1: 特征提取模块 - 插件管线（extract-only，无QC）

Author: QCcalc Team
Date: 2026-01-09
Updated: 2026-01-28 (Session #V6: 重构为插件管线)
"""

import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Tuple, List

from rph_core.utils.log_manager import LoggerMixin
from rph_core.utils.file_io import read_xyz
from .status import FeatureStatus

logger = logging.getLogger(__name__)


# Constants
HARTREE_TO_KCAL = 627.5094740631  # Hartree to kcal/mol


class FeatureMiner(LoggerMixin):
    """
    V6.1: 特征提取器 (Step 4) - 插件管线版本（extract-only）

    V6.1 变更:
    - 完全移除QC调用（FragmentExtractor, GaussianLogParser等）
    - 集成插件系统（10个extractors）
    - 输出3文件契约：features_raw.csv, features_mlr.csv, feature_meta.json

    输入:
    - ts_final: TS 结构 (来自 Step 3)
    - reactant: 底物 (来自 Step 2)
    - product: 产物 (来自 Step 1)
    - forming_bonds: 形成键的原子索引 ((i, j), (k, l))
    - ts_fchk, ts_orca_out: TS fchk/ORCA输出（用于NBO/interaction_analysis）
    - reactant_fchk, reactant_orca_out: 底物fchk/ORCA输出
    - product_fchk, product_orca_out: 产物fchk/ORCA输出

    输出:
    - features_raw.csv: 完整特征表（所有可用特征）
    - features_mlr.csv: MLR就绪特征表（固定≤10列）
    - feature_meta.json: 元数据（schema, trace, status）
    """

    def __init__(self, config: dict):
        """
        初始化特征提取器

        Args:
            config: 配置字典
        """
        self.config = config
        self.logger.info("FeatureMiner 初始化完成 (V6.1 extract-only plugin pipeline)")

    def run(
        self,
        ts_final: Path,
        reactant: Path,
        product: Path,
        output_dir: Path,
        forming_bonds: Optional[Tuple[Tuple[int, int], Tuple[int, int]]] = None,
        fragment_indices: Optional[Tuple[List[int], List[int]]] = None,
        sp_matrix_report: Optional[object] = None,  # 保留参数兼容性
        ts_fchk: Optional[Path] = None,  # V6.1: fchk input
        ts_orca_out: Optional[Path] = None,  # V6.1: ORCA output
        reactant_fchk: Optional[Path] = None,  # V6.1: reactant fchk
        reactant_orca_out: Optional[Path] = None,  # V6.1: reactant ORCA
        product_fchk: Optional[Path] = None,  # V6.1: product fchk
        product_orca_out: Optional[Path] = None,  # V6.1: product ORCA
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
            ts_fchk: TS fchk 文件（用于NBO/interaction_analysis）
            ts_orca_out: TS ORCA 输出（用于NBO）
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
            DEFAULT_MLR_COLUMNS
        )
        from pathlib import Path
        import time

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        self.logger.info("Step 4 开始: 特征提取 (V6.1 extract-only plugin pipeline)")

        work_dir = output_dir.parent
        s3_dir = (work_dir / "S3_TransitionAnalysis" if (work_dir / "S3_TransitionAnalysis").exists()
                  else work_dir / "S3_TS" if (work_dir / "S3_TS").exists()
                  else None)

        context = FeatureContext(
            ts_xyz=ts_final,
            reactant_xyz=reactant,
            product_xyz=product,
            ts_fchk=ts_fchk,
            reactant_fchk=reactant_fchk,
            product_fchk=product_fchk,
            ts_orca_out=ts_orca_out,
            s3_dir=s3_dir,
            forming_bonds=forming_bonds,
            fragment_indices=fragment_indices,
            sp_report=sp_matrix_report,
            job_run_policy="disallow"
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
        all_warnings = []

        for extractor in extractors:
            plugin_name = extractor.get_plugin_name()
            self.logger.debug(f"Running plugin: {plugin_name}")

            trace = extractor.run(context)

            # Aggregate features
            if trace._extracted_features:
                all_features.update(trace._extracted_features)

            plugin_traces[plugin_name] = trace

            # Collect warnings
            if hasattr(trace, 'warnings') and trace.warnings:
                all_warnings.extend(trace.warnings)

        feature_result = FeatureResult(
            features=all_features,
            schema_version="6.1",
            schema_signature="",
            feature_status=self._aggregate_plugin_status(plugin_traces),
            method="extract-only",
            solvent="",
            temperature_K=298.15,
            enabled_plugins=[e.get_plugin_name() for e in extractors],
            plugin_traces=plugin_traces
        )

        features_raw_csv = output_dir / "features_raw.csv"
        features_mlr_csv = output_dir / "features_mlr.csv"
        feature_meta_json = output_dir / "feature_meta.json"

        self.logger.info("Writing 3-file output contract...")

        write_features_raw_csv(features_raw_csv, all_features)

        write_features_mlr_csv(features_mlr_csv, all_features, DEFAULT_MLR_COLUMNS)

        feature_result.to_json(feature_meta_json)

        if all_warnings:
            self.logger.warning(f"Total warnings from plugins: {len(all_warnings)}")
            for warn in all_warnings[:5]:
                self.logger.warning(f"  - {warn}")

        if feature_result.feature_status == FeatureStatus.OK:
            overall_status = "ok"
        elif feature_result.feature_status == FeatureStatus.PARTIAL:
            overall_status = "partial"
        elif feature_result.feature_status == FeatureStatus.FAILED:
            overall_status = "failed"
        else:
            overall_status = "invalid_inputs"

        self.logger.info(f"✓ 特征提取完成: {features_raw_csv} (status: {overall_status})")
        self.logger.info(f"✓ features_mlr.csv: {features_mlr_csv}")
        self.logger.info(f"✓ feature_meta.json: {feature_meta_json}")

        return features_raw_csv

    def _aggregate_plugin_status(self, traces: dict) -> FeatureStatus:
        if not traces:
            return FeatureStatus.MISSING_INPUTS

        statuses = [trace.status for trace in traces.values()]
        if any(status == FeatureStatus.FAILED for status in statuses):
            return FeatureStatus.FAILED
        if any(status in (FeatureStatus.PARTIAL, FeatureStatus.SKIPPED) for status in statuses):
            return FeatureStatus.PARTIAL
        if all(status == FeatureStatus.OK for status in statuses):
            return FeatureStatus.OK
        return FeatureStatus.INVALID_INPUTS
