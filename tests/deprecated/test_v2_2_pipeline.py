import pytest
import os
import sys
from pathlib import Path
import logging

# 将项目根目录添加到 sys.path
sys.path.append(str(Path(__file__).parent.parent))

from rph_core.utils.config_loader import load_config
from rph_core.steps.anchor.handler import AnchorPhase
from rph_core.steps.step2_retro.retro_scanner import RetroScanner
from rph_core.steps.step3_opt.ts_optimizer import TSOptimizer
from rph_core.steps.step4_features.feature_miner import FeatureMiner

def test_v2_2_data_flow():
    """
    验证 Diels-Alder 反应的数据流:
    Ethylene + Butadiene -> Cyclohexene
    """
    config_path = Path("tests/test_fast_config.yaml")
    config = load_config(config_path)
    work_dir = Path("tests/tmp_v2_2_test/da_reaction")
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. AnchorPhase: 处理底物池和产物
    mols = {
        "butadiene": "C=CC=C",
        "ethylene": "C=C",
        "product": "C1=CCCCC1"
    }
    anchor = AnchorPhase(config)
    anchor_res = anchor.run(mols, work_dir / "S1_Anchor")
    
    assert anchor_res.success, f"AnchorPhase 失败: {anchor_res.error_message}"
    assert "butadiene" in anchor_res.anchored_molecules
    assert "ethylene" in anchor_res.anchored_molecules
    assert "product" in anchor_res.anchored_molecules
    
    e_substrate_pool_l2 = (
        anchor_res.anchored_molecules["butadiene"]["e_l2"] + 
        anchor_res.anchored_molecules["ethylene"]["e_l2"]
    )
    print(f"\n[DEBUG] 底物池能量 (L2): {e_substrate_pool_l2:.8f} Hartree")

    # 2. RetroScanner: 生成 TS 初猜
    retro = RetroScanner(config)
    ts_guess, reactant_complex, bonds = retro.run(
        product_xyz=anchor_res.anchored_molecules["product"]["xyz"],
        output_dir=work_dir / "S2_Retro"
    )
    assert ts_guess.exists()
    assert reactant_complex.exists()

    # 3. TSOptimizer: 验证 QCTaskRunner
    # 注意：这里需要 TSOptimizer 支持 run_with_qctaskrunner 并返回正确的 sp_report
    optimizer = TSOptimizer(config)
    ts_analysis_res = optimizer.run_with_qctaskrunner(
        ts_guess=ts_guess,
        reactant=reactant_complex,
        product=anchor_res.anchored_molecules["product"]["xyz"],
        output_dir=work_dir / "S3_TS",
        e_product_l2=anchor_res.anchored_molecules["product"]["e_l2"],
        forming_bonds=None
    )
    
    # 验证 SPMatrixReport
    report = ts_analysis_res.sp_report
    assert report is not None, "TSOptimizer 未返回 SPMatrixReport"
    
    # 注入底物池能量（模拟 Orchestrator 的行为）
    # 在真实流程中，Orchestrator 应该负责校准 report 中的 e_reactant
    report.e_reactant = e_substrate_pool_l2
    print(f"[DEBUG] 校准后的底物能量: {report.e_reactant:.8f} Hartree")

    # 4. FeatureMiner: 特征提取
    miner = FeatureMiner(config)
    features_csv = miner.run(
        ts_final=ts_analysis_res.ts_final_xyz,
        reactant=reactant_complex,
        product=anchor_res.anchored_molecules["product"]["xyz"],
        output_dir=work_dir / "S4_Features",
        forming_bonds=None,
        sp_matrix_report=report
    )
    
    assert features_csv.exists()
    import pandas as pd
    df = pd.read_csv(features_csv)
    print("\n[RESULT] 提取的特征:")
    print(df.to_string())

    # 验证关键物理量
    assert "dG_activation_L2" in df.columns
    assert "E_distortion_total_L2" in df.columns
    
    # 验证是否使用了 L2 能量
    assert df["L2_available"].iloc[0] == 1.0

if __name__ == "__main__":
    test_v2_2_data_flow()
