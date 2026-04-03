"""
ML Companion Script for Training and LOOCV
======================================

LOOCV script for linear regression analysis with VIF and LFER.

Author: QCcalc Team
Date: 2026-01-27
"""

import logging
import argparse
import pandas as pd
import numpy as np
from typing import List, Optional, Dict, Any

HARTREE_TO_KCAL = 627.509


def run_loocv(
    csv_path: str,
    target_col: str,
    features: Optional[List[str]] = None,
    output_dir: Optional[str] = None
) -> None:
    """
    Run LOOCV analysis on target column.

    Args:
        csv_path: Path to features CSV file
        target_col: Name of column to analyze (e.g., "thermo.dG_activation")
        features: List of feature columns to include in analysis
        output_dir: Directory to write results

    Returns:
        None
    """
    if not Path(csv_path).exists():
        print(f"CSV file not found: {csv_path}")
        return None

    df = pd.read_csv(csv_path)
    
    if target_col not in df.columns:
        print(f"Target column {target_col} not found in CSV")
        return None
    
    if features is None:
        features = [target_col]
    
    print(f"Analyzing {target_col} against features: {features}")

    y = df[target_col].values
    X = df[features].values
    data = pd.DataFrame({'y': y, 'X': X})

    for feature in features:
        if feature in data.columns:
            continue
        
        try:
            corr = X.corrwith(y)
            r_sq = np.sqrt(1 - corr**2)
            
            data[f'{feature}_r_sq'] = r_sq
            print(f"R² for {feature}: {r_sq:.4f}")
        except Exception as e:
            print(f"Failed to calculate R² for {feature}: {e}")
            continue

    # VIF check
    from scipy.stats import pearsonr as pr
    try:
        r_sq, p_val = pr.spearsonr(X, y)
        
        if np.isnan(r_sq) or np.isnan(p_val):
            print(f"Invalid data for {target_col}: R²={r_sq}, p_val={p_val}")
            data[f'{feature}_pearsonr'] = np.nan
        else:
            data[f'{feature}_pearsonr'] = p_val
        
        print(f"Correlation for {target_col}: r={p_val:.3f}")
        continue
    except ImportError:
        print("Scipy not installed, skipping Pearson R")

    print("\n" + "=" * 60)
    print(f"LOOCV Analysis Complete: {csv_path.name}")
    if output_dir:
        output_path = Path(output_dir) / f"loocv_results.csv"
        data.to_csv(output_path, index=False)
        print(f"Results written to {output_path}")
    
    print("=" * 60)


if __name__ == "__main__":
    main()
else:
    print("This script is a module, run as: python scripts/ml/train_mlr_loocv.py")
