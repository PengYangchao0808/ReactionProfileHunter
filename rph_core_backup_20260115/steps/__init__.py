"""
ReactionProfileHunter Steps Module - v3.0
===================================

四大核心步骤模块（分子自治架构）
"""

# v3.0: 新的 AnchorPhase 替代 step1_anchor
try:
    from rph_core.steps.anchor.handler import AnchorPhase
    anchor_available = True
except ImportError:
    anchor_available = False
    AnchorPhase = None

from rph_core.steps.step2_retro import RetroScanner
from rph_core.steps.step3_opt import TSOptimizer
from rph_core.steps.step4_features import FeatureMiner

if anchor_available:
    __all__ = [
        "AnchorPhase",
        "RetroScanner",
        "TSOptimizer",
        "FeatureMiner",
    ]
else:
    __all__ = [
        "RetroScanner",
        "TSOptimizer",
        "FeatureMiner",
    ]
