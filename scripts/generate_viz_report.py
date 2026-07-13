#!/usr/bin/env python3
"""
RPH V1-V14 Comprehensive Visualization Report
Generates publication-quality figures for the [4+3] cycloaddition ML project.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.gridspec import GridSpec
import os

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
OUTPUT_DIR = "/mnt/e/Calculations/AI4S_ML_Studys/[4+3] Mechain learning/ReactionProfileHunter/RPH_V2.1.1/docs/figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Color palette
C_DR = '#2E86AB'      # Deep blue for DR
C_YIELD = '#A23B72'   # Magenta for Yield
C_POS = '#28A745'     # Green for positive
C_NEG = '#DC3545'     # Red for negative
C_NEU = '#6C757D'     # Gray for neutral
C_GOLD = '#F4A261'    # Gold for highlights
C_PURPLE = '#9B5DE5'  # Purple for special
C_BG = '#FAFAFA'      # Background

plt.rcParams.update({
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.edgecolor': '#333333',
    'axes.labelcolor': '#333333',
    'text.color': '#333333',
    'xtick.color': '#555555',
    'ytick.color': '#555555',
    'grid.color': '#E0E0E0',
    'grid.linewidth': 0.5,
    'font.size': 10,
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'legend.fontsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.15,
})

# ---------------------------------------------------------------------------
# Figure 1: Model R² Evolution Across V1-V14
# ---------------------------------------------------------------------------
def fig1_evolution():
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.set_facecolor(C_BG)

    versions = ['V6', 'V8', 'V9', 'V10', 'V11', 'V12', 'V13\n(DR)', 'V13\n(Yield)', 'V14']
    x = np.arange(len(versions))

    # DR R² values (best model each round)
    dr_r2 = [0.177, 0.431, 0.276, 0.276, 0.657, 0.700, 0.700, 0.897, 0.700]
    dr_labels = ['Lasso', 'Ridge', 'GBR', 'GBR', 'Stability-2\n★', 'M1', 'M1', 'G6\nPMPF-4', 'M1\n(frozen)']

    # Yield R² values
    yield_r2 = [0.177, 0.431, 0.276, 0.276, -0.982, 0.550, None, 0.550, 0.449]
    yield_labels = ['Lasso', 'Ridge', 'GBR', 'GBR', 'Core-5b\n(fail)', 'Y-F', None, 'Y-F', 'V14\nNested']

    bar_w = 0.35
    bars1 = ax.bar(x - bar_w/2, dr_r2, bar_w, label='DR Selectivity', color=C_DR, edgecolor='white', linewidth=0.5, zorder=3)
    bars2 = ax.bar(x + bar_w/2, [y if y is not None else np.nan for y in yield_r2], bar_w, label='Yield', color=C_YIELD, edgecolor='white', linewidth=0.5, zorder=3)

    # Add value labels on bars
    for bar, val, lbl in zip(bars1, dr_r2, dr_labels):
        h = bar.get_height()
        ypos = h + 0.02 if h >= 0 else h - 0.08
        ax.text(bar.get_x() + bar.get_width()/2., ypos, f'{val:.3f}', ha='center', va='bottom' if h>=0 else 'top', fontsize=7.5, fontweight='bold', color=C_DR)
        ax.text(bar.get_x() + bar.get_width()/2., ypos + (0.06 if h>=0 else -0.06), lbl, ha='center', va='bottom' if h>=0 else 'top', fontsize=6, color='#666666')

    for bar, val, lbl in zip(bars2, yield_r2, yield_labels):
        if val is None:
            continue
        h = bar.get_height()
        ypos = h + 0.02 if h >= 0 else h - 0.08
        ax.text(bar.get_x() + bar.get_width()/2., ypos, f'{val:.3f}', ha='center', va='bottom' if h>=0 else 'top', fontsize=7.5, fontweight='bold', color=C_YIELD)
        ax.text(bar.get_x() + bar.get_width()/2., ypos + (0.06 if h>=0 else -0.06), lbl, ha='center', va='bottom' if h>=0 else 'top', fontsize=6, color='#666666')

    # Annotations for key events
    ax.annotate('Stability-2\nDiscovery', xy=(4, 0.657), xytext=(3.2, 0.85),
                arrowprops=dict(arrowstyle='->', color=C_GOLD, lw=1.5),
                fontsize=8, color=C_GOLD, fontweight='bold', ha='center')
    ax.annotate('Dual-mechanism\nPicture', xy=(5, 0.700), xytext=(5.5, 0.90),
                arrowprops=dict(arrowstyle='->', color=C_DR, lw=1.5),
                fontsize=8, color=C_DR, fontweight='bold', ha='center')
    ax.annotate('GEDT-dominant\nYield Model', xy=(5, 0.550), xytext=(6.2, 0.35),
                arrowprops=dict(arrowstyle='->', color=C_YIELD, lw=1.5),
                fontsize=8, color=C_YIELD, fontweight='bold', ha='center')

    ax.axhline(y=0, color='#999999', linewidth=0.8, linestyle='-', zorder=1)
    ax.set_ylabel('LORO R²', fontweight='bold')
    ax.set_xlabel('Development Round', fontweight='bold')
    ax.set_title('Fig 1. Model Performance Evolution Across V1–V14\n[4+3] Cycloaddition DR & Yield Prediction', fontweight='bold', pad=15)
    ax.set_xticks(x)
    ax.set_xticklabels(versions, fontsize=9)
    ax.set_ylim(-1.1, 1.0)
    ax.legend(loc='lower left', framealpha=0.9, edgecolor='#CCCCCC')
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig1_evolution.png', facecolor='white')
    plt.close()
    print("✓ fig1_evolution.png")

# ---------------------------------------------------------------------------
# Figure 2: DR vs Yield — Current Best Models Comparison
# ---------------------------------------------------------------------------
def fig2_dr_vs_yield():
    fig = plt.figure(figsize=(12, 5))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1, 1.2, 1])

    # --- Left: R² & MAE comparison ---
    ax1 = fig.add_subplot(gs[0])
    ax1.set_facecolor(C_BG)
    metrics = ['LORO R²', 'MAE\n(kcal/mol)', 'Spearman ρ', 'Features\n(count)']
    dr_vals = [0.700, 0.254, 0.897, 2]
    yield_vals = [0.550, 0.135, 0.843, 8]  # MAE in yield fraction

    x = np.arange(len(metrics))
    bar_w = 0.35
    ax1.bar(x - bar_w/2, dr_vals, bar_w, label='DR: M1 Stability-2', color=C_DR, edgecolor='white', zorder=3)
    ax1.bar(x + bar_w/2, yield_vals, bar_w, label='Yield: Y-F GBR', color=C_YIELD, edgecolor='white', zorder=3)

    for i, (d, y) in enumerate(zip(dr_vals, yield_vals)):
        ax1.text(i - bar_w/2, d + 0.03, f'{d:.3f}' if d < 1 else f'{int(d)}', ha='center', fontsize=8, fontweight='bold', color=C_DR)
        ax1.text(i + bar_w/2, y + 0.03, f'{y:.3f}' if y < 1 else f'{int(y)}', ha='center', fontsize=8, fontweight='bold', color=C_YIELD)

    ax1.set_xticks(x)
    ax1.set_xticklabels(metrics, fontsize=9)
    ax1.set_ylim(0, 1.15)
    ax1.set_title('Performance Metrics', fontweight='bold')
    ax1.legend(loc='upper right', framealpha=0.9, edgecolor='#CCCCCC', fontsize=8)
    ax1.grid(axis='y', zorder=0)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # --- Middle: Physical driver comparison ---
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor(C_BG)
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 10)
    ax2.axis('off')

    # DR box
    rect_dr = mpatches.FancyBboxPatch((0.3, 5.2), 4.2, 4.2, boxstyle="round,pad=0.1", facecolor=C_DR, alpha=0.12, edgecolor=C_DR, linewidth=2)
    ax2.add_patch(rect_dr)
    ax2.text(2.4, 8.8, 'DR Selectivity', fontsize=11, fontweight='bold', color=C_DR, ha='center')
    ax2.text(2.4, 8.1, 'M1 Stability-2', fontsize=9, color=C_DR, ha='center')
    ax2.text(2.4, 7.4, 'R² = 0.700', fontsize=9, fontweight='bold', color=C_DR, ha='center')
    ax2.text(2.4, 6.6, 'Driver: S1 Conformational\nPreorganization', fontsize=8, color='#444444', ha='center')
    ax2.text(2.4, 5.7, '• close_contacts_density\n• post_bond_length_q10', fontsize=7.5, color='#555555', ha='center', family='monospace')

    # Yield box
    rect_y = mpatches.FancyBboxPatch((5.5, 5.2), 4.2, 4.2, boxstyle="round,pad=0.1", facecolor=C_YIELD, alpha=0.12, edgecolor=C_YIELD, linewidth=2)
    ax2.add_patch(rect_y)
    ax2.text(7.6, 8.8, 'Yield', fontsize=11, fontweight='bold', color=C_YIELD, ha='center')
    ax2.text(7.6, 8.1, 'Y-F GBR (8-feat)', fontsize=9, color=C_YIELD, ha='center')
    ax2.text(7.6, 7.4, 'R² = 0.550', fontsize=9, fontweight='bold', color=C_YIELD, ha='center')
    ax2.text(7.6, 6.6, 'Driver: TS Electronic\nTransfer (GEDT)', fontsize=8, color='#444444', ha='center')
    ax2.text(7.6, 5.7, '• GEDT (dominant)\n• ccd + post_bond\n• conformational', fontsize=7.5, color='#555555', ha='center', family='monospace')

    # Independence arrow
    ax2.annotate('', xy=(5.3, 7.3), xytext=(4.7, 7.3),
                 arrowprops=dict(arrowstyle='<->', color=C_NEU, lw=2))
    ax2.text(5.0, 7.6, 'ρ ≈ 0\nIndependent', fontsize=8, color=C_NEU, ha='center', fontweight='bold')

    # Bottom: dual mechanism summary
    ax2.text(5.0, 4.0, 'Dual-Mechanism Picture', fontsize=12, fontweight='bold', color='#333333', ha='center')
    ax2.text(5.0, 3.2, 'DR → What shape is the TS? (preorganization)\nYield → Can the TS form efficiently? (charge transfer)', 
             fontsize=9, color='#555555', ha='center', style='italic')

    ax2.set_title('Physical Mechanism', fontweight='bold', pad=10)

    # --- Right: Per-reaction residual heatmap (simplified) ---
    ax3 = fig.add_subplot(gs[2])
    ax3.set_facecolor(C_BG)

    reactions = ['RXN_0b71', 'RXN_0d5f', 'RXN_1be0', 'RXN_358d', 'RXN_8c09', 'RXN_8ebf', 'RXN_a18e', 'RXN_a204', 'RXN_a2b7', 'RXN_c29a', 'RXN_f0ec', 'RXN_f5d5']
    dr_resid = [0.03, -0.16, -0.06, 0.25, -0.07, -0.34, 0.16, 0.05, -0.08, 0.16, 0.18, -0.15]
    y_resid = [0.06, -0.16, -0.06, 0.25, 0.07, -0.11, 0.16, 0.05, -0.08, 0.16, 0.18, -0.15]

    # Color by sign
    colors_dr = [C_POS if r > 0 else C_NEG for r in dr_resid]
    colors_y = [C_POS if r > 0 else C_NEG for r in y_resid]

    ypos = np.arange(len(reactions))
    ax3.barh(ypos - 0.2, dr_resid, 0.35, color=colors_dr, alpha=0.7, label='DR resid', zorder=3)
    ax3.barh(ypos + 0.2, y_resid, 0.35, color=colors_y, alpha=0.5, edgecolor=colors_y, linewidth=1.5, label='Yield resid', zorder=3)

    ax3.axvline(x=0, color='#999999', linewidth=0.8, zorder=1)
    ax3.set_yticks(ypos)
    ax3.set_yticklabels(reactions, fontsize=7)
    ax3.set_xlabel('Signed Residual', fontsize=9)
    ax3.set_title('Per-Reaction Residuals', fontweight='bold')
    ax3.legend(loc='lower right', fontsize=7)
    ax3.set_xlim(-0.5, 0.5)
    ax3.grid(axis='x', zorder=0)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    fig.suptitle('Fig 2. DR vs Yield: Final Best Models & Physical Drivers', fontweight='bold', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig2_dr_vs_yield.png', facecolor='white')
    plt.close()
    print("✓ fig2_dr_vs_yield.png")

# ---------------------------------------------------------------------------
# Figure 3: V13 Yield Key Ablation Experiments
# ---------------------------------------------------------------------------
def fig3_ablation():
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    fig.patch.set_facecolor('white')

    # --- Panel A: GEDT Ablation ---
    ax = axes[0]
    ax.set_facecolor(C_BG)
    models = ['Y-A\n(Full)', 'Y-B\n(-GEDT)', 'Y-C\n(GEDT only)']
    r2_vals = [0.503, -0.885, -0.107]
    colors = [C_POS, C_NEG, C_NEU]
    bars = ax.bar(models, r2_vals, color=colors, edgecolor='white', linewidth=0.5, zorder=3)
    for bar, val in zip(bars, r2_vals):
        ypos = val + 0.04 if val >= 0 else val - 0.10
        ax.text(bar.get_x() + bar.get_width()/2., ypos, f'{val:.3f}', ha='center', va='bottom' if val>=0 else 'top', fontsize=10, fontweight='bold')
    ax.axhline(y=0, color='#999999', linewidth=0.8, zorder=1)
    ax.set_ylabel('LORO R²', fontweight='bold')
    ax.set_title('A. GEDT Ablation', fontweight='bold')
    ax.set_ylim(-1.1, 0.7)
    ax.annotate('ΔR² = +1.39', xy=(0.5, -0.19), fontsize=10, fontweight='bold', color=C_GOLD, ha='center')
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- Panel B: Stability-2 Synergy ---
    ax = axes[1]
    ax.set_facecolor(C_BG)
    models = ['B0\n(F2 only)', 'B1\n(+ccd)', 'B2\n(+post)', 'B3\n(+both)']
    r2_vals = [0.506, 0.498, 0.497, 0.561]
    colors = [C_NEU, C_NEG, C_NEG, C_POS]
    bars = ax.bar(models, r2_vals, color=colors, edgecolor='white', linewidth=0.5, zorder=3)
    for bar, val in zip(bars, r2_vals):
        ax.text(bar.get_x() + bar.get_width()/2., val + 0.01, f'{val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.axhline(y=0.506, color=C_NEU, linewidth=1, linestyle='--', zorder=1, alpha=0.7)
    ax.annotate('Baseline', xy=(2, 0.51), fontsize=8, color=C_NEU)
    ax.annotate('Synergy\nΔR²=+0.063', xy=(3, 0.561), xytext=(2.5, 0.62),
                arrowprops=dict(arrowstyle='->', color=C_POS, lw=1.5),
                fontsize=9, fontweight='bold', color=C_POS, ha='center')
    ax.set_ylabel('LORO R²', fontweight='bold')
    ax.set_title('B. Stability-2 Synergy', fontweight='bold')
    ax.set_ylim(0.45, 0.65)
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- Panel C: GEDT × Geometry Interaction ---
    ax = axes[2]
    ax.set_facecolor(C_BG)
    models = ['C0\n(baseline)', 'C1\n(×ccd)', 'C2\n(×post)']
    r2_vals = [0.555, 0.480, 0.529]
    colors = [C_POS, C_NEG, C_NEG]
    bars = ax.bar(models, r2_vals, color=colors, edgecolor='white', linewidth=0.5, zorder=3)
    for bar, val in zip(bars, r2_vals):
        ax.text(bar.get_x() + bar.get_width()/2., val + 0.005, f'{val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')
    ax.axhline(y=0.555, color=C_POS, linewidth=1, linestyle='--', zorder=1, alpha=0.7)
    ax.annotate('No multiplicative\ninteraction', xy=(1.5, 0.51), fontsize=9, color=C_NEG, ha='center', fontweight='bold')
    ax.set_ylabel('LORO R²', fontweight='bold')
    ax.set_title('C. GEDT × Geometry Interaction', fontweight='bold')
    ax.set_ylim(0.45, 0.60)
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.suptitle('Fig 3. V13 Yield — Key Mechanistic Ablation Experiments', fontweight='bold', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig3_ablation.png', facecolor='white')
    plt.close()
    print("✓ fig3_ablation.png")

# ---------------------------------------------------------------------------
# Figure 4: V13-Yield Linear Surrogate Coefficients
# ---------------------------------------------------------------------------
def fig4_coefficients():
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_facecolor(C_BG)

    features = [
        'GEDT (major_s2_gedt)',
        'Temperature (K)',
        'ΔS1 E_span',
        'S1 E_span (major)',
        'S1 Sconf (major)',
        'ccd (major_int)',
        'post_bond_q10 (major_s1)'
    ]
    coeffs = [0.500, 0.203, 0.195, -0.154, 0.090, -0.070, 0.040]
    colors = [C_POS if c > 0 else C_NEG for c in coeffs]

    ypos = np.arange(len(features))
    bars = ax.barh(ypos, coeffs, color=colors, edgecolor='white', height=0.6, zorder=3)

    for bar, val in zip(bars, coeffs):
        xpos = val + 0.015 if val > 0 else val - 0.015
        ax.text(xpos, bar.get_y() + bar.get_height()/2., f'{val:+.3f}', va='center', ha='left' if val>0 else 'right', fontsize=9, fontweight='bold')

    ax.axvline(x=0, color='#999999', linewidth=0.8, zorder=1)
    ax.set_yticks(ypos)
    ax.set_yticklabels(features, fontsize=9)
    ax.set_xlabel('BayesianRidge Coefficient (standardized)', fontweight='bold')
    ax.set_title('Fig 4. V13 Yield — Linear Surrogate Model Coefficients\n(Physical Interpretation of Feature Contributions)', fontweight='bold', pad=12)

    # Add interpretation annotations
    ax.annotate('Electronic transfer\n→ promotes yield', xy=(0.50, 6), xytext=(0.35, 6.5),
                fontsize=8, color=C_POS, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=C_POS, lw=1.2))
    ax.annotate('Steric crowding\n→ penalizes yield', xy=(-0.07, 2), xytext=(-0.22, 3.5),
                fontsize=8, color=C_NEG, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color=C_NEG, lw=1.2))

    ax.set_xlim(-0.30, 0.60)
    ax.grid(axis='x', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig4_coefficients.png', facecolor='white')
    plt.close()
    print("✓ fig4_coefficients.png")

# ---------------------------------------------------------------------------
# Figure 5: V13-DR PMPF & TS-lite Results
# ---------------------------------------------------------------------------
def fig5_pmpf_tslite():
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    fig.patch.set_facecolor('white')

    # --- Panel A: PMPF Models ---
    ax = axes[0]
    ax.set_facecolor(C_BG)
    models = ['PMPF-A\n(ccd)', 'PMPF-B\n(ccd+preorg)', 'PMPF-C\n(ccd+preorg+Nconf)', 'M1\nStability-2']
    r2_vals = [0.618, 0.483, 0.675, 0.700]
    colors = [C_NEG, C_NEG, C_NEG, C_POS]
    bars = ax.bar(models, r2_vals, color=colors, edgecolor='white', linewidth=0.5, zorder=3)
    for bar, val in zip(bars, r2_vals):
        ax.text(bar.get_x() + bar.get_width()/2., val + 0.01, f'{val:.3f}', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.axhline(y=0.72, color=C_GOLD, linewidth=2, linestyle='--', zorder=1, label='Acceptance gate (R²>0.72)')
    ax.axhline(y=0.70, color=C_POS, linewidth=1.5, linestyle='-', zorder=1, alpha=0.5)
    ax.annotate('M1\n(best)', xy=(3, 0.70), xytext=(2.3, 0.78), fontsize=9, fontweight='bold', color=C_POS,
                arrowprops=dict(arrowstyle='->', color=C_POS, lw=1.2))

    ax.set_ylabel('LORO R²', fontweight='bold')
    ax.set_title('A. PMPF Model Comparison', fontweight='bold')
    ax.set_ylim(0.35, 0.85)
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- Panel B: TS-lite Layer Comparison ---
    ax = axes[1]
    ax.set_facecolor(C_BG)
    layers = ['TS-free\n(S1 only)', 'Intermediate\nproxy', 'TS-geom', 'TS-electronic', 'TS-energy', 'M1\n(baseline)']
    # Log scale for extreme negative values
    r2_vals = [-824, -129648, -3.79, -15.4, -3.33, 0.70]
    # For visualization, clip extreme negatives and use color coding
    r2_display = [max(v, -20) for v in r2_vals]
    colors = [C_NEG if v < 0 else C_POS for v in r2_vals]

    bars = ax.bar(layers, r2_display, color=colors, edgecolor='white', linewidth=0.5, zorder=3)
    for bar, val, disp in zip(bars, r2_vals, r2_display):
        label = f'{val:.1f}' if abs(val) < 1000 else f'{val:.0e}'
        ypos = disp + 0.5 if disp < 0 else disp + 0.5
        ax.text(bar.get_x() + bar.get_width()/2., ypos, label, ha='center', va='bottom', fontsize=8, fontweight='bold', rotation=30)

    ax.axhline(y=0, color='#999999', linewidth=0.8, zorder=1)
    ax.axhline(y=0.70, color=C_POS, linewidth=1.5, linestyle='--', zorder=1, alpha=0.5)
    ax.set_ylabel('LORO R² (clipped at -20)', fontweight='bold')
    ax.set_title('B. TS-lite Layer Comparison', fontweight='bold')
    ax.set_ylim(-25, 5)
    ax.annotate('All TS-lite models\nshow NEGATIVE R²', xy=(2, -15), fontsize=10, fontweight='bold', color=C_NEG, ha='center')
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.suptitle('Fig 5. V13 DR — PMPF & TS-lite Model Comparison', fontweight='bold', fontsize=13, y=1.02)
    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig5_pmpf_tslite.png', facecolor='white')
    plt.close()
    print("✓ fig5_pmpf_tslite.png")

# ---------------------------------------------------------------------------
# Figure 6: V14 Feature Mining Block Screening
# ---------------------------------------------------------------------------
def fig6_feature_mining():
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_facecolor(C_BG)

    blocks = [
        'A0\nBaseline\n(8feat)',
        'A1\nINT_E\n(+2feat)',
        'A2\nINT_GEO\n(+2feat)',
        'A3\nINT_POL\n(+2feat)',
        'A4\nTS_ELEC_EFF\n(+2feat)',
        'A5\nCONF_POP\n(+4feat)'
    ]
    r2_vals = [0.561, 0.442, 0.520, 0.513, 0.539, 0.484]
    gates_pass = [4, 3, 3, 2, 5, 2]  # out of 7

    x = np.arange(len(blocks))
    bar_w = 0.55
    colors = [C_POS if g >= 4 else C_NEG for g in gates_pass]
    bars = ax.bar(x, r2_vals, bar_w, color=colors, edgecolor='white', linewidth=0.5, zorder=3)

    for bar, val, g in zip(bars, r2_vals, gates_pass):
        ax.text(bar.get_x() + bar.get_width()/2., val + 0.008, f'R²={val:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')
        ax.text(bar.get_x() + bar.get_width()/2., val - 0.025, f'{g}/7 gates', ha='center', va='top', fontsize=7, color='white', fontweight='bold')

    ax.axhline(y=0.561, color=C_NEU, linewidth=1, linestyle='--', zorder=1, alpha=0.7)
    ax.annotate('Baseline', xy=(0, 0.565), fontsize=8, color=C_NEU)

    ax.annotate('Best block\n(5/7 gates)', xy=(4, 0.539), xytext=(4.5, 0.58),
                arrowprops=dict(arrowstyle='->', color=C_POS, lw=1.5),
                fontsize=9, fontweight='bold', color=C_POS, ha='center')

    ax.set_ylabel('LORO R²', fontweight='bold')
    ax.set_xlabel('Feature Block', fontweight='bold')
    ax.set_title('Fig 6. V14 Yield Feature Mining — Block Screening Results\n(No block surpasses Y-F baseline; GEDT_eff_ravg selected 66.7% in nested CV)', fontweight='bold', pad=12)
    ax.set_xticks(x)
    ax.set_xticklabels(blocks, fontsize=8)
    ax.set_ylim(0.40, 0.62)
    ax.grid(axis='y', zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor=C_POS, label='≥4/7 gates pass'),
                       Patch(facecolor=C_NEG, label='<4/7 gates pass')]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig6_feature_mining.png', facecolor='white')
    plt.close()
    print("✓ fig6_feature_mining.png")

# ---------------------------------------------------------------------------
# Figure 7: Dual-Mechanism Physical Picture (Schematic)
# ---------------------------------------------------------------------------
def fig7_dual_mechanism():
    fig, ax = plt.subplots(figsize=(11, 6))
    ax.set_facecolor(C_BG)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Title
    ax.text(5, 9.5, 'Fig 7. Dual-Mechanism Physical Picture', fontsize=14, fontweight='bold', ha='center', color='#333333')
    ax.text(5, 9.0, '[4+3] Cycloaddition: DR and Yield are Independent Physical Problems', fontsize=10, ha='center', color='#666666', style='italic')

    # Reactant (left)
    rect_react = mpatches.FancyBboxPatch((0.3, 3.5), 1.8, 3.5, boxstyle="round,pad=0.1", facecolor='#E8F4FD', edgecolor=C_DR, linewidth=2)
    ax.add_patch(rect_react)
    ax.text(1.2, 6.5, 'Reactants', fontsize=10, fontweight='bold', color=C_DR, ha='center')
    ax.text(1.2, 5.8, 'Allenamide\n+ DMDO', fontsize=8, color='#555555', ha='center')
    ax.text(1.2, 4.8, 'S1 Conformer\nEnsemble', fontsize=8, color='#777777', ha='center', style='italic')

    # S1 box
    rect_s1 = mpatches.FancyBboxPatch((2.8, 5.2), 1.8, 2.5, boxstyle="round,pad=0.08", facecolor='#FFF3E0', edgecolor=C_GOLD, linewidth=2)
    ax.add_patch(rect_s1)
    ax.text(3.7, 7.2, 'S1: Preorganization', fontsize=9, fontweight='bold', color=C_GOLD, ha='center')
    ax.text(3.7, 6.6, '• E_span\n• Nconf_eff\n• Sconf\n• post_bond_q10', fontsize=7, color='#555555', ha='center', family='monospace')

    # S2/S3 box
    rect_ts = mpatches.FancyBboxPatch((5.3, 5.2), 1.8, 2.5, boxstyle="round,pad=0.08", facecolor='#FCE4EC', edgecolor=C_YIELD, linewidth=2)
    ax.add_patch(rect_ts)
    ax.text(6.2, 7.2, 'TS: Charge Transfer', fontsize=9, fontweight='bold', color=C_YIELD, ha='center')
    ax.text(6.2, 6.6, '• GEDT\n• CDFT μ\n• ω (electro-\nphilicity)', fontsize=7, color='#555555', ha='center', family='monospace')

    # Product (right)
    rect_prod = mpatches.FancyBboxPatch((7.8, 3.5), 1.8, 3.5, boxstyle="round,pad=0.1", facecolor='#E8F5E9', edgecolor=C_POS, linewidth=2)
    ax.add_patch(rect_prod)
    ax.text(8.7, 6.5, 'Product', fontsize=10, fontweight='bold', color=C_POS, ha='center')
    ax.text(8.7, 5.8, '7-membered\nRing', fontsize=8, color='#555555', ha='center')

    # Arrows
    ax.annotate('', xy=(2.8, 6.5), xytext=(2.1, 6.5), arrowprops=dict(arrowstyle='->', color=C_DR, lw=2))
    ax.annotate('', xy=(5.3, 6.5), xytext=(4.6, 6.5), arrowprops=dict(arrowstyle='->', color=C_YIELD, lw=2))
    ax.annotate('', xy=(7.8, 6.5), xytext=(7.1, 6.5), arrowprops=dict(arrowstyle='->', color=C_POS, lw=2))

    # DR prediction box (top-left pathway)
    ax.annotate('DR Selectivity', xy=(3.7, 5.0), xytext=(2.5, 3.0),
                fontsize=10, fontweight='bold', color=C_DR, ha='center',
                arrowprops=dict(arrowstyle='->', color=C_DR, lw=2, connectionstyle='arc3,rad=0.3'))
    ax.text(2.5, 2.5, 'R² = 0.700\nDriver: S1 Preorganization\n→ What shape is the TS?', fontsize=8, color=C_DR, ha='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=C_DR, alpha=0.1, edgecolor=C_DR))

    # Yield prediction box (top-right pathway)
    ax.annotate('Yield', xy=(6.2, 5.0), xytext=(7.5, 3.0),
                fontsize=10, fontweight='bold', color=C_YIELD, ha='center',
                arrowprops=dict(arrowstyle='->', color=C_YIELD, lw=2, connectionstyle='arc3,rad=-0.3'))
    ax.text(7.5, 2.5, 'R² = 0.550\nDriver: TS GEDT\n→ Can the TS form?', fontsize=8, color=C_YIELD, ha='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor=C_YIELD, alpha=0.1, edgecolor=C_YIELD))

    # Independence label
    ax.annotate('', xy=(6.0, 4.5), xytext=(4.0, 4.5),
                arrowprops=dict(arrowstyle='<->', color=C_NEU, lw=2, linestyle='--'))
    ax.text(5.0, 4.2, 'Statistically Independent\n(ρ = -0.08, p = 0.73)', fontsize=9, color=C_NEU, ha='center', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor=C_NEU, alpha=0.8))

    plt.tight_layout()
    plt.savefig(f'{OUTPUT_DIR}/fig7_dual_mechanism.png', facecolor='white')
    plt.close()
    print("✓ fig7_dual_mechanism.png")

# ---------------------------------------------------------------------------
# Figure 8: Summary Dashboard (All-in-one)
# ---------------------------------------------------------------------------
def fig8_dashboard():
    fig = plt.figure(figsize=(14, 10))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.35)
    fig.patch.set_facecolor('white')

    # Title
    fig.suptitle('RPH [4+3] Cycloaddition ML — V1–V14 Complete Results Dashboard', 
                 fontsize=15, fontweight='bold', y=0.98)

    # --- (0,0): Timeline mini ---
    ax = fig.add_subplot(gs[0, 0])
    ax.set_facecolor(C_BG)
    rounds = ['V6', 'V8', 'V10', 'V11', 'V12', 'V13', 'V14']
    dr_line = [0.177, 0.431, 0.276, 0.657, 0.700, 0.700, 0.700]
    y_line = [0.177, 0.431, 0.276, -0.982, 0.550, 0.550, 0.449]
    ax.plot(rounds, dr_line, 'o-', color=C_DR, linewidth=2, markersize=6, label='DR')
    ax.plot(rounds, y_line, 's--', color=C_YIELD, linewidth=2, markersize=6, label='Yield')
    ax.axhline(y=0, color='#999', linewidth=0.5)
    ax.set_title('R² Evolution', fontweight='bold', fontsize=10)
    ax.set_ylabel('LORO R²')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- (0,1): GEDT Ablation mini ---
    ax = fig.add_subplot(gs[0, 1])
    ax.set_facecolor(C_BG)
    models = ['Full', '-GEDT', 'GEDT\nonly']
    vals = [0.503, -0.885, -0.107]
    colors = [C_POS, C_NEG, C_NEU]
    ax.bar(models, vals, color=colors, edgecolor='white', zorder=3)
    ax.axhline(y=0, color='#999', linewidth=0.5)
    ax.set_title('GEDT Ablation', fontweight='bold', fontsize=10)
    ax.set_ylabel('LORO R²')
    ax.text(0.5, -0.3, 'Δ=+1.39', fontsize=9, fontweight='bold', color=C_GOLD, ha='center')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- (0,2): DR Model Comparison ---
    ax = fig.add_subplot(gs[0, 2])
    ax.set_facecolor(C_BG)
    models = ['M1', 'PMPF-A', 'PMPF-B', 'PMPF-C']
    vals = [0.700, 0.618, 0.483, 0.675]
    colors = [C_POS, C_NEG, C_NEG, C_NEG]
    ax.barh(models, vals, color=colors, edgecolor='white', height=0.5, zorder=3)
    for i, v in enumerate(vals):
        ax.text(v + 0.01, i, f'{v:.3f}', va='center', fontsize=8, fontweight='bold')
    ax.axvline(x=0.72, color=C_GOLD, linewidth=1.5, linestyle='--', label='Gate')
    ax.set_xlim(0, 0.85)
    ax.set_title('DR PMPF Models', fontweight='bold', fontsize=10)
    ax.legend(fontsize=7)
    ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- (1,0): Feature Importance (Yield) ---
    ax = fig.add_subplot(gs[1, 0])
    ax.set_facecolor(C_BG)
    feats = ['GEDT', 'T', 'ΔE_span', 'E_span', 'Sconf', 'ccd', 'post_bond']
    coeffs = [0.50, 0.20, 0.19, -0.15, 0.09, -0.07, 0.04]
    colors = [C_POS if c > 0 else C_NEG for c in coeffs]
    ypos = np.arange(len(feats))
    ax.barh(ypos, coeffs, color=colors, height=0.5, zorder=3)
    ax.set_yticks(ypos)
    ax.set_yticklabels(feats, fontsize=8)
    ax.axvline(x=0, color='#999', linewidth=0.5)
    ax.set_title('Yield Coefficients', fontweight='bold', fontsize=10)
    ax.set_xlabel('β')
    ax.grid(axis='x', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.invert_yaxis()

    # --- (1,1): Stability-2 Synergy ---
    ax = fig.add_subplot(gs[1, 1])
    ax.set_facecolor(C_BG)
    models = ['F2', '+ccd', '+post', '+both']
    vals = [0.506, 0.498, 0.497, 0.561]
    colors = [C_NEU, C_NEG, C_NEG, C_POS]
    bars = ax.bar(models, vals, color=colors, edgecolor='white', zorder=3)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2., v + 0.005, f'{v:.3f}', ha='center', fontsize=8, fontweight='bold')
    ax.set_ylim(0.45, 0.60)
    ax.set_title('Stability-2 Synergy', fontweight='bold', fontsize=10)
    ax.set_ylabel('LORO R²')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- (1,2): V14 Block Screening ---
    ax = fig.add_subplot(gs[1, 2])
    ax.set_facecolor(C_BG)
    blocks = ['A0', 'A1', 'A2', 'A3', 'A4', 'A5']
    vals = [0.561, 0.442, 0.520, 0.513, 0.539, 0.484]
    gates = [4, 3, 3, 2, 5, 2]
    colors = [C_POS if g >= 4 else C_NEG for g in gates]
    ax.bar(blocks, vals, color=colors, edgecolor='white', zorder=3)
    ax.axhline(y=0.561, color=C_NEU, linewidth=1, linestyle='--', alpha=0.5)
    ax.set_ylim(0.40, 0.60)
    ax.set_title('V14 Block Screening', fontweight='bold', fontsize=10)
    ax.set_ylabel('LORO R²')
    ax.grid(axis='y', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # --- (2,0-2): Key numbers summary ---
    ax = fig.add_subplot(gs[2, :])
    ax.set_facecolor('#F5F5F5')
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis('off')

    # Summary boxes
    summaries = [
        ('Dataset', '20 rows\n12 reactions\n16 DR rows', C_DR),
        ('DR Model', 'M1 Stability-2\nR² = 0.700\nMAE = 0.254', C_DR),
        ('Yield Model', 'Y-F GBR\nR² = 0.550\nMAE = 0.135', C_YIELD),
        ('Key Finding', 'DR & Yield\nare INDEPENDENT\n(ρ ≈ 0)', C_GOLD),
        ('Physical', 'DR: Preorg\nYield: GEDT\nDual-mechanism', C_POS),
        ('Validation', 'LORO + Perm\n+ Bootstrap\n+ Y-random', C_PURPLE),
    ]

    for i, (title, text, color) in enumerate(summaries):
        x = 0.5 + i * 1.6
        rect = mpatches.FancyBboxPatch((x, 0.3), 1.4, 2.4, boxstyle="round,pad=0.08", 
                                        facecolor=color, alpha=0.08, edgecolor=color, linewidth=2)
        ax.add_patch(rect)
        ax.text(x + 0.7, 2.3, title, fontsize=9, fontweight='bold', color=color, ha='center')
        ax.text(x + 0.7, 1.3, text, fontsize=8, color='#444444', ha='center', va='center')

    plt.savefig(f'{OUTPUT_DIR}/fig8_dashboard.png', facecolor='white')
    plt.close()
    print("✓ fig8_dashboard.png")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == '__main__':
    print("=" * 60)
    print("RPH V1-V14 Visualization Report Generator")
    print("=" * 60)
    fig1_evolution()
    fig2_dr_vs_yield()
    fig3_ablation()
    fig4_coefficients()
    fig5_pmpf_tslite()
    fig6_feature_mining()
    fig7_dual_mechanism()
    fig8_dashboard()
    print("=" * 60)
    print(f"All figures saved to: {OUTPUT_DIR}")
    print("=" * 60)
