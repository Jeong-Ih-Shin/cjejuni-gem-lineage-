"""
Step 11: Focused ST-443 vs ST-658 comparison on human gut diet media (Figure 4).

Inputs:
    intermediate/human_gut_simulation_results.csv  (from step 07)

Outputs:
    intermediate/st443_vs_st658_statistics.csv
    output/figures/Figure4_st443_vs_st658.{png, svg}

Methods:
    Compare lineage-mean mu between ST-443 (n = 7) and ST-658 (n = 4) using
    a one-sided Mann-Whitney U test (a priori hypothesis: ST-443 < ST-658,
    based on ST-443's pantothenate/CoA biosynthesis deficit identified in
    step 02). Effect size is quantified by Cohen's d.

    Two diet media (from step 07):
      - Western diet
      - High-fiber diet
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # WSL2 / headless
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats as scipy_stats

# ============================================================
# Paths
# ============================================================
REPO_ROOT        = Path(__file__).resolve().parent.parent
INTERMEDIATE_DIR = REPO_ROOT / 'intermediate'
FIGURES_DIR      = REPO_ROOT / 'output' / 'figures'
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SIM_CSV   = INTERMEDIATE_DIR / 'human_gut_simulation_results.csv'
STATS_CSV = INTERMEDIATE_DIR / 'st443_vs_st658_statistics.csv'
FIG_PNG   = FIGURES_DIR / 'Figure4_st443_vs_st658.png'
FIG_SVG   = FIGURES_DIR / 'Figure4_st443_vs_st658.svg'


# ============================================================
# Configuration  (matches final approved figure styling)
# ============================================================
MEDIA = [
    ('mu_western_diet',    'Western diet'),
    ('mu_high_fiber_diet', 'High-fiber diet'),
]

LINEAGES = [
    ('ST-443 complex', 'ST-443', '#e67e22'),  # orange
    ('ST-658 complex', 'ST-658', '#27ae60'),  # green
]

FIG_W = 1.67    # width per panel (final approved)
FIG_H = 3.5     # height


# ============================================================
# Statistics
# ============================================================
def stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'


def cohen_d(x, y):
    nx, ny = len(x), len(y)
    if nx < 2 or ny < 2:
        return np.nan
    s = np.sqrt(((nx - 1) * np.var(x, ddof=1) +
                 (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    if s == 0:
        return np.nan
    return (np.mean(x) - np.mean(y)) / s


def compute_stats(df):
    rows = []
    for col, label in MEDIA:
        x443 = df[df['CC'] == 'ST-443 complex'][col].dropna().values
        x658 = df[df['CC'] == 'ST-658 complex'][col].dropna().values
        u, p = scipy_stats.mannwhitneyu(x443, x658, alternative='less')
        d = cohen_d(x443, x658)
        rows.append({
            'medium':           label,
            'st443_n':          len(x443),
            'st443_mean':       float(np.mean(x443)),
            'st443_sd':         float(np.std(x443, ddof=1)) if len(x443) > 1 else 0.0,
            'st658_n':          len(x658),
            'st658_mean':       float(np.mean(x658)),
            'st658_sd':         float(np.std(x658, ddof=1)) if len(x658) > 1 else 0.0,
            'mannwhitney_u':    float(u),
            'p_value_onesided': float(p),
            'cohens_d':         float(d),
            'significance':     stars(p),
        })
    stat_df = pd.DataFrame(rows)
    stat_df.to_csv(STATS_CSV, index=False)
    print(f"✓ saved {STATS_CSV}")
    return stat_df


# ============================================================
# Figure 4 — 2-panel bar chart (final approved styling)
# ============================================================
def annotate(ax, sub, row):
    """Bracket + significance stars + p/d annotation."""
    ymax = sub['mu'].max()
    ymin = sub['mu'].min()
    yrange = ymax - ymin if ymax > ymin else 0.01
    bracket_y = ymax + yrange * 0.18
    ax.plot([0, 0, 1, 1],
            [bracket_y, bracket_y + yrange * 0.04,
             bracket_y + yrange * 0.04, bracket_y],
            color='black', lw=1)
    ax.text(0.5, bracket_y + yrange * 0.07, row['significance'],
            ha='center', va='bottom', fontsize=11, fontweight='bold')
    ax.text(0.5, bracket_y + yrange * 0.55,
            f"p = {row['p_value_onesided']:.3f}\nd = {row['cohens_d']:+.2f}",
            ha='center', va='bottom', fontsize=8)
    ax.set_ylim(bottom=0.27, top=ymax + yrange * 2.5)


def plot_figure4(df, stat_df):
    fig, axes = plt.subplots(1, 2, figsize=(FIG_W * 2, FIG_H), sharey=True)

    for ax, (col, label) in zip(axes, MEDIA):
        # Long-form sub-DataFrame for this medium
        sub_rows = []
        for cc_full, _, _ in LINEAGES:
            vals = df[df['CC'] == cc_full][col].dropna().values
            for v in vals:
                sub_rows.append({'CC': cc_full, 'mu': v})
        sub = pd.DataFrame(sub_rows)

        # Bar chart (mean ± SD)
        means = sub.groupby('CC')['mu'].mean().reindex(
            [cc for cc, _, _ in LINEAGES])
        sds = sub.groupby('CC')['mu'].std(ddof=1).reindex(
            [cc for cc, _, _ in LINEAGES]).fillna(0)
        colors = [c for _, _, c in LINEAGES]
        ax.bar(range(len(LINEAGES)), means.values, yerr=sds.values,
               color=colors, alpha=0.6, capsize=6,
               edgecolor='black', linewidth=1, width=0.6)

        # Individual strains (swarm, black)
        sns.swarmplot(data=sub, x='CC', y='mu',
                      order=[cc for cc, _, _ in LINEAGES],
                      palette={cc: 'black' for cc, _, _ in LINEAGES},
                      ax=ax, size=5, alpha=0.85)

        # Annotation
        row = stat_df[stat_df['medium'] == label].iloc[0]
        annotate(ax, sub, row)

        # n labels
        ns = {cc: int((sub['CC'] == cc).sum()) for cc, _, _ in LINEAGES}
        ax.set_xticklabels(
            [f'{short}\n(n={ns[cc]})' for cc, short, _ in LINEAGES]
        )
        ax.set_xlabel('')
        ax.set_ylabel('Simulated μ (1/h)' if ax is axes[0] else '')
        ax.set_title(label, fontsize=9, pad=22)
        ax.spines[['top', 'right']].set_visible(False)

    plt.tight_layout()
    fig.savefig(FIG_PNG, dpi=300, bbox_inches='tight')
    fig.savefig(FIG_SVG, bbox_inches='tight')
    print(f"✓ saved {FIG_PNG}")
    print(f"✓ saved {FIG_SVG}")
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main():
    df = pd.read_csv(SIM_CSV)
    print(f"Loaded {len(df)} strains from {SIM_CSV.name}")
    for col, label in MEDIA:
        if col not in df.columns:
            raise ValueError(
                f"Column '{col}' missing from {SIM_CSV.name}. "
                f"Did step 07 finish for both diets?"
            )
        valid = df[col].dropna()
        print(f"  {label}: {len(valid)} valid, mu {valid.min():.4f} – {valid.max():.4f}")
    print()

    stat_df = compute_stats(df)
    print("\nStatistics:")
    print(stat_df.to_string(index=False))
    print()

    plot_figure4(df, stat_df)


if __name__ == '__main__':
    main()
