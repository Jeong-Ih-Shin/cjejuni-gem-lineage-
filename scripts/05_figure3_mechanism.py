"""
Step 05 — Figure 3: Pathway-level GEM mechanism (ST-443 LOST vs ST-658 GAINED).

Inputs:
    intermediate/ST443_vs_others_enrichment.csv  (from step 02)
    intermediate/ST658_vs_others_enrichment.csv  (from step 02)

Outputs:
    output/figures/Figure3_mechanism.{png,svg}

Reactions are classified into functional pathway categories by keyword
matching on the reaction name. This is a heuristic — verify against
ModelSEED/KEGG annotation for any specific claim.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import textwrap

sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import INTERMEDIATE_DIR, NPG, SEM, apply_style, save_fig

apply_style()


def classify_443(name):
    n = name.lower()
    if any(k in n for k in ['pantoate', 'pantothen', 'methylenetetrahydrofolate',
                            '3-methyl-2-oxobutanoate', 'beta-alanine', '3-oxopropanoate']):
        return 'Pantothenate / CoA biosynthesis'
    if 'betaine' in n or 'bet-' in n:
        return 'Glycine betaine / osmolyte uptake'
    if 'glycerol' in n:
        return 'Glycerol metabolism / transport'
    if any(k in n for k in ['fucose', 'gdp-4-dehydro-6-deoxy', 'gdp-mannose']):
        return 'Fucose / GDP-mannose biosynthesis'
    if ('oxidoreductase(deaminating)' in n or 'oxidoreductase (deaminating)' in n
            or 'putrescine' in n):
        return 'Monoamine / polyamine oxidase'
    if 'udp-' in n:
        return 'Nucleotide sugar metabolism'
    return 'Other (unclassified)'


def classify_658(name):
    n = name.lower()
    if 'acetylneuraminate' in n or 'neuraminate' in n:
        return 'Sialic acid (N-acetylneuraminate)'
    if 'sulfide' in n or 'thiosulfate' in n:
        return 'Anaerobic sulfur respiration'
    if any(k in n for k in ['phosphatidyl', 'phospholipase', 'cardiolipin', 'acyl',
                            'ketoacyl', 'malonyl', 'palmitoyl', 'hexenoyl', 'fatty acid']):
        return 'Fatty acid / phospholipid biosynthesis'
    if 'sodium' in n or 'na+' in n or 'antiport' in n:
        return 'Sodium-coupled transport'
    if 'glucosamine' in n or 'ribose' in n:
        return 'Sugar / amino-sugar transport'
    return 'Other (unclassified)'


def aggregate(df):
    a = df.groupby('cat').agg(n=('reaction', 'count'),
                              min_p=('p', 'min'),
                              mean_pct=('pct_diff', 'mean')).reset_index()
    a = a.sort_values('n', ascending=True)
    other = a[a['cat'].str.startswith('Other')]
    named = a[~a['cat'].str.startswith('Other')]
    return pd.concat([other, named], ignore_index=True)


def fmt_p(p):
    for thr, label in [(1e-6,'$p<10^{-6}$'),(1e-5,'$p<10^{-5}$'),(1e-4,'$p<10^{-4}$'),
                       (1e-3,'$p<10^{-3}$'),(1e-2,'p<0.01')]:
        if p < thr: return label
    return f'p={p:.2g}'


def main():
    st443 = pd.read_csv(INTERMEDIATE_DIR / 'ST443_vs_others_enrichment.csv')
    st658 = pd.read_csv(INTERMEDIATE_DIR / 'ST658_vs_others_enrichment.csv')

    st443_sig = st443[(st443['pct_diff'] < -50) & (st443['p'] < 0.01)].copy()
    st658_sig = st658[(st658['pct_diff'] > +50) & (st658['p'] < 0.01)].copy()
    print(f'ST-443 LOST reactions (sig): {len(st443_sig)}')
    print(f'ST-658 GAINED reactions (sig): {len(st658_sig)}')

    st443_sig['cat'] = st443_sig['name'].apply(classify_443)
    st658_sig['cat'] = st658_sig['name'].apply(classify_658)
    agg443 = aggregate(st443_sig)
    agg658 = aggregate(st658_sig)
    print('\nST-443:'); print(agg443.to_string(index=False))
    print('\nST-658:'); print(agg658.to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(6.69, 3.5))
    plt.subplots_adjust(top=0.86, bottom=0.13, left=0.04, right=0.98, wspace=0.55)

    # A — ST-443 LOST
    ax = axes[0]
    ax.grid(axis='x', alpha=0.5, color=SEM['grid'], linewidth=0.5)
    y_pos = np.arange(len(agg443))
    colors = [SEM['lost'] if not c.startswith('Other') else '#D5D5D5' for c in agg443['cat']]
    ax.barh(y_pos, agg443['n'], color=colors, edgecolor='white', linewidth=1.5, height=0.72)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([textwrap.fill(s, width=20) for s in agg443['cat']], fontsize=7)
    ax.set_xlabel('Number of reactions LOST (vs other strains)')
    ax.set_title('A', loc='left', x=-0.55, fontsize=11, fontweight='bold')
    for i, (n, p) in enumerate(zip(agg443['n'], agg443['min_p'])):
        ax.text(n + 0.3, i, f'{n}  ({fmt_p(p)})', va='center', fontsize=8, color=SEM['text'])
    ax.set_xlim(0, agg443['n'].max() * 1.45 if len(agg443) else 1)

    # B — ST-658 GAINED with sulfur respiration highlight
    ax = axes[1]
    ax.grid(axis='x', alpha=0.5, color=SEM['grid'], linewidth=0.5)
    y_pos = np.arange(len(agg658))
    sresp_idx = None
    colors = [SEM['gained'] if not c.startswith('Other') else '#D5D5D5' for c in agg658['cat']]
    edge_colors = [NPG['red'] if i == sresp_idx else 'white' for i in range(len(agg658))]
    linewidths  = [2.5      if i == sresp_idx else 1.5     for i in range(len(agg658))]
    ax.barh(y_pos, agg658['n'], color=colors, edgecolor=edge_colors,
            linewidth=linewidths, height=0.72)
    ax.set_yticks(y_pos)
    ax.set_yticklabels([textwrap.fill(s, width=20) for s in agg658['cat']], fontsize=7)
    for i, label in enumerate(ax.get_yticklabels()):
        if i == sresp_idx:
            label.set_color(NPG['red'])
            label.set_weight('bold')
    ax.set_xlabel('Number of reactions GAINED (vs other strains)')
    ax.set_title('B', loc='left', x=-0.45, fontsize=11, fontweight='bold')
    for i, (n, p) in enumerate(zip(agg658['n'], agg658['min_p'])):
        bold = (i == sresp_idx)
        color = NPG['red'] if bold else SEM['text']
        ax.text(n + 0.15, i, f'{n}  ({fmt_p(p)})', va='center', fontsize=8,
                color=color, fontweight='bold' if bold else 'normal')
    ax.set_xlim(0, agg658['n'].max() * 1.55 if len(agg658) else 1)

    if sresp_idx is not None:
        n = agg658['n'].iloc[sresp_idx]
        xmax = agg658['n'].max() * 1.55
        ax.annotate('Direct mechanism for\nin vitro sulfite (z = +0.95)\n— central finding',
                    xy=(n + 0.15, sresp_idx),
                    xytext=(xmax * 0.62, sresp_idx + 1.2),
                    fontsize=8.5, color=NPG['red'], fontweight='bold',
                    ha='center', va='center',
                    bbox=dict(boxstyle='round,pad=0.4', facecolor='#FFF4F0',
                              edgecolor=NPG['red'], linewidth=1.2),
                    arrowprops=dict(arrowstyle='->', color=NPG['red'], lw=1.5,
                                    connectionstyle='arc3,rad=-0.2'))

    save_fig(fig, 'Figure3_mechanism')
    plt.close(fig)
    print('Done.')


if __name__ == '__main__':
    main()
