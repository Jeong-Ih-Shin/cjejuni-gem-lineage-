"""
Step 08 - Revision analyses (response to reviewer comments).

Adds the analyses requested during peer review:

  1. Strain-level assembly and model quality metrics
     (CheckM2 completeness/contamination, QUAST metrics, cgMLST called alleles,
      Roary gene counts, gapseq gap-filling) and their association with
      reaction counts.
  2. Partial Mantel test: reaction-content distance vs CC membership,
     controlling for a quality-metric distance matrix.
  3. Sensitivity analysis excluding the one contaminated isolate.
  4. Isolate-level bootstrap confidence intervals for within-/between-CC
     distances (GSMM and in vitro).
  5. Direct CC-443 vs CC-658 comparison, all isolates and chicken only.
  6. Reaction-level nominal and FDR-adjusted enrichment for all CCs.
  7. Two-sided test of simulated growth with bootstrap effect-size interval.

Inputs:
    data/metadata.xlsx                              (Cluster_assignments, Long_format)
    intermediate/reaction_presence_absence.csv      (from step 01)
    intermediate/human_gut_simulation_results.csv   (from step 07)
    intermediate/quality_report.tsv                 (CheckM2 output)
    intermediate/gapfill_counts.csv                 (gapseq draft vs final reactions)
    intermediate/assembly_qc.csv                    (QUAST + cgMLST + gene counts)

Outputs (all written to intermediate/):
    revision_strain_qc.csv
    revision_qc_correlations.csv
    revision_partial_mantel.csv
    revision_bootstrap_ci_gsmm.csv
    revision_bootstrap_ci_invitro.csv
    revision_direct_cc443_vs_cc658.csv
    revision_reaction_enrichment_fdr.csv
    revision_allcc_enrichment.csv
    revision_growth_twosided.csv
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform, cdist
from itertools import combinations

from scipy.stats import spearmanr, rankdata, fisher_exact

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / 'data'
INTERMEDIATE_DIR = REPO_ROOT / 'intermediate'

CHEMICALS = ['Lactate', 'Serine', 'Proline', 'Asparagine',
             'Formate', 'Sulfite', 'Gluconate']
CONTAMINATED = 'cmccj_1273'          # CheckM2 contamination 21.8%
QC_METRICS = ['Completeness_pct', 'Contamination_pct', 'Contigs', 'N50',
              'Genome_size_bp', 'GC_pct', 'cgMLST_called_pct',
              'Gene_count', 'Gapfilled_reactions']
N_BOOT = 2000
N_PERM = 9999
SEED = 42


# ============================================================
# Helpers
# ============================================================
def cc_label(cc):
    """'ST-21 complex' -> 'CC-21' (display label only)."""
    return str(cc).replace(' complex', '').replace('ST-', 'CC-')


def bh_fdr(p):
    """Benjamini-Hochberg FDR across all tested reactions."""
    p = np.asarray(p, float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1)
    return out


def load_reactions():
    pa = pd.read_csv(INTERMEDIATE_DIR / 'reaction_presence_absence.csv',
                     index_col=0).T
    pa = pa[~pa.index.str.contains('KCTC|ATCC', case=False, na=False)]
    return pa.astype(int)


def load_cc():
    meta = pd.read_excel(DATA_DIR / 'metadata.xlsx',
                         sheet_name='Cluster_assignments')
    meta = meta[~meta['sample'].astype(str)
                .str.contains('KCTC|ATCC', case=False, na=False)]
    return meta


def load_qc():
    """Merge assembly QC, CheckM2 and gap-filling into one strain-level table."""
    qc = pd.read_csv(INTERMEDIATE_DIR / 'assembly_qc.csv', index_col=0)
    checkm = pd.read_csv(INTERMEDIATE_DIR / 'quality_report.tsv', sep='\t')
    gap = pd.read_csv(INTERMEDIATE_DIR / 'gapfill_counts.csv', index_col=0)
    qc = qc.join(checkm.set_index('Strain')[['Completeness_pct', 'Contamination_pct']]) \
           .join(gap[['Draft_reactions', 'Gapfilled_reactions']])
    return qc


# ============================================================
# 1. Quality metrics vs reaction counts
# ============================================================
def qc_correlations(pa, qc, samples):
    rc = pa.loc[samples].sum(axis=1)
    rows = []
    for m in QC_METRICS:
        v = pd.to_numeric(qc.loc[samples, m], errors='coerce')
        keep = v.notna()
        r, p = spearmanr(v[keep], rc[keep])
        rows.append({'metric': m, 'spearman_rho': r, 'p_value': p,
                     'n': int(keep.sum())})
    out = pd.DataFrame(rows)
    out.to_csv(INTERMEDIATE_DIR / 'revision_qc_correlations.csv', index=False)
    print('\n[1] Quality metrics vs reaction count')
    print(out.to_string(index=False))
    return out


# ============================================================
# 2-3. Partial Mantel, with and without the contaminated isolate
# ============================================================
def _partial_corr(a, b, c):
    rab = np.corrcoef(a, b)[0, 1]
    rac = np.corrcoef(a, c)[0, 1]
    rbc = np.corrcoef(b, c)[0, 1]
    return (rab - rac * rbc) / np.sqrt((1 - rac ** 2) * (1 - rbc ** 2))


def partial_mantel(pa, qc, cc_map, samples, label, n_perm=N_PERM, seed=SEED):
    """Reaction distance ~ CC membership, controlling for quality distance."""
    sub = pa.loc[samples]
    n = len(samples)
    ccv = np.array([str(cc_map[s]) for s in samples])

    Q = qc.loc[samples, QC_METRICS].apply(pd.to_numeric, errors='coerce')
    Q = ((Q - Q.mean()) / Q.std()).fillna(0)

    A = squareform(pdist(sub.values, metric='jaccard'))   # reaction content
    B = (ccv[:, None] != ccv[None, :]).astype(float)      # different CC = 1
    C = squareform(pdist(Q.values))                       # quality distance
    iu = np.triu_indices(n, 1)

    a, b, c = rankdata(A[iu]), rankdata(B[iu]), rankdata(C[iu])
    r_cc = _partial_corr(a, b, c)     # CC effect, quality controlled
    r_qc = _partial_corr(a, c, b)     # quality effect, CC controlled

    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_perm):
        perm = rng.permutation(n)
        a_p = rankdata(A[np.ix_(perm, perm)][iu])
        if abs(_partial_corr(a_p, b, c)) >= abs(r_cc):
            hits += 1
    p_cc = (hits + 1) / (n_perm + 1)

    # how many CCs keep within < between
    ser = pd.Series(ccv, index=samples)
    counts = ser.value_counts()
    coherent = 0
    total = 0
    for this_cc in counts[counts >= 2].index:
        ins = ser[ser == this_cc].index
        outs = ser[ser != this_cc].index
        within = pdist(sub.loc[ins].values, metric='jaccard').mean()
        between = cdist(sub.loc[ins].values, sub.loc[outs].values,
                        metric='jaccard').mean()
        total += 1
        coherent += int(within < between)

    return {'analysis': label, 'n_isolates': n,
            'partial_r_CC_given_quality': r_cc, 'p_value': p_cc,
            'partial_r_quality_given_CC': r_qc,
            'CCs_with_within_lt_between': f'{coherent}/{total}'}


# ============================================================
# 4. Isolate-level bootstrap CIs
# ============================================================
def bootstrap_ci(in_mat, out_mat, metric='jaccard', n_boot=N_BOOT, seed=SEED):
    """Resample isolates (not pairs), because pairwise distances are not independent."""
    rng = np.random.default_rng(seed)
    ws, bs = [], []
    n_in, n_out = len(in_mat), len(out_mat)
    for _ in range(n_boot):
        i_idx = rng.integers(0, n_in, n_in)
        o_idx = rng.integers(0, n_out, n_out)
        if len(np.unique(i_idx)) < 2:
            continue
        bi, bo = in_mat[i_idx], out_mat[o_idx]
        ws.append(pdist(bi, metric=metric).mean())
        bs.append(cdist(bi, bo, metric=metric).mean())
    if not ws:
        return (np.nan,) * 4
    return (np.percentile(ws, 2.5), np.percentile(ws, 97.5),
            np.percentile(bs, 2.5), np.percentile(bs, 97.5))


def per_cc_bootstrap(mat_df, cc_series, metric, out_name):
    counts = cc_series.value_counts()
    rows = []
    for this_cc in counts[counts >= 2].index:
        ins = cc_series[cc_series == this_cc].index
        outs = cc_series[cc_series != this_cc].index
        A = mat_df.loc[ins].values
        B = mat_df.loc[outs].values
        within = pdist(A, metric=metric).mean()
        between = cdist(A, B, metric=metric).mean()
        w_lo, w_hi, b_lo, b_hi = bootstrap_ci(A, B, metric=metric)
        rows.append({'CC': cc_label(this_cc), 'n': int(counts[this_cc]),
                     'within': within, 'within_lo': w_lo, 'within_hi': w_hi,
                     'between': between, 'between_lo': b_lo, 'between_hi': b_hi,
                     'within_lt_between': within < between})
    out = pd.DataFrame(rows).sort_values('n', ascending=False)
    out.to_csv(INTERMEDIATE_DIR / out_name, index=False)
    return out


# ============================================================
# 5-6. Enrichment: direct comparison, nominal vs FDR, all CCs
# ============================================================
def enrichment(pa, in_ids, out_ids, min_pct_diff=50):
    """Fisher's exact test per reaction between two isolate groups."""
    n_in, n_out = len(in_ids), len(out_ids)
    rows = []
    for col in pa.columns:
        a = int(pa.loc[in_ids, col].sum())
        b = int(pa.loc[out_ids, col].sum())
        if a == 0 and b == 0:
            continue
        pct_in, pct_out = a / n_in * 100, b / n_out * 100
        diff = pct_in - pct_out
        _, p = fisher_exact([[a, n_in - a], [b, n_out - b]])
        rows.append({'reaction': col, 'pct_in': pct_in, 'pct_out': pct_out,
                     'pct_diff': diff, 'p_nominal': p})
    df = pd.DataFrame(rows)
    # Benjamini-Hochberg correction is applied across every reaction tested,
    # before the effect-size filter, so that the reported FDR reflects the full
    # search rather than the pre-selected subset.
    df['fdr_bh_all_tested'] = bh_fdr(df.p_nominal.values)
    df.attrs['n_tested'] = len(df)
    return df[abs(df.pct_diff) >= min_pct_diff].copy()


def pathway_summary(pa, cc_series, repo_root):
    """Representative KEGG pathways among the FDR-significant reactions of each CC."""
    from collections import Counter
    ann = pd.read_csv(repo_root / 'intermediate' / 'all_reactions_annotations.csv')
    r2p = {}
    with open(repo_root / 'intermediate' / 'kegg_reaction_pathway.tsv') as fh:
        for line in fh:
            a, b = line.strip().split('\t')
            r2p.setdefault(a.replace('rn:', ''), []).append(b.replace('path:', ''))
    with open(repo_root / 'intermediate' / 'kegg_pathway_names.tsv') as fh:
        pname = dict(l.strip().split('\t') for l in fh)
    # very broad maps carry little functional meaning here
    generic = {'map01100', 'map01110', 'map01120', 'map01200',
               'map01230', 'map01250', 'map02010', 'map01240'}
    rxn2paths = {}
    for _, r in ann.dropna(subset=['kegg']).iterrows():
        ps = [p for k in str(r['kegg']).split(';')
              for p in r2p.get(k.strip(), []) if p not in generic]
        if ps:
            rxn2paths[r['reaction']] = ps

    counts = cc_series.value_counts()
    rows = []
    for this_cc in counts[counts >= 2].index:
        ins = cc_series[cc_series == this_cc].index
        outs = cc_series[cc_series != this_cc].index
        e = enrichment(pa, ins, outs)
        sig = e[(e.p_nominal < 0.01) & (e.fdr_bh_all_tested < 0.05)]
        for direction, part in [('depleted', sig[sig.pct_diff < 0]),
                                ('enriched', sig[sig.pct_diff > 0])]:
            cnt = Counter(p for r in part.reaction for p in rxn2paths.get(r, []))
            top = '; '.join(f'{pname.get(p, p)} ({n})' for p, n in cnt.most_common(3))
            rows.append({'CC': cc_label(this_cc), 'n': int(counts[this_cc]),
                         'direction': direction, 'n_reactions': len(part),
                         'top_pathways': top or '(no KEGG pathway annotation)'})
    out = pd.DataFrame(rows)
    out.to_csv(INTERMEDIATE_DIR / 'revision_cc_pathway_summary.csv', index=False)
    return out


def all_cc_summary(pa, cc_series):
    counts = cc_series.value_counts()
    rows = []
    for this_cc in counts[counts >= 2].index:
        ins = cc_series[cc_series == this_cc].index
        outs = cc_series[cc_series != this_cc].index
        e = enrichment(pa, ins, outs)
        nom = e[e.p_nominal < 0.01]
        fdr = nom[nom.fdr_bh_all_tested < 0.05]
        rows.append({'CC': cc_label(this_cc), 'n': int(counts[this_cc]),
                     'depleted_nominal': int((nom.pct_diff < 0).sum()),
                     'enriched_nominal': int((nom.pct_diff > 0).sum()),
                     'depleted_fdr': int((fdr.pct_diff < 0).sum()),
                     'enriched_fdr': int((fdr.pct_diff > 0).sum())})
    out = pd.DataFrame(rows).sort_values('n', ascending=False)
    out.to_csv(INTERMEDIATE_DIR / 'revision_allcc_enrichment.csv', index=False)
    return out


# ============================================================
# 7. Growth: two-sided test and bootstrap effect size
# ============================================================
def cohen_d(x, y):
    nx, ny = len(x), len(y)
    s = np.sqrt(((nx - 1) * np.var(x, ddof=1) +
                 (ny - 1) * np.var(y, ddof=1)) / (nx + ny - 2))
    return (np.mean(x) - np.mean(y)) / s if s > 1e-12 else np.nan


def exact_permutation_p(x, y, decimals=10):
    """Two-sided exact permutation test on rank sums.

    FBA predictions differ only at solver precision (about 1e-15), so values are
    rounded before ranking; tied values receive the average rank. All possible
    group assignments are enumerated, which is feasible at these sample sizes.
    """
    xr = np.round(x, decimals)
    yr = np.round(y, decimals)
    pooled = np.concatenate([xr, yr])
    n, nx = len(pooled), len(xr)
    ranks = rankdata(pooled)
    obs_u = ranks[:nx].sum() - nx * (nx + 1) / 2
    expected = nx * len(yr) / 2
    us = np.array([ranks[list(idx)].sum() - nx * (nx + 1) / 2
                   for idx in combinations(range(n), nx)])
    return float(np.mean(np.abs(us - expected) >= abs(obs_u - expected) - 1e-12)), len(us)


def growth_tests(sim, cc443, cc658, cc658_chicken):
    rows = []
    for col, label in [('mu_western_diet', 'Western diet'),
                       ('mu_high_fiber_diet', 'High-fiber diet')]:
        for g658, gl in [(cc658, 'all CC-658'), (cc658_chicken, 'CC-658 chicken only')]:
            x = sim[sim.strain.isin(cc443)][col].dropna().values
            y = sim[sim.strain.isin(g658)][col].dropna().values
            p_two, n_perm = exact_permutation_p(x, y)
            d = cohen_d(x, y)
            rng = np.random.default_rng(SEED)
            ds = []
            for _ in range(10000):
                di = cohen_d(rng.choice(x, len(x), replace=True),
                             rng.choice(y, len(y), replace=True))
                if np.isfinite(di):
                    ds.append(di)
            lo, hi = (np.percentile(ds, 2.5), np.percentile(ds, 97.5)) if ds else (np.nan, np.nan)
            rows.append({'diet': label, 'comparison': gl,
                         'n_CC443': len(x), 'n_CC658': len(y),
                         'mean_CC443': x.mean(), 'sd_CC443': x.std(ddof=1),
                         'mean_CC658': y.mean(),
                         'sd_CC658': y.std(ddof=1) if len(y) > 1 else 0.0,
                         'p_exact_two_sided': p_two, 'n_permutations': n_perm, 'cohens_d': d,
                         'd_boot_lo': lo, 'd_boot_hi': hi,
                         'boot_valid_frac': len(ds) / 10000,
                         'boot_excluded_zero_sd_frac': 1 - len(ds) / 10000})
    out = pd.DataFrame(rows)
    out.to_csv(INTERMEDIATE_DIR / 'revision_growth_twosided.csv', index=False)
    return out


# ============================================================
# Main
# ============================================================
def main():
    pa = load_reactions()
    meta = load_cc()
    cc_map = dict(zip(meta['sample'], meta['CC']))
    qc = load_qc()

    samples = [s for s in pa.index if s in cc_map and s in qc.index]
    cc_series = pd.Series([cc_map[s] for s in samples], index=samples)
    print(f'Isolates analysed: {len(samples)}')

    qc.loc[samples].to_csv(INTERMEDIATE_DIR / 'revision_strain_qc.csv')

    # 1. quality metrics vs reaction counts
    qc_correlations(pa, qc, samples)

    # 2-3. partial Mantel, full cohort and excluding the contaminated isolate
    pm_rows = [partial_mantel(pa, qc, cc_map, samples, 'all isolates')]
    keep = [s for s in samples if s != CONTAMINATED]
    pm_rows.append(partial_mantel(pa, qc, cc_map, keep,
                                  f'excluding {CONTAMINATED}'))
    pm = pd.DataFrame(pm_rows)
    pm.to_csv(INTERMEDIATE_DIR / 'revision_partial_mantel.csv', index=False)
    print('\n[2-3] Partial Mantel (quality controlled)')
    print(pm.to_string(index=False))

    # 4. bootstrap CIs - GSMM
    gsmm_ci = per_cc_bootstrap(pa.loc[samples], cc_series, 'jaccard',
                               'revision_bootstrap_ci_gsmm.csv')
    print('\n[4] Bootstrap CI (GSMM Jaccard)')
    print(gsmm_ci.to_string(index=False))

    # 4. bootstrap CIs - in vitro phenotypes
    phen = pd.read_excel(DATA_DIR / 'metadata.xlsx', sheet_name='Long_format')
    phen['fc_trimmed'] = pd.to_numeric(phen['fc_trimmed'], errors='coerce')
    phen['timepoint'] = pd.to_numeric(phen['timepoint'], errors='coerce')
    phen = phen.dropna(subset=['fc_trimmed', 'timepoint'])
    phen = phen[~phen['sample'].astype(str)
                .str.contains('KCTC|ATCC', case=False, na=False)]
    wide = (phen[phen.timepoint == 24]
            .pivot_table(index='sample', columns='chemical', values='fc_trimmed')
            .reindex(columns=CHEMICALS).dropna())
    # population SD (ddof=0), matching StandardScaler used in the figure scripts
    z = (wide - wide.mean()) / wide.std(ddof=0)
    z_samples = [s for s in z.index if s in cc_map]
    invitro_ci = per_cc_bootstrap(z.loc[z_samples],
                                  pd.Series([cc_map[s] for s in z_samples],
                                            index=z_samples),
                                  'euclidean',
                                  'revision_bootstrap_ci_invitro.csv')
    print('\n[4] Bootstrap CI (in vitro Euclidean)')
    print(invitro_ci.to_string(index=False))

    # 5. direct CC-443 vs CC-658 comparison
    cc443 = [s for s in samples if str(cc_map[s]) == 'ST-443 complex']
    cc658 = [s for s in samples if str(cc_map[s]) == 'ST-658 complex']
    cc658_chicken = [s for s in cc658
                     if str(meta.set_index('sample').loc[s, 'Origin']).lower() == 'chicken']
    direct = enrichment(pa, cc443, cc658)
    direct = direct.rename(columns={'pct_in': 'pct_CC443', 'pct_out': 'pct_CC658'})
    chicken = enrichment(pa, cc443, cc658_chicken)
    direct = direct.merge(
        chicken[['reaction', 'p_nominal']].rename(
            columns={'p_nominal': 'p_chicken_only'}),
        on='reaction', how='left')
    direct.sort_values('p_nominal').to_csv(
        INTERMEDIATE_DIR / 'revision_direct_cc443_vs_cc658.csv', index=False)
    print(f'\n[5] Direct CC-443 vs CC-658: {len(direct)} reactions with |diff| >= 50 points')
    print(f'    CC-658 chicken isolates: {len(cc658_chicken)}')

    # 6. reaction-level nominal vs FDR for the two focal lineages
    parts = []
    for ids, name in [(cc443, 'CC-443'), (cc658, 'CC-658')]:
        others = [s for s in samples if s not in ids]
        e = enrichment(pa, ids, others)
        e.insert(0, 'focal_CC', name)
        parts.append(e[e.p_nominal < 0.01])
    focal = pd.concat(parts)
    focal.to_csv(INTERMEDIATE_DIR / 'revision_reaction_enrichment_fdr.csv', index=False)
    print('\n[6] Focal-lineage reactions (nominal P < 0.01)')
    print(focal.groupby(['focal_CC', focal.pct_diff > 0])
          .agg(n=('reaction', 'size'), n_fdr=('fdr_bh_all_tested', lambda s: int((s < 0.05).sum())))
          .to_string())

    allcc = all_cc_summary(pa, cc_series)
    print('\n[6] All CCs')
    print(allcc.to_string(index=False))

    paths = pathway_summary(pa, cc_series, REPO_ROOT)
    print('\n[6] Representative KEGG pathways per CC')
    print(paths.to_string(index=False))

    # 7. growth, two-sided
    sim = pd.read_csv(INTERMEDIATE_DIR / 'human_gut_simulation_results.csv')
    growth = growth_tests(sim, cc443, cc658, cc658_chicken)
    print('\n[7] Simulated growth (two-sided)')
    print(growth.to_string(index=False))

    print('\nDone. Outputs written to intermediate/.')


if __name__ == '__main__':
    main()
