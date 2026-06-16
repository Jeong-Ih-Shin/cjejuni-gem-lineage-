"""
Extract EC numbers and KEGG reaction IDs from ALL 86 SBML files (union)
for ST-443 LOST and ST-658 GAINED differential reactions.

Uses regex parsing (fast, ~10 sec total) instead of cobrapy loading.

Outputs (in intermediate/):
  - all_reactions_annotations.csv    : reaction → EC, KEGG mapping from full union
  - ST443_diff_with_annotations.csv  : ST-443 LOST reactions with annotations
  - ST658_diff_with_annotations.csv  : ST-658 GAINED reactions with annotations
  - ST443_EC_list.txt / ST658_EC_list.txt    : EC numbers only
  - ST443_KEGG_list.txt / ST658_KEGG_list.txt : KEGG reaction IDs only

Usage (from repo root):
    python scripts/extract_annotations_for_DAVID.py
"""
import re
import pandas as pd
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
INT  = REPO / 'intermediate'
GEMS = REPO / 'gems'

# ---- 1. Scan all SBML files for reaction annotations -------------------
RXN_HEAD_RE = re.compile(r'<reaction\s+metaid="R_(rxn\d+_c0)"')
EC_RE       = re.compile(r'identifiers\.org/ec-code:([\d\.\-]+)"')
KEGG_RE     = re.compile(r'identifiers\.org/kegg\.reaction:(R\d+)"')

mapping = {}  # rxn_id → {'ec': set, 'kegg': set}

sbml_files = sorted(GEMS.glob('cmccj_*.xml'))
print(f'Scanning {len(sbml_files)} SBML files...')
for i, sbml in enumerate(sbml_files, 1):
    text = sbml.read_text()
    chunks = re.split(r'(?=<reaction\s+metaid="R_rxn)', text)
    for ch in chunks:
        m = RXN_HEAD_RE.match(ch)
        if not m:
            continue
        rxn_id = m.group(1)
        block = ch.split('</reaction>', 1)[0]
        ecs   = set(EC_RE.findall(block))
        keggs = set(KEGG_RE.findall(block))
        if rxn_id not in mapping:
            mapping[rxn_id] = {'ec': ecs, 'kegg': keggs}
        else:
            mapping[rxn_id]['ec']   |= ecs
            mapping[rxn_id]['kegg'] |= keggs
    if i % 20 == 0:
        print(f'  [{i:3d}/{len(sbml_files)}] processed')

print(f'  Total unique reactions: {len(mapping)}')
print(f'  With EC:   {sum(1 for v in mapping.values() if v["ec"])}')
print(f'  With KEGG: {sum(1 for v in mapping.values() if v["kegg"])}')

map_df = pd.DataFrame([
    {'reaction': k,
     'ec':   ';'.join(sorted(v['ec'])),
     'kegg': ';'.join(sorted(v['kegg']))}
    for k, v in mapping.items()
])
map_df.to_csv(INT / 'all_reactions_annotations.csv', index=False)


# ---- 2. Process one CC's differential reactions ------------------------
def process(label, csv_in, direction, csv_out, ec_out, kegg_out):
    df = pd.read_csv(csv_in)
    sig = df[(df['pct_diff'].apply(direction)) & (df['p'] < 0.01)].copy()
    merged = sig.merge(map_df, on='reaction', how='left').fillna('')
    merged.to_csv(csv_out, index=False)

    ec_list = (merged['ec'].astype(str)
               .replace('', pd.NA).dropna()
               .str.split(';').explode().unique())
    with open(ec_out, 'w') as f:
        f.write('\n'.join(sorted(ec_list)) + '\n')

    kegg_list = (merged['kegg'].astype(str)
                 .replace('', pd.NA).dropna()
                 .str.split(';').explode().unique())
    with open(kegg_out, 'w') as f:
        f.write('\n'.join(sorted(kegg_list)) + '\n')

    print(f'\n{label}: {len(merged)} differential reactions')
    print(f'  with EC:   {(merged["ec"]   != "").sum():3d}')
    print(f'  with KEGG: {(merged["kegg"] != "").sum():3d}')
    print(f'  unique EC:   {len(ec_list):3d}  →  {ec_out.name}')
    print(f'  unique KEGG: {len(kegg_list):3d}  →  {kegg_out.name}')


# ---- 3. ST-443 LOST ---------------------------------------------------
process(
    label='ST-443 LOST',
    csv_in=INT / 'ST443_vs_others_enrichment.csv',
    direction=lambda x: x < -50,
    csv_out=INT / 'ST443_diff_with_annotations.csv',
    ec_out=INT / 'ST443_EC_list.txt',
    kegg_out=INT / 'ST443_KEGG_list.txt',
)

# ---- 4. ST-658 GAINED -------------------------------------------------
process(
    label='ST-658 GAINED',
    csv_in=INT / 'ST658_vs_others_enrichment.csv',
    direction=lambda x: x > 50,
    csv_out=INT / 'ST658_diff_with_annotations.csv',
    ec_out=INT / 'ST658_EC_list.txt',
    kegg_out=INT / 'ST658_KEGG_list.txt',
)

print('\nDone.')
print('\nNext: paste the EC list into:')
print('  • KEGG Mapper: https://www.kegg.jp/kegg/mapper/search.html')
