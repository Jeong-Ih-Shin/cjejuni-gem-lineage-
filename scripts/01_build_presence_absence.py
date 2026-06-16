"""
Step 01 — Build reaction presence/absence matrix from GEM SBML files.

This is the entry point of the pipeline. Run this first to produce the
intermediate matrix that downstream steps depend on.

Inputs:
    gems/*.xml (or *.sbml)
        Gapseq-reconstructed strain-specific GEMs.
        Filename = strain ID (recommended: unified cmccj_xxxx.xml form).
        The optional gapseq "_target_jejuni" suffix is stripped automatically.

Outputs:
    intermediate/reaction_presence_absence.csv
        Pan-reactome × strains, binary 0/1.
        Rows: reaction IDs (union across all GEMs)
        Cols: strain IDs (from filename, with suffix stripped)

    intermediate/reaction_names.csv
        Two-column lookup: reaction, name (extracted from SBML <reaction name="...">).
        Used downstream for human-readable enrichment output.

Dependency:
    pip install cobra
"""
from pathlib import Path
import sys
import pandas as pd

# Allow direct execution (python scripts/01_...py) by ensuring scripts/ is on path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from utils import GEMS_DIR, INTERMEDIATE_DIR

try:
    import cobra
except ImportError:
    sys.exit('cobra not installed. Run: pip install cobra')


def extract_gem(sbml_path):
    """Return (reaction_id → name) dict for one GEM."""
    model = cobra.io.read_sbml_model(str(sbml_path))
    return {r.id: (r.name or r.id) for r in model.reactions}


def main():
    gem_files = sorted(list(GEMS_DIR.glob('*.xml')) + list(GEMS_DIR.glob('*.sbml')))
    if not gem_files:
        sys.exit(f'No GEM files (.xml or .sbml) found in {GEMS_DIR}. '
                 f'Place gapseq-reconstructed GEMs there before running.')

    print(f'Found {len(gem_files)} GEM files in {GEMS_DIR}')

    strain_reactions = {}     # strain_id → set of reaction IDs
    reaction_names = {}       # rxn_id → human-readable name (first encountered)
    failed = []

    for i, gem_path in enumerate(gem_files, 1):
        # Strip gapseq suffix "_target_jejuni" if present
        strain_id = gem_path.stem
        if strain_id.endswith('_target_jejuni'):
            strain_id = strain_id[:-len('_target_jejuni')]
        try:
            rxn_to_name = extract_gem(gem_path)
            strain_reactions[strain_id] = set(rxn_to_name.keys())
            for rid, rname in rxn_to_name.items():
                if rid not in reaction_names:
                    reaction_names[rid] = rname
            print(f'  [{i:>3}/{len(gem_files)}] {strain_id}: {len(rxn_to_name)} reactions')
        except Exception as e:
            failed.append((strain_id, str(e)))
            print(f'  [{i:>3}/{len(gem_files)}] {strain_id}: FAILED ({e})')

    if failed:
        print(f'\n{len(failed)} files failed to load. They will be excluded.')

    # If duplicates after canonicalization (e.g., both FG_0568.xml and MFDS2010853.xml),
    # keep only the first encountered — they should be the same strain.
    pan_reactions = sorted(set().union(*strain_reactions.values()))
    strains = sorted(strain_reactions.keys())
    print(f'\nPan-reactome size: {len(pan_reactions)} unique reactions')
    print(f'Strains in matrix: {len(strains)}')

    pa = pd.DataFrame(0, index=pan_reactions, columns=strains, dtype='int8')
    for strain, rxns in strain_reactions.items():
        pa.loc[list(rxns & set(pan_reactions)), strain] = 1

    out_matrix = INTERMEDIATE_DIR / 'reaction_presence_absence.csv'
    pa.to_csv(out_matrix)
    print(f'Saved: {out_matrix.relative_to(INTERMEDIATE_DIR.parent)}  (shape: {pa.shape})')

    out_names = INTERMEDIATE_DIR / 'reaction_names.csv'
    pd.DataFrame([{'reaction': r, 'name': n} for r, n in sorted(reaction_names.items())]).to_csv(out_names, index=False)
    print(f'Saved: {out_names.relative_to(INTERMEDIATE_DIR.parent)}  ({len(reaction_names)} entries)')


if __name__ == '__main__':
    main()
