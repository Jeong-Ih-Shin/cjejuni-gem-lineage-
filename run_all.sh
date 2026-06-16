#!/usr/bin/env bash
# Run the full pipeline end-to-end: GEM SBMLs → presence/absence → enrichment → figures.
# Run from the repository root:  bash run_all.sh
set -e
cd "$(dirname "$0")/scripts"

echo "=== Step 01: build presence/absence matrix from GEM files ==="
python 01_build_presence_absence.py
echo

echo "=== Step 02: compute per-CC Fisher exact enrichment ==="
python 02_compute_enrichment.py
echo

echo "=== Step 03: Figure 1 (GEM lineage signal) ==="
python 03_figure1_gem_lineage.py
echo

echo "=== Step 04: Figure 2 (in vitro phenotype lineage signal) ==="
python 04_figure2_invitro_profiles.py
echo

echo "=== Step 05: Figure 3 (ST-443 vs ST-658 mechanism) ==="
python 05_figure3_mechanism.py
echo

echo "=== Step 06: in silico growth simulation on Western + High-fiber diets ==="
python 06_human_gut_simulation.py
echo

echo "=== Step 07: Figure 4 (ST-443 vs ST-658 comparison) ==="
python 07_figure4_st443_vs_st658.py
echo

echo "=== Done. Outputs are in output/figures/ ==="
