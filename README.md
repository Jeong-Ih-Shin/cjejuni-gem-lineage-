# *C. jejuni* GEM × in vitro phenotype analysis — reproducibility repository

Code and minimal data for the manuscript:

> **Lineage-organized metabolic configuration and lineage-restricted GEM–phenotype mapping in 86 Korean *Campylobacter jejuni* clinical isolates.**
> (Manuscript in preparation)

The full analysis pipeline runs end-to-end **from gapseq-reconstructed GEM SBML files** to the four main-text figures.

---

## Repository layout

```
.
├── README.md
├── requirements.txt
├── run_all.sh                ← runs the full pipeline (steps 01–06)
│
├── data/                     ← raw inputs only (user-supplied)
│   └── metadata.xlsx     in vitro phenotype dataset
│
├── gems/                     ← GEM SBML files (user-supplied; one per strain)
│   └── README.md             instructions for populating this directory
│
├── scripts/                  ← analysis pipeline (numbered, run in order)
│   ├── utils.py              shared paths, palette, plot style, ID mapping
│   ├── 01_build_presence_absence.py   gems/*.xml      → intermediate/reaction_presence_absence.csv
│   ├── 02_compute_enrichment.py       presence/absence + CC labels → enrichment CSVs
│   ├── 03_figure1_gem_lineage.py      Jaccard, PCA, per-CC bars
│   ├── 04_figure2_invitro_profiles.py 24h fc z-score, heatmap, per-CC bars
│   ├── 05_figure3_mechanism.py        ST-443 LOST + ST-658 GAINED reactions, manual pathway grouping
│   ├── 06_human_gut_simulation.py     FBA on Western + High-fiber diets → mu per strain
│   ├── 07_figure4_st443_vs_st658.py   one-sided MW + Cohen d, Figure 4
│   └── validate_with_KEGG.py        optional: KEGG pathway/module validation of Figure 3 groupings
│
├── intermediate/             ← created at runtime
│   ├── reaction_presence_absence.csv
│   ├── reaction_names.csv
│   ├── ST443_vs_others_enrichment.csv
│   ├── ST658_vs_others_enrichment.csv
│   ├── human_gut_simulation_results.csv
│   └── st443_vs_st658_statistics.csv
│
└── output/                   ← created at runtime
    ├── figures/              Figure1–Figure4 (PNG + SVG, 300 dpi)
```

---

## Quick start

```bash
# 1. Install Python dependencies
pip install -r requirements.txt

# 2. Populate the inputs:
#    - data/metadata.xlsx     (provided)
#    - gems/*.xml                                    (provide your own SBML files; cmccj_xxxx.xml form)

# 3. Run the full pipeline
bash run_all.sh
```

That's it. Outputs land in `output/figures/`. To re-run only a single step, execute the corresponding `scripts/0N_*.py` directly.

---

## Inputs

### `data/metadata.xlsx`
In-house OD600 substrate-utilization assay across 86 strains × 7 substrates × 7 timepoints. Sheets:
- `Long_format` — one row per (sample, chemical, timepoint), value = trimmed fold change relative to inoculum (`fc_trimmed`)
- `Cluster_assignments` — sample → CC (Clonal Complex) mapping

The single per-(strain, chemical) metric used in all downstream analyses is the **24-hour fold change** (24h fc), z-score normalized per chemical.

### `gems/*.xml` (or `*.sbml`)
Gapseq-reconstructed strain-specific GEMs, one file per strain, named with the unified `cmccj_xxxx.xml` convention. Strain IDs in the GEM filenames must match the `sample` IDs used in `metadata.xlsx`. See `gems/README.md` for the upstream gapseq command.

---

## Pipeline steps

| Step | Script | Inputs | Outputs |
|---|---|---|---|
| 01 | `01_build_presence_absence.py` | `gems/*.xml`, `data/strain_id_map.xlsx` | `intermediate/reaction_presence_absence.csv`, `intermediate/reaction_names.csv` |
| 02 | `02_compute_enrichment.py` | presence/absence + `Cluster_assignments` | `intermediate/ST{443,658}_vs_others_enrichment.csv` |
| 03 | `03_figure1_gem_lineage.py` | presence/absence + CCs | `output/figures/Figure1_gem_lineage.{png,svg}` |
| 04 | `04_figure2_invitro_profiles.py` | `metadata.xlsx` (Long_format + Cluster_assignments) | `output/figures/Figure2_invitro_profiles.{png,svg}` |
| 05 | `05_figure3_mechanism.py` | step-02 enrichment CSVs | `output/figures/Figure3_mechanism.{png,svg}` |
| 06 | `06_human_gut_simulation.py` | `gems/*.xml`, `data/metadata.xlsx` (Cluster_assignments) | `intermediate/human_gut_simulation_results.csv` |
| 07 | `07_figure4_st443_vs_st658.py` | step-06 results | `intermediate/st443_vs_st658_statistics.csv`, `output/figures/Figure4_st443_vs_st658.{png,svg}` |

`run_all.sh` executes them in order.

---

## Upstream pipeline (from raw sequences to GEMs)

This repository starts from reconstructed GEMs. The upstream steps were:

1. **Sequencing & assembly** — Illumina paired-end → SPAdes / Unicycler
2. **MLST / CC assignment** — `mlst` (T. Seemann) → PubMLST *C. jejuni / coli* scheme
3. **GEM reconstruction** — [gapseq v1.2 (commit 6c2ff0e9)](https://github.com/jotech/gapseq):

   ```bash
   gapseq doall <strain>.fna -m Bacteria -b 100
   ```

   The resulting `<strain>.xml` is placed in `gems/`.

---

## Notes on the pathway classification (Figure 3)

Step 05 classifies significantly differential reactions into functional pathway categories by **keyword matching on the SBML reaction name** (extracted in step 01). This is a heuristic. For any specific reaction-pathway claim, verify the assignment against ModelSEED / KEGG annotation.

The classification keywords are inline at the top of `05_figure3_mechanism.py` (`classify_443`, `classify_658`) and can be edited.

---

## Environment

Tested on Python 3.11, Ubuntu 24.04 (WSL2). gapseq was run separately on Linux (Ubuntu 22.04).

---

## License

MIT.

---

## Citation

Citation will be added upon publication.

---

## Optional: External validation against KEGG

Pathway groupings used in Figure 3 (manual keyword matching on reaction names)
were independently validated against KEGG pathway and module annotations.

```bash
python scripts/validate_with_KEGG.py
```

This extracts EC numbers and KEGG reaction IDs from SBML annotations (union
scan across all 86 strain models). Outputs in `intermediate/` enable direct
cross-reference with KEGG Mapper (https://www.kegg.jp/kegg/mapper/). The major
Figure 3 categories — Pantothenate/CoA biosynthesis, GDP-mannose/fucose
biosynthesis, anaerobic sulfur respiration, sialic acid biosynthesis, and
fatty acid biosynthesis — were independently confirmed.

Requires internet access for KEGG REST tables (~5 MB, cached on first run).
