# Code availability

Code for the figures in:

Salmonowicz H. et al., Therapy-induced senescence preserves dormant OXPHOS to promote
metabolic adaptability (manuscript).

This repository reads the manuscript's EV tables and produces each published figure:
pre-ranked GSEA on the proteomics data (Fig 4J, 5E), a cross-model permutation test for
pathway generality across 30 published senescence datasets (Fig 1G/1H, EV1O), GO
Biological Process enrichment (EV7B), an ANOVA + Bonferroni comparison on the
isolated-mitochondria proteomics (Fig 2B), and metabolomics quantification/statistics
plotting (Fig 6E and the boxplot panels). Each script reads its inputs from `data/` and
writes its output next to itself.

Every script's only inputs are the manuscript's own EV tables — nothing else to download,
no separate bundled files. Every non-EV-table dependency this repository used to have is
now either part of an EV table or hardcoded directly as a Python set literal inside the
script that uses it: five gene lists are hardcoded, the Payea et al. 2024 dataset is now
EV_Table_12, and the GO Biological Process annotation used by the EV7B bubble plots is now
its own column, `GO_Names_P`, inside EV_Table_8 itself (see "GO annotation" below).

## Figures

| Script | Figure | Output | EV table(s) needed |
|---|---|---|---|
| `fig1G_generality_matrix.py` | Fig 1G | 30-model generality dot-matrix | EV_Table_8, EV_Table_10, EV_Table_12 |
| `fig1H_generality_summary.py` | Fig 1H | Generality summary (diverging stacked bars) | EV_Table_8, EV_Table_10, EV_Table_12 |
| `fig2B_isolated_mito_pathway_ANOVA.py` | Fig 2B | Isolated-mito pathway vs. mitoproteome (ANOVA + Bonferroni) | EV_Table_2 |
| `fig4J_GSEA_SenMayo.py` | Fig 4J | SenMayo GSEA, GAL contrast | EV_Table_5 |
| `fig5E_GSEA_ATF4_ISR.py` | Fig 5E | GSEA, ATF4/ISR gene list | EV_Table_8 |
| `fig5G_ISR_heatmap.py` | Fig 5G | ISR/ATF4 heatmap, 10 genes | EV_Table_8, EV_Table_5 |
| `fig6E_pathway_scheme.py` | Fig 6E | CAM-track pathway scheme | EV_Table_11 |
| `figEV1O_OXPHOS_correlation.py` | EV1O | OXPHOS-vs-pathways correlation (3-panel scatter) | EV_Table_8, EV_Table_10, EV_Table_12 |
| `metabolite_boxplots_ALL_PANELS.py` | Fig 2G/H, 4C–G, 5C, 6C/D, EV6G/H, EV10F/G–H | 21 metabolite boxplot panels | EV_Table_3 (EV_Table_11 optional — brackets read "FDR pending" without it) |
| `trajectory_wholecell_proteomics.py` | Fig 1E | Time-resolved whole-cell trajectory plot | EV_Table_1 |
| `GOenrichment_CAM_vs_CTRL_bubbleplots.py` | EV7B | GO Biological Process enrichment, 3 bubble plots | EV_Table_8 (needs its `GO_Names_P` column) |

`fig2B`'s output reproduces the manuscript's reported n=4 (1C metabolism) and n=6
(Sulfur metabolism) exactly.

**Method in brief.** GSEA (Fig 4J, Fig 5E) is pre-ranked, run with `gseapy`, ranking
metric sign(log2FC) x -log10(p-value). Category-level significance across models
(Fig 1H, EV1O) is assessed by a competitive gene-resampling permutation test on Log2FC
values, not by the EV tables' own per-protein p/q-values. Fig 2B uses an ordinary ANOVA
with Bonferroni correction. Metabolomics panels (Fig 6E and the boxplot script) report
the facility's own FDR/q-values as-is; nothing is recomputed except the ratio panels,
which use a Welch's t-test on the log2 ratio.

EV_Table_10 (30-model meta-analysis, Anerillas et al.), EV_Table_11 (consolidated
metabolomics statistical results), and EV_Table_12 (Payea et al. 2024 dataset) are
matched by glob pattern (`EV_Table_10*.xlsx`,
`EV_Table_*_Metabolomics_Statistical_Results*.xlsx`, `EV_Table_12*.xlsx`), so any final
EV table number works without a code change.

## Running

Place the manuscript's EV tables under `figures/data/`, then run any script from the
`figures/` directory, e.g.:

```
cd figures
python3 fig1G_generality_matrix.py
```

Each script checks `./data/` first, then the current folder, and raises a clear
`FileNotFoundError` naming exactly what's missing.

## GO annotation (used by the EV7B bubble plots)

`GOenrichment_CAM_vs_CTRL_bubbleplots.py` needs one extra column on EV_Table_8,
`GO_Names_P`: for each row, the semicolon-joined set of GO Biological Process term names
for that row's gene(s) (unioned across every constituent symbol when `Gene names` is a
semicolon-joined protein group, e.g. `UQCRFS1;UQCRFS1P1`). It was generated once from the
original, unmodified QuickGO Biological Process annotation export (GO release 2025-06-01,
exported 2025-06-29), keeping only human (`TAXON ID` 9606) Biological Process
(`GO ASPECT` `'P'`) rows — see `verification/PROVENANCE.md` for the exact derivation and
match-rate numbers. Nothing else needs to be downloaded or bundled for this figure: the
annotation ships as part of the EV table itself, the same as every other input in this
repository.

## verification/

Not figure-producing code, and nothing here is meant to be run — `PROVENANCE.md` is a
plain-text reference document, not a script. It records where every non-EV-table input
(the five gene lists hardcoded directly in the scripts, the GO_Names_P annotation column
now part of EV_Table_8, and the Payea dataset before it became EV_Table_12) came from and
how it was derived, for the audit trail.

## Data

EV tables and other input data were generated for this study and are available as
described in the Data Availability statement of the manuscript. They are not part of
this repository.

## Environment

Python >= 3.10 with pandas, numpy, matplotlib, seaborn, scipy, openpyxl and gseapy (see
`requirements.txt`).
Tested with Python 3.11.15, pandas 3.0.2, numpy 2.4.4, matplotlib 3.10.9, seaborn
0.13.2, scipy 1.17.1, openpyxl 3.1.5, gseapy 1.3.1 — and separately with gseapy 1.1.9,
the version the GSEA scripts were originally written against.

## References

Rath S, Sharma R, Gupta R, Ast T, Chan C, Durham TJ, et al. (2021) MitoCarta3.0: an
updated mitochondrial proteome now with sub-organelle localization and pathway
annotations. Nucleic Acids Res 49(D1):D1541–D1547. doi:10.1093/nar/gkaa1011

Saul D, Kosinsky RL, Atkinson EJ, Doolittle ML, Zhang X, LeBrasseur NK, et al. (2022) A
new gene set identifies senescent cells and predicts senescence-associated pathways
across tissues. Nat Commun 13:4827. doi:10.1038/s41467-022-32552-1

## Author

Hanna Salmonowicz

## License

MIT, see LICENSE.
