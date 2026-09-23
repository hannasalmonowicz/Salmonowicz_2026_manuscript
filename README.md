# Code availability

Code for the figures in:

Salmonowicz H. et al., Therapy-induced senescence preserves dormant OXPHOS to promote
metabolic adaptability (manuscript).

This repository contains the code used to generate the figures listed below
(Fig 1E, 1G, 1H, 2B, 2G/H, 4C–G, 4J, 5C, 5E, 5G, 6C/D, 6E, EV1O, EV6G/H, EV7B, EV10F/G–H).
Each script reads its input files from `data/` and writes its output next to itself.

Each script's only required inputs are the manuscript's EV tables. A small number of
fixed gene lists are defined directly within the scripts that use them, and the GO
annotation required for Figure EV7B is included as a column in EV_Table_8 (see "GO
annotation" below).

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

**Statistics.** GSEA (Fig 4J, 5E) is pre-ranked and run with the `gseapy` package.
Fig 1H and EV1O assess pathway-category significance across all 30 models using a
resampling-based permutation test, independent of the EV tables' own per-protein
p-values. Fig 2B uses an ANOVA with Bonferroni correction. The metabolomics panels
(Fig 6E and the boxplot script) report the source facility's FDR/q-values without
recalculation, except the ratio panels, which use a Welch's t-test.

EV_Table_10, EV_Table_11, and EV_Table_12 are matched by filename pattern rather than
an exact name (e.g. any file starting with `EV_Table_10` and ending `.xlsx`), so the
scripts still work even if the final table numbers in the published manuscript end up
different from what's used here.

## Running

Place the manuscript's EV tables under `figures/data/`, then run any script from the
`figures/` directory, e.g.:

```bash
cd figures
python3 fig1G_generality_matrix.py
```

Each script checks `./data/` first, then the current folder, and raises a clear
`FileNotFoundError` naming exactly what's missing.

## GO annotation

EV_Table_8 contains the GO Biological Process annotation column (`GO_Names_P`) required
to reproduce Figure EV7B (`GOenrichment_CAM_vs_CTRL_bubbleplots.py`).

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

Hanna Salmonowicz. Code written with the assistance of Claude (Anthropic).

## License

MIT, see LICENSE.
