# Code availability

Code for the figures in:

Salmonowicz H. et al., Therapy-induced senescence preserves dormant OXPHOS to promote
metabolic adaptability (manuscript).

This repository covers the proteomics, metabolomics and pathway-enrichment figures:
reading the manuscript's EV tables, differential abundance and enrichment testing
(GSEA, permutation testing, ANOVA), and the plotting code for each figure. Each script
reads its inputs from `data/` and writes its output next to itself.

## Figures

| Script | Figure | Output | Also needs |
|---|---|---|---|
| `fig1G_generality_matrix.py` | Fig 1G | 30-model generality dot-matrix (2 panels) | `Supplementary_PanelE_Data_1.xlsx`, `analMito_Anerillas2026.xlsx`, `analMito_Payea2024.xlsx`, `5_Annotations.xlsx` |
| `fig1H_generality_summary.py` | Fig 1H | Generality summary (diverging stacked bars) | `Supplementary_PanelE_Data_1.xlsx`, `analMito_Anerillas2026.xlsx`, `analMito_Payea2024.xlsx` |
| `fig2B_isolated_mito_pathway_ANOVA.py` | Fig 2B | Isolated-mito pathway vs. mitoproteome (ANOVA + Bonferroni) | `MitoCarta_pathway_annotation.xlsx`, `mitoproteome_isolated_mito.txt` (built by `verification/build_isolated_mito_annotation_from_curated_lists.py`) |
| `fig4J_GSEA_SenMayo.py` | Fig 4J | SenMayo GSEA, GAL contrast | `SAUL_SEN_MAYO.v2025.1.Hs.tsv` (SenMayo gene set, Saul et al. 2022) |
| `fig5E_GSEA_ATF4_ISR.py` | Fig 5E | GSEA, ATF4/ISR gene list | — |
| `fig5G_ISR_heatmap.py` | Fig 5G | ISR/ATF4 heatmap, 10 genes | — |
| `fig6E_pathway_scheme.py` | Fig 6E | CAM-track pathway scheme | `EV_Table_X_Metabolomics_Statistical_Results.xlsx` |
| `figEV1O_OXPHOS_correlation.py` | EV1O | OXPHOS-vs-pathways correlation (3-panel scatter) | `Supplementary_PanelE_Data_1.xlsx`, `analMito_Anerillas2026.xlsx`, `analMito_Payea2024.xlsx` |
| `metabolite_boxplots_ALL_PANELS.py` | Fig 2G/H, 4C–G, 5C, 6C/D, EV6G/H, EV10F/G–H | 21 metabolite boxplot panels | — |
| `trajectory_wholecell_proteomics.py` | Fig 1E | Time-resolved whole-cell trajectory plot | `Figure_1_protein_lists_Whole_Cell.xlsx` |
| `GOenrichment_CAM_vs_CTRL_bubbleplots.py` | EV7B | GO Biological Process enrichment, 3 bubble plots | `GOBP_annotations_CAM_IMT.xlsx` (built by `verification/build_GOBP_annotations_CAM_IMT.py`) |

`fig2B`'s output reproduces the manuscript's reported n=4 (1C metabolism) and n=6
(Sulfur metabolism) exactly.

**Method in brief.** GSEA (Fig 4J, Fig 5E) is pre-ranked, run with `gseapy`, ranking
metric sign(log2FC) x -log10(p-value). Category-level significance across models
(Fig 1H, EV1O) is assessed by a competitive gene-resampling permutation test on Log2FC
values, not by the EV tables' own per-protein p/q-values. Fig 2B uses an ordinary ANOVA
with Bonferroni correction. Metabolomics panels (Fig 6E and the boxplot script) report
the facility's own FDR/q-values as-is; nothing is recomputed except the ratio panels,
which use a Welch's t-test on the log2 ratio.

## Running

Place the manuscript's EV tables and the other files listed above under `figures/data/`,
then run any script from the `figures/` directory, e.g.:

```
cd figures
python3 fig1G_generality_matrix.py
```

Each script checks `./data/` first, then the current folder, and raises a clear
`FileNotFoundError` naming exactly what's missing.

## verification/

Supporting scripts — not figure-producing.

| Script | Purpose |
|---|---|
| `build_GOBP_annotations_CAM_IMT.py` | GO Biological Process annotation for the bubble plots, from a QuickGO export. |
| `build_isolated_mito_annotation_from_curated_lists.py` | MitoCarta annotation + reference gene list for Fig 2B, from a manually curated protein list. |
| `build_MitoCarta_annotation_isolated_mito.py` | Earlier version of the script above, derived from the raw public MitoCarta3.0 database instead. Kept for the audit trail. |

## Data

EV tables and other input data were generated for this study and are available as
described in the Data Availability statement of the manuscript. They are not part of
this repository.

## Environment

Python >= 3.10 with pandas, numpy, matplotlib, seaborn, scipy, openpyxl, xlrd and
gseapy (see `requirements.txt`).
Tested with Python 3.11.15, pandas 3.0.2, numpy 2.4.4, matplotlib 3.10.9, seaborn
0.13.2, scipy 1.17.1, openpyxl 3.1.5, xlrd 2.0.2, gseapy 1.3.1 — and separately with
gseapy 1.1.9, the version the GSEA scripts were originally written against.

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
