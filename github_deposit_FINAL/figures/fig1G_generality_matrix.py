#!/usr/bin/env python3
"""
Fig 1G -- two versions of the main 30-model dot-matrix, using a "1C metabolism (folate)"
category that includes GART and ATIC, differing only in which dataset backs the MRC5
column:

  Time-resolved -- MRC5 from Supplementary_PanelE_Data_1.xlsx (doxorubicin/Day 9 dataset).
  CAM/IMT       -- MRC5 from EV_Table_8, Sen CTRL vs Prolif CTRL comparison (CAM and IMT
                    are separate processing arms in that file; CTRL is the arm plotted
                    here, as the direct analog of "senescent vs proliferating" used for
                    every other model).

Every other column (IMR90/Payea, all 14 Anerillas SenCat cell lines x CTIS/IRIS) is
identical between the two versions -- same gene lists, same source data -- so any visual
difference between them comes only from the MRC5 dataset used.

Inputs -- two kinds (same split as the EV1O correlation script):

(1) EV tables:
      EV_Table_8*.xlsx   -- MRC5 CAM/IMT proteomics
      EV_Table_10*.xlsx  -- all 14 SenCat lines

(2) Not part of any EV table -- annotation / category-definition files, bundled as their
    own standalone files in the deposit:
      Supplementary_PanelE_Data_1.xlsx  -- gene lists + MRC5 (time-resolved) proteome
      5_Annotations.xlsx                -- Mito_gene_expression flags for MRC5 (time-resolved) only
      analMito_Anerillas2026.xlsx       -- MitoCarta3.0 pathway annotations (SenCat + reused for CAM/IMT)
      analMito_Payea2024.xlsx           -- Payea et al. 2024 published dataset

Note (same as EV1O): EV_Table_8 has no 'MC3.0 Mito Pathways' column of its own, so the
Mito_gene_expression category for the MRC5 (CAM/IMT) column reuses the MitoCarta3.0
annotation loaded from analMito_Anerillas2026.xlsx instead (the same universal gene ->
pathway mapping used for the other models, not CAM/IMT-specific). The time-resolved MRC5
column keeps its own, more curated 5_Annotations.xlsx boolean-flag annotation, unchanged.
"""
import glob
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import textwrap

RENAME = {
    'ATP5A1': 'ATP5F1A', 'ATP5B': 'ATP5F1B', 'ATP5C1': 'ATP5F1C', 'ATP5D': 'ATP5F1D',
    'ATP5F1': 'ATP5PB', 'ATP5H': 'ATP5PD', 'ATP5I': 'ATP5ME', 'ATP5J': 'ATP5PF',
    'ATP5J2': 'ATP5MF', 'ATP5L': 'ATP5MG', 'ATP5O': 'ATP5PO',
    'UQCRFS1;UQCRFS1P1': 'UQCRFS1', 'SQRDL': 'SQOR',
}
rn = lambda g: RENAME.get(g, g)
MITO_GE_PATTERN = re.compile(r'Mitochondrial central dogma > (?:Translation|mtDNA maintenance|mtRNA metabolism)')
N_PERM = 2000
DATA_DIR = 'data'


def find_file(patterns, label):
    """Find a required input by glob pattern(s), checking ./data/ first, then the
    current folder. Raises loudly and by name if not found -- never silently skips
    a required input."""
    for pat in patterns:
        matches = glob.glob(f'{DATA_DIR}/{pat}') or glob.glob(pat)
        if matches:
            return matches[0]
    raise FileNotFoundError(
        f"Missing required input for {label}: none of {patterns} found in "
        f"./{DATA_DIR}/ or next to this script. Check the filename/location and rerun."
    )


# ---- EV tables ----
CAM_IMT_PATH = find_file(['EV_Table_8*.xlsx'], 'MRC5 (CAM/IMT)')
ANERILLAS_PATH = find_file(['EV_Table_10*.xlsx'], '30-model TIS meta-analysis (Anerillas)')

# ---- annotation / category-definition files -- not EV tables, must be added ----
PANELE_PATH = find_file(
    ['Supplementary_PanelE_Data_1.xlsx'],
    'canonical per-pathway gene lists + MRC5 (time-resolved) proteome '
    '-- not an EV table, needs to be added as a standalone deposited file')
ANNOTATIONS_PATH = find_file(
    ['5_Annotations.xlsx', '5_Annotations_copy.xlsx'],
    'Mito_gene_expression pathway flags for MRC5 (time-resolved) only '
    '-- not an EV table, needs to be added as a standalone deposited file')
PAYEA_PATH = find_file(
    ['analMito_Payea2024.xlsx'],
    'Payea et al. 2024 IMR90 (etoposide) published dataset '
    '-- NOT an EV table (external published data), needs to be added as a standalone '
    'deposited file')
ANERILLAS_MITOCARTA_PATH = find_file(
    ['analMito_Anerillas2026.xlsx'],
    'MitoCarta3.0 pathway annotations for the 30-model gene universe '
    '-- NOT an EV table, needs to be added as a standalone deposited file')

print('All required inputs found:')
for _label, _path in [('CAM_IMT_PATH', CAM_IMT_PATH), ('ANERILLAS_PATH', ANERILLAS_PATH),
                       ('PANELE_PATH', PANELE_PATH), ('ANNOTATIONS_PATH', ANNOTATIONS_PATH),
                       ('PAYEA_PATH', PAYEA_PATH),
                       ('ANERILLAS_MITOCARTA_PATH', ANERILLAS_MITOCARTA_PATH)]:
    print(f'  {_label}: {_path}')
print()

# Shared MitoCarta3.0 pathway lookup, loaded once, reused for both the SenCat set and
# the MRC5 CAM/IMT column, since EV_Table_8 has no such column itself (see note above).
mito_col = 'MitoCarta3.0_MitoPathways'
_mitocarta_annot = pd.read_excel(ANERILLAS_MITOCARTA_PATH, sheet_name='MitoCarta_annotations')
_mitocarta_annot.columns = _mitocarta_annot.columns.str.strip()
_mitocarta_annot['Gene Symbol'] = _mitocarta_annot['Gene Symbol'].map(rn)
mito_pathway_lookup = (_mitocarta_annot.set_index('Gene Symbol')[mito_col]
                        if mito_col in _mitocarta_annot.columns else None)
if mito_pathway_lookup is not None:
    mito_pathway_lookup = mito_pathway_lookup[~mito_pathway_lookup.index.duplicated()]

CATS_SIMPLE = ['OXPHOS_subunits', 'Fatty_acid_oxidation', 'BCAA_metabolism',
               'NEAA_metabolism', 'Sulfur_metabolism']
ONEC_FOLATE_GENES = ['ALDH1L1', 'MTHFD1', 'MTHFR', 'SHMT1', 'DHFR', 'TYMS',
                      'MTHFD1L', 'MTHFD2', 'ALDH1L2', 'SHMT2', 'GART', 'ATIC']  # +GART, +ATIC per PI

gene_lists = {}
for cat in CATS_SIMPLE:
    d = pd.read_excel(PANELE_PATH, sheet_name=cat, header=3)
    gene_lists[cat] = [rn(g) for g in d['Gene'].dropna().tolist() if g not in ('Mean', 'SEM')]
gene_lists['OneC_folate_metabolism'] = [rn(g) for g in ONEC_FOLATE_GENES]
print('1C-folate category now has', len(gene_lists['OneC_folate_metabolism']), 'genes:', gene_lists['OneC_folate_metabolism'])


def records_for_series(fc_series, model_label, group):
    """fc_series: index=gene symbol (already renamed), values=raw log2FC. Returns
    (records list, full-proteome background array) using the same background-correction
    logic as build_data.py."""
    fc_series = fc_series.dropna()
    bg_mean = fc_series.mean()
    corrected = fc_series - bg_mean
    recs = []
    for cat, genes in gene_lists.items():
        for g in genes:
            if g in corrected.index:
                recs.append(dict(model=model_label, group=group, category=cat, gene=g, log2FC_rel=corrected[g]))
    return recs, corrected.values


def add_mito_ge(recs, fc_series_corrected_full, model_label, group, ge_genes):
    for g in ge_genes:
        if g in fc_series_corrected_full.index:
            recs.append(dict(model=model_label, group=group, category='Mito_gene_expression',
                              gene=g, log2FC_rel=fc_series_corrected_full[g]))


all_records = []
model_background = {}

# ---------------------------------------------------------- MRC5, time-resolved (PanelE)
all_prot = pd.read_excel(PANELE_PATH, sheet_name='All_proteins', header=3)
all_prot = all_prot[~all_prot['Gene'].isin(['Mean', 'SEM'])]
all_prot['Gene'] = all_prot['Gene'].map(rn)
own_fc = all_prot.set_index('Gene')['FC_day9_Pro (log2FC)']
MRC5_TR_LABEL = 'MRC5 (time-resolved)'
recs, bg = records_for_series(own_fc, MRC5_TR_LABEL, 'Reference')
model_background[MRC5_TR_LABEL] = bg
own_bg_mean = own_fc.dropna().mean()
own_corrected = own_fc.dropna() - own_bg_mean
ann = pd.read_excel(ANNOTATIONS_PATH, sheet_name='ALL')
ann['Gene_names'] = ann['Gene_names'].map(rn)
mask = (ann['mtDNA_maintenance'] == 1) | (ann['mtRNA_metabolism'] == 1) | (ann['Translation'] == 1)
ge_genes_mrc5 = set(ann.loc[mask, 'Gene_names'].dropna().map(rn))
add_mito_ge(recs, own_corrected, MRC5_TR_LABEL, 'Reference', ge_genes_mrc5)
all_records += recs
print(MRC5_TR_LABEL, 'full-proteome background size:', len(bg))

# ---------------------------------------------------------- MRC5, CAM/IMT (EV_Table_8)
# sheet_name=0 (by position) + column-whitespace strip: prior EV tables in this project
# have shown up with renamed sheets / stray leading spaces in different exported copies.
cam = pd.read_excel(CAM_IMT_PATH, sheet_name=0)
cam.columns = cam.columns.str.strip()
cam = cam.dropna(subset=['Gene names']).drop_duplicates(subset=['Gene names'], keep='first')
cam['Gene names'] = cam['Gene names'].map(rn)
cam = cam.set_index('Gene names')
cam_fc = cam['Log2FC Sen CTRL vs Prolif CTRL']
MRC5_CI_LABEL = 'MRC5 (CAM/IMT)'
recs, bg = records_for_series(cam_fc, MRC5_CI_LABEL, 'Reference')
model_background[MRC5_CI_LABEL] = bg
cam_bg_mean = cam_fc.dropna().mean()
cam_corrected = cam_fc.dropna() - cam_bg_mean
# EV_Table_8 has no 'MC3.0 Mito Pathways' column -- reusing the shared
# mito_pathway_lookup (from analMito_Anerillas2026.xlsx) loaded above instead.
ge_genes_cam = set()
if mito_pathway_lookup is not None:
    ge_genes_cam = set(g for g in cam_fc.dropna().index
                        if g in mito_pathway_lookup.index and isinstance(mito_pathway_lookup[g], str)
                        and MITO_GE_PATTERN.search(mito_pathway_lookup[g]))
add_mito_ge(recs, cam_corrected, MRC5_CI_LABEL, 'Reference', ge_genes_cam)
all_records += recs
print(MRC5_CI_LABEL, 'full-proteome background size:', len(bg), '| Mito_gene_expression genes:', len(ge_genes_cam))

# ---------------------------------------------------------- IMR90 (Payea 2024) -- shared, unchanged source
payea = pd.read_excel(PAYEA_PATH, sheet_name='Data')
payea['Symbol'] = payea['Symbol'].map(rn)
cyc_mean = payea[['Cyc_1', 'Cyc_2', 'Cyc_3']].mean(axis=1)
etop_mean = payea[['Etop_1', 'Etop_2', 'Etop_3']].mean(axis=1)
valid = (cyc_mean > 0) & (etop_mean > 0)
payea_fc = pd.Series(np.log2(etop_mean[valid] / cyc_mean[valid]).values, index=payea.loc[valid, 'Symbol'].values)
payea_fc = payea_fc[~payea_fc.index.duplicated()]
PAYEA_LABEL = 'IMR90 (etoposide, Payea 2024)'
recs, bg = records_for_series(payea_fc, PAYEA_LABEL, 'Reference')
model_background[PAYEA_LABEL] = bg
payea_bg_mean = payea_fc.dropna().mean()
payea_corrected = payea_fc.dropna() - payea_bg_mean
payea_paths = payea.set_index('Symbol')['MitoCarta3.0_MitoPathways']
payea_paths = payea_paths[~payea_paths.index.duplicated()]
ge_genes_payea = set(g for g in payea_fc.index
                      if g in payea_paths.index and isinstance(payea_paths[g], str) and MITO_GE_PATTERN.search(payea_paths[g]))
add_mito_ge(recs, payea_corrected, PAYEA_LABEL, 'Reference', ge_genes_payea)
all_records += recs

# ---------------------------------------------------------- Anerillas SenCat -- shared, unchanged source
CELLTYPE_SHEETS = {'WI38': 'WI38', 'BJ': 'BJ', 'HSAEC': 'HSAEC', 'HEKn': 'HEKn', 'HCAEC': 'HCAEC',
                    'HUVEC': 'HUVEC', 'BMMSC': 'BMMSC', 'HVSMC': 'HVSMC', 'HSKM': 'HSKM', 'PBMC': 'PBMC',
                    'PreAdipo': 'PreAdipo', 'NHO': 'NHO', 'NHA': 'NHA', 'HEMn': 'HEM'}
# mito_pathway_lookup was already loaded once above (shared with the MRC5/CAM-IMT
# block) -- reused here for the 30-model set, same as the original script's sen_paths.

for ct, sheet in CELLTYPE_SHEETS.items():
    d = pd.read_excel(ANERILLAS_PATH, sheet_name=sheet)
    d.columns = d.columns.str.strip()
    d = d.dropna(subset=['Gene Symbol']).drop_duplicates(subset=['Gene Symbol'], keep='first')
    d['Gene Symbol'] = d['Gene Symbol'].map(rn)
    d = d.set_index('Gene Symbol')
    for cond, group in [('CTIS', 'CTIS'), ('IRIS', 'IRIS')]:
        diff_col = f"Student's T-test Difference {ct}_{cond}_{ct}_P"
        if diff_col not in d.columns:
            continue
        diff = d[diff_col]
        model_label = f'{ct} {cond}'
        recs, bg = records_for_series(diff, model_label, group)
        model_background[model_label] = bg
        bg_mean = diff.dropna().mean()
        corrected = diff.dropna() - bg_mean
        ge_genes = set()
        if mito_pathway_lookup is not None:
            ge_genes = set(g for g in diff.dropna().index
                            if g in mito_pathway_lookup.index and isinstance(mito_pathway_lookup[g], str)
                            and MITO_GE_PATTERN.search(mito_pathway_lookup[g]))
        add_mito_ge(recs, corrected, model_label, group, ge_genes)
        all_records += recs

gene_df = pd.DataFrame.from_records(all_records)
print('total gene-level records:', len(gene_df))

# ---------------------------------------------------------- competitive permutation test (all models)
rng = np.random.default_rng(0)
perm_rows = []
for (model, cat), sub in gene_df.groupby(['model', 'category']):
    vals = sub['log2FC_rel'].values
    n = len(vals)
    obs_mean = vals.mean()
    bg = model_background[model]
    B = len(bg)
    if n == 0 or n > B:
        perm_rows.append(dict(model=model, category=cat, n=n, obs_mean=np.nan, perm_p=np.nan))
        continue
    keys = rng.random((N_PERM, B))
    idx = np.argpartition(keys, n, axis=1)[:, :n]
    null_means = bg[idx].mean(axis=1)
    p = (np.sum(np.abs(null_means) >= abs(obs_mean)) + 1) / (N_PERM + 1)
    perm_rows.append(dict(model=model, category=cat, n=n, obs_mean=obs_mean, perm_p=p))
perm = pd.DataFrame(perm_rows)
print(perm[perm.category == 'OneC_folate_metabolism'].sort_values('model').to_string(index=False))

# ================================================================== plotting (same spec as build_generality_matrix.py)
CATS_ORDER = ['OXPHOS_subunits', 'Fatty_acid_oxidation', 'BCAA_metabolism',
              'NEAA_metabolism', 'OneC_folate_metabolism', 'Sulfur_metabolism', 'Mito_gene_expression']
CAT_LABELS = {
    'OXPHOS_subunits': 'OXPHOS subunits', 'Fatty_acid_oxidation': 'Fatty acid oxidation',
    'BCAA_metabolism': 'BCAA metabolism', 'NEAA_metabolism': 'NEAA metabolism',
    'OneC_folate_metabolism': '1C metabolism (folate)', 'Sulfur_metabolism': 'Sulfur metabolism',
    'Mito_gene_expression': 'Mito translation,\nmtDNA/mtRNA metab.',
}
CELLTYPE_ORDER = ['WI38', 'BJ', 'HSAEC', 'HEKn', 'HCAEC', 'HUVEC', 'BMMSC', 'HVSMC',
                  'HSKM', 'PBMC', 'PreAdipo', 'NHO', 'NHA', 'HEMn']
INK, MUTED, GRID, SURFACE = '#0b0b0b', '#52514e', '#e1e0d9', '#fcfcfb'
DIV_BLUE, DIV_MID, DIV_RED = '#2a78d6', '#f0efec', '#e34948'
cmap = LinearSegmentedColormap.from_list('div_blue_red', [DIV_BLUE, DIV_MID, DIV_RED], N=256)
VLIM = 1.0
norm = Normalize(vmin=-VLIM, vmax=VLIM)


def size_for_p(p):
    if pd.isna(p):
        return 18
    neglogp = -np.log10(max(p, 1 / 2001))
    return np.clip(18 + neglogp * 55, 18, 260)


def alpha_for_p(p):
    return 1.0 if (pd.notna(p) and p < 0.05) else 0.42


def build_figure(mrc5_label, mrc5_tag, title_tag, out_prefix):
    n_ref = 2
    n_ct = len(CELLTYPE_ORDER)
    col_order = [mrc5_label, 'IMR90 (etoposide, Payea 2024)'] + \
                [f'{ct} CTIS' for ct in CELLTYPE_ORDER] + [f'{ct} IRIS' for ct in CELLTYPE_ORDER]
    col_labels = ['MRC5', 'IMR90 (Payea et al. 2024)'] + CELLTYPE_ORDER + CELLTYPE_ORDER

    fig, ax = plt.subplots(figsize=(12.6, 5.6))
    row_y = {cat: len(CATS_ORDER) - 1 - i for i, cat in enumerate(CATS_ORDER)}
    for cat, y in row_y.items():
        if (len(CATS_ORDER) - 1 - y) % 2 == 0:
            ax.axhspan(y - 0.5, y + 0.5, color='#f7f7f5', zorder=0)

    for x, model in enumerate(col_order):
        for cat in CATS_ORDER:
            row = perm[(perm.model == model) & (perm.category == cat)]
            if row.empty or pd.isna(row.iloc[0].obs_mean):
                continue
            r = row.iloc[0]
            y = row_y[cat]
            color = cmap(norm(np.clip(r.obs_mean, -VLIM, VLIM)))
            ax.scatter(x, y, s=size_for_p(r.perm_p), color=color, alpha=alpha_for_p(r.perm_p),
                       edgecolor=SURFACE, linewidth=0.9, zorder=3)

    for boundary in (0.5, n_ref - 0.5, n_ref + n_ct - 0.5):
        ax.axvline(boundary, color=GRID, lw=1.1, zorder=1)
    SOURCE_Y, INDUCTION_Y = len(CATS_ORDER) + 0.55, len(CATS_ORDER) - 0.20
    ax.text((n_ref - 0.5 + n_ref + 2 * n_ct - 0.5) / 2, SOURCE_Y, 'Anerillas et al. 2026',
            fontsize=9.5, fontweight='bold', color=INK, ha='center', va='bottom')
    ax.text(0.5 - 0.15, INDUCTION_Y, mrc5_tag, fontsize=8.3, fontweight='bold', color=MUTED, ha='right', va='bottom')
    ax.text(0.5 + 0.15, INDUCTION_Y, 'Etoposide', fontsize=8.3, fontweight='bold', color=MUTED, ha='left', va='bottom')
    for label, x0, x1 in [('CTIS (DOX)', n_ref - 0.5, n_ref + n_ct - 0.5),
                           ('IRIS', n_ref + n_ct - 0.5, n_ref + 2 * n_ct - 0.5)]:
        ax.text((x0 + x1) / 2, INDUCTION_Y, label, fontsize=8.3, fontweight='bold', color=MUTED, ha='center', va='bottom')

    ax.set_xlim(-0.7, len(col_order) + 6.3)
    ax.set_ylim(-0.7, len(CATS_ORDER) + 1.05)
    ax.set_xticks(range(len(col_order)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=6.6, color=MUTED)
    ax.set_yticks(list(row_y.values()))
    ax.set_yticklabels([CAT_LABELS[c] for c in row_y.keys()], fontsize=9.3, color=INK)
    ax.tick_params(axis='both', length=0)
    for s in ('top', 'right', 'left', 'bottom'):
        ax.spines[s].set_visible(False)
    ax.grid(axis='x', color=GRID, lw=0.4, zorder=0, alpha=0.5)

    for cat, y in row_y.items():
        sub = perm[perm.category == cat]
        n = len(sub)
        sig_up = ((sub.obs_mean > 0) & (sub.perm_p < 0.05)).sum()
        sig_down = ((sub.obs_mean < 0) & (sub.perm_p < 0.05)).sum()
        if sig_up >= sig_down:
            txt, col = f'{100 * sig_up / n:.0f}% sig. up', DIV_RED
        else:
            txt, col = f'{100 * sig_down / n:.0f}% sig. down', DIV_BLUE
        ax.text(len(col_order) + 0.6, y, txt, fontsize=8.6, color=col, fontweight='bold', ha='left', va='center')

    cbar_ax = fig.add_axes([0.935, 0.30, 0.012, 0.32])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_ticks([-VLIM, 0, VLIM])
    cbar.set_ticklabels([f'≤-{VLIM:.0f}', '0', f'≥+{VLIM:.0f}'])
    cbar.ax.tick_params(labelsize=7.5, colors=MUTED, length=0)
    cbar.outline.set_visible(False)
    cbar_ax.set_title('log$_2$FC', fontsize=8, color=INK, pad=6)

    leg_handles = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=MUTED, markeredgecolor=SURFACE,
               markersize=11, alpha=1.0, label='significant (p<0.05)'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor=MUTED, markeredgecolor=SURFACE,
               markersize=6, alpha=0.42, label='not significant'),
    ]
    ax.legend(handles=leg_handles, loc='lower center', bbox_to_anchor=(0.30, -0.30), ncol=2,
              frameon=False, fontsize=8, handletextpad=0.6, columnspacing=1.4)

    fig.text(0.01, 0.99, f'Direction, magnitude, and significance -- MRC5 dataset: {title_tag}', fontsize=14,
              fontweight='bold', color=INK, ha='left', va='top')
    subtitle = ('All 29 non-MRC5 columns are identical between the two figure variants -- only the MRC5 column '
                'changes.')
    y_sub = 0.955
    for ln in textwrap.wrap(subtitle, width=140):
        fig.text(0.01, y_sub, ln, fontsize=8.7, color=MUTED, ha='left', va='top')
        y_sub -= 0.028

    fig.subplots_adjust(left=0.125, right=0.90, top=0.865, bottom=0.24)
    fig.savefig(f'{out_prefix}.png', dpi=300, facecolor=SURFACE)
    fig.savefig(f'{out_prefix}.svg', facecolor=SURFACE)
    plt.close(fig)
    print('saved', f'{out_prefix}.png/svg')


build_figure(MRC5_TR_LABEL, 'DOX', 'time-resolved (PanelE_v1)', 'Generality_matrix_time_resolved')
build_figure(MRC5_CI_LABEL, 'DOX', 'CAM/IMT (Sen CTRL vs Prolif CTRL)', 'Generality_matrix_CAM_IMT')