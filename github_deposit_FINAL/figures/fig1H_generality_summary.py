#!/usr/bin/env python3
"""
Fig 1H -- generality-summary figure, single self-contained script (no intermediate
cache files).

Inputs -- two kinds (same split as the EV1O correlation script and Fig 1G matrix script):

(1) EV tables:
      EV_Table_8*.xlsx   -- MRC5 CAM/IMT proteomics
      EV_Table_10*.xlsx  -- 30-model TIS meta-analysis (Anerillas)

(2) Not part of any EV table -- annotation / category-definition files, must be bundled
    as their own standalone files in the deposit:
      Supplementary_PanelE_Data_1.xlsx
      analMito_Anerillas2026.xlsx
      analMito_Payea2024.xlsx

Note (same as EV1O and Fig 1G): EV_Table_8 has no 'MC3.0 Mito Pathways' column of its
own, so the Mito_gene_expression category for MRC5 reuses the MitoCarta3.0 annotation
loaded from analMito_Anerillas2026.xlsx instead (the same universal gene -> pathway
mapping used for the other models, not model-specific).

Output:
  Generality_summary.png, Generality_summary.svg, generality_summary_table.csv

For each of the 8 curated categories, across all 30 models (2 reference + 14 SenCat
cell lines x CTIS/IRIS), shows a diverging stacked bar: what fraction of models trend
up vs down, and what fraction of those are individually significant (permutation
test, p<0.05).

MRC5 is sourced from 'Log2FC Sen CTRL vs Prolif CTRL' (EV_Table_8), the direct analog of
"senescent vs proliferating" used for every other model; CAM and IMT are additional
processing arms in that file that are not used here. The table's own per-protein
p-value/q-value columns are not used either -- category-level significance always comes
from a competitive gene-resampling permutation test on the Log2FC values, applied
identically to every model (MRC5, Payea, SenCat alike), so MRC5 is on the same
statistical footing as everything else in the figure.
"""
import glob
import os
import re
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.getcwd()
print('Running from:', HERE)
DATA_DIR = 'data'


def find_file(patterns, label):
    """Find a required input by glob pattern(s): checks ./data/ first, then
    the current folder, and raises FileNotFoundError naming the missing file."""
    for pat in patterns:
        matches = glob.glob(f'{DATA_DIR}/{pat}') or glob.glob(pat)
        if matches:
            return matches[0]
    raise FileNotFoundError(
        f"Missing required input for {label}: none of {patterns} found in "
        f"./{DATA_DIR}/ or next to this script. Check the filename/location and rerun."
    )


# ---- EV tables ----
CAM_IMT_PATH = find_file(
    ['EV_Table_8*.xlsx'],
    'MRC5 CAM/IMT proteomics (Sen CTRL vs Prolif CTRL)')
ANERILLAS_PATH = find_file(
    ['EV_Table_10*.xlsx'],
    '30-model TIS meta-analysis (Anerillas)')

# ---- (2) annotation / category-definition files -- NOT EV tables, must be added ---
PANELE_PATH = find_file(
    ['Supplementary_PanelE_Data_1.xlsx'],
    'canonical per-pathway gene lists -- NOT an EV table, needs to be added as a '
    'standalone deposited file')
ANERILLAS_MITOCARTA_PATH = find_file(
    ['analMito_Anerillas2026.xlsx'],
    'MitoCarta3.0 pathway annotations for the 30-model gene universe '
    '-- NOT an EV table, needs to be added as a standalone deposited file')
PAYEA_PATH = find_file(
    ['analMito_Payea2024.xlsx'],
    'Payea et al. 2024 IMR90 (etoposide) published dataset '
    '-- NOT an EV table (external published data), needs to be added as a standalone '
    'deposited file')

print('All 5 required inputs found:')
for _label, _path in [('CAM_IMT_PATH', CAM_IMT_PATH), ('ANERILLAS_PATH', ANERILLAS_PATH),
                       ('PANELE_PATH', PANELE_PATH),
                       ('ANERILLAS_MITOCARTA_PATH', ANERILLAS_MITOCARTA_PATH),
                       ('PAYEA_PATH', PAYEA_PATH)]:
    print(f'  {_label}: {_path}')
print()

# ============================================================== 1. build gene-level data
RENAME = {
    'ATP5A1': 'ATP5F1A', 'ATP5B': 'ATP5F1B', 'ATP5C1': 'ATP5F1C', 'ATP5D': 'ATP5F1D',
    'ATP5F1': 'ATP5PB', 'ATP5H': 'ATP5PD', 'ATP5I': 'ATP5ME', 'ATP5J': 'ATP5PF',
    'ATP5J2': 'ATP5MF', 'ATP5L': 'ATP5MG', 'ATP5O': 'ATP5PO',
    'UQCRFS1;UQCRFS1P1': 'UQCRFS1',
    'SQRDL': 'SQOR',
}

ONEC_FOLATE_GENES = ['ALDH1L1', 'MTHFD1', 'MTHFR', 'SHMT1', 'DHFR', 'TYMS',
                      'MTHFD1L', 'MTHFD2', 'ALDH1L2', 'SHMT2', 'GART', 'ATIC']

CATS_SIMPLE = ['OXPHOS_subunits', 'Fatty_acid_oxidation', 'BCAA_metabolism',
               'NEAA_metabolism', 'Sulfur_metabolism', 'Methionine_metabolism_GOBP']

MITO_GE_PATTERN = re.compile(
    r'Mitochondrial central dogma > (?:Translation|mtDNA maintenance|mtRNA metabolism)')


def rn(g):
    return RENAME.get(g, g)


records = []
model_background = {}

# ---- canonical gene lists
gene_lists = {}
for cat in CATS_SIMPLE:
    d = pd.read_excel(PANELE_PATH, sheet_name=cat, header=3)
    genes = [rn(g) for g in d['Gene'].dropna().tolist() if g not in ('Mean', 'SEM')]
    gene_lists[cat] = genes
gene_lists['OneC_folate_metabolism'] = [rn(g) for g in ONEC_FOLATE_GENES]
print('Canonical gene-list sizes:')
for k, v in gene_lists.items():
    print(' ', k, len(v))

# ---- MRC5, from EV_Table_8 (CAM/IMT proteomics), Sen CTRL vs Prolif CTRL only
# sheet_name=0 (by position) + column-whitespace strip: prior EV tables in this project
# have shown up with renamed sheets / stray leading spaces in different exported copies.
cam = pd.read_excel(CAM_IMT_PATH, sheet_name=0)
cam.columns = cam.columns.str.strip()
cam = cam.dropna(subset=['Gene names']).drop_duplicates(subset=['Gene names'], keep='first')
cam['Gene names'] = cam['Gene names'].map(rn)
cam = cam.set_index('Gene names')
own_fc = cam['Log2FC Sen CTRL vs Prolif CTRL']
own_bg_mean = own_fc.dropna().mean()

# Shared MitoCarta3.0 pathway lookup, loaded once here so it can be reused below for
# the 30-model set too -- see the note in this script's docstring: EV_Table_8 itself
# does not carry this column, so the same universal (not model-specific) annotation
# is reused for it.
mito_col = 'MitoCarta3.0_MitoPathways'
mitocarta_annot = pd.read_excel(ANERILLAS_MITOCARTA_PATH, sheet_name='MitoCarta_annotations')
mitocarta_annot.columns = mitocarta_annot.columns.str.strip()
mitocarta_annot['Gene Symbol'] = mitocarta_annot['Gene Symbol'].map(rn)
mito_pathway_lookup = (mitocarta_annot.set_index('Gene Symbol')[mito_col]
                        if mito_col in mitocarta_annot.columns else None)
if mito_pathway_lookup is not None:
    mito_pathway_lookup = mito_pathway_lookup[~mito_pathway_lookup.index.duplicated()]

MRC5_LABEL = 'MRC5 (CAM/IMT)'
model_background[MRC5_LABEL] = (own_fc.dropna() - own_bg_mean).values
print('\nMRC5 full-proteome background size:', len(model_background[MRC5_LABEL]))
for cat, genes in gene_lists.items():
    for g in genes:
        if g in own_fc.index and pd.notna(own_fc[g]):
            records.append(dict(model=MRC5_LABEL, category=cat, gene=g,
                                 log2FC_rel=own_fc[g] - own_bg_mean))

ge_genes_mrc5 = []
if mito_pathway_lookup is not None:
    ge_genes_mrc5 = [g for g in own_fc.dropna().index
                     if g in mito_pathway_lookup.index and isinstance(mito_pathway_lookup[g], str)
                     and MITO_GE_PATTERN.search(mito_pathway_lookup[g])]
for g in ge_genes_mrc5:
    records.append(dict(model=MRC5_LABEL, category='Mito_gene_expression',
                         gene=g, log2FC_rel=own_fc[g] - own_bg_mean))
print('MRC5 Mito_gene_expression genes:', len(ge_genes_mrc5),
      '(sourced from analMito_Anerillas2026.xlsx, reused for the CAM/IMT gene '
      'universe -- EV_Table_8 has no MitoCarta pathway column of its own.)')

# ---- Payea 2024 (etoposide-TIS vs cycling, IMR90)
payea = pd.read_excel(PAYEA_PATH, sheet_name='Data')
payea['Symbol'] = payea['Symbol'].map(rn)
cyc_cols = ['Cyc_1', 'Cyc_2', 'Cyc_3']
etop_cols = ['Etop_1', 'Etop_2', 'Etop_3']
cyc_mean = payea[cyc_cols].mean(axis=1)
etop_mean = payea[etop_cols].mean(axis=1)
valid = (cyc_mean > 0) & (etop_mean > 0)
payea_fc_raw = pd.Series(np.log2(etop_mean[valid] / cyc_mean[valid]).values, index=payea.loc[valid, 'Symbol'].values)
payea_fc_raw = payea_fc_raw[~payea_fc_raw.index.duplicated()]
payea_bg_mean = payea_fc_raw.mean()

PAYEA_LABEL = 'IMR90 (etoposide, Payea 2024)'
model_background[PAYEA_LABEL] = (payea_fc_raw.dropna() - payea_bg_mean).values
print('Payea full-proteome background size:', len(model_background[PAYEA_LABEL]))
for cat, genes in gene_lists.items():
    for g in genes:
        if g in payea_fc_raw.index and pd.notna(payea_fc_raw[g]):
            records.append(dict(model=PAYEA_LABEL, category=cat, gene=g,
                                 log2FC_rel=payea_fc_raw[g] - payea_bg_mean))

mito_col = 'MitoCarta3.0_MitoPathways'
payea_paths = payea.set_index('Symbol')[mito_col]
payea_paths = payea_paths[~payea_paths.index.duplicated()]
ge_genes = [g for g in payea_fc_raw.index
            if g in payea_paths.index and isinstance(payea_paths[g], str) and MITO_GE_PATTERN.search(payea_paths[g])]
for g in ge_genes:
    records.append(dict(model=PAYEA_LABEL, category='Mito_gene_expression', gene=g,
                         log2FC_rel=payea_fc_raw[g] - payea_bg_mean))
print('Payea Mito_gene_expression genes:', len(ge_genes))

# ---- SenCat / Anerillas (14 cell lines x CTIS/IRIS)
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
    for cond in ['CTIS', 'IRIS']:
        diff_col = f"Student's T-test Difference {ct}_{cond}_{ct}_P"
        if diff_col not in d.columns:
            continue
        diff = d[diff_col]
        bg_mean = diff.mean()
        model_label = f'{ct} {cond}'
        model_background[model_label] = (diff.dropna() - bg_mean).values
        for cat, genes in gene_lists.items():
            for g in genes:
                if g in diff.index and pd.notna(diff[g]):
                    records.append(dict(model=model_label, category=cat, gene=g,
                                         log2FC_rel=diff[g] - bg_mean))
        if mito_pathway_lookup is not None:
            ge_genes_ct = [g for g in diff.index
                           if g in mito_pathway_lookup.index and isinstance(mito_pathway_lookup[g], str)
                           and MITO_GE_PATTERN.search(mito_pathway_lookup[g])]
            for g in ge_genes_ct:
                if pd.notna(diff[g]):
                    records.append(dict(model=model_label, category='Mito_gene_expression', gene=g,
                                         log2FC_rel=diff[g] - bg_mean))

gene_df = pd.DataFrame.from_records(records)
print('\ntotal gene-level records:', len(gene_df))
print('models built:', sorted(gene_df.model.unique()))
assert MRC5_LABEL in gene_df.model.unique(), f'MRC5 (CAM/IMT) missing from gene_df -- check {CAM_IMT_PATH}'

# ============================================================== 2. competitive permutation test
N_PERM = 2000
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
print('\npermutation test done:', perm.shape)
print(perm[perm.model == MRC5_LABEL].sort_values('category').to_string(index=False))

# ============================================================== 3. summary figure
CATS_ORDER = ['OXPHOS_subunits', 'Fatty_acid_oxidation', 'BCAA_metabolism',
              'NEAA_metabolism', 'OneC_folate_metabolism', 'Sulfur_metabolism',
              'Methionine_metabolism_GOBP', 'Mito_gene_expression']
CAT_LABELS = {
    'OXPHOS_subunits': 'OXPHOS subunits',
    'Fatty_acid_oxidation': 'Fatty acid oxidation',
    'BCAA_metabolism': 'BCAA metabolism',
    'NEAA_metabolism': 'NEAA metabolism',
    'OneC_folate_metabolism': '1C metabolism (folate)',
    'Sulfur_metabolism': 'Sulfur metabolism',
    'Methionine_metabolism_GOBP': 'Methionine metabolism',
    'Mito_gene_expression': 'Mito translation, mtDNA/mtRNA metab.',
}

rows = []
for cat in CATS_ORDER:
    sub = perm[perm.category == cat]
    n = len(sub)
    sig_up = ((sub.obs_mean > 0) & (sub.perm_p < 0.05)).sum()
    ns_up = ((sub.obs_mean > 0) & (sub.perm_p >= 0.05)).sum()
    ns_down = ((sub.obs_mean < 0) & (sub.perm_p >= 0.05)).sum()
    sig_down = ((sub.obs_mean < 0) & (sub.perm_p < 0.05)).sum()
    rows.append(dict(category=cat, n=n, sig_up=sig_up, ns_up=ns_up, ns_down=ns_down, sig_down=sig_down,
                      pct_positive=100 * (sig_up + ns_up) / n, pct_sig_up=100 * sig_up / n))
summary = pd.DataFrame(rows)
summary.to_csv('generality_summary_table.csv', index=False)
print('\n' + summary.to_string(index=False))

INK, MUTED, GRID = '#0b0b0b', '#6b6a66', '#e1e0d9'
DOWN_COLOR_SIG, DOWN_COLOR_NS = '#1a4d8f', '#a9c3e3'
UP_COLOR_SIG, UP_COLOR_NS = '#b5341f', '#eab8ac'

fig, ax = plt.subplots(figsize=(9.8, 5.2))
y_pos = np.arange(len(CATS_ORDER))[::-1]

for y, cat in zip(y_pos, CATS_ORDER):
    r = summary[summary.category == cat].iloc[0]
    n = r.n
    sig_down_frac = 100 * r.sig_down / n
    ns_down_frac = 100 * r.ns_down / n
    ns_up_frac = 100 * r.ns_up / n
    sig_up_frac = 100 * r.sig_up / n
    ax.barh(y, -ns_down_frac, left=0, color=DOWN_COLOR_NS, height=0.62, zorder=3)
    ax.barh(y, -sig_down_frac, left=-ns_down_frac, color=DOWN_COLOR_SIG, height=0.62, zorder=3)
    ax.barh(y, ns_up_frac, left=0, color=UP_COLOR_NS, height=0.62, zorder=2)
    ax.barh(y, sig_up_frac, left=ns_up_frac, color=UP_COLOR_SIG, height=0.62, zorder=3)
    ax.text(105, y, f'{r.pct_sig_up:.0f}% sig. up  (n={n})', fontsize=8.3, color=INK, va='center', ha='left')

ax.axvline(0, color=INK, lw=1.0, zorder=4)
ax.set_yticks(y_pos)
ax.set_yticklabels([CAT_LABELS[c] for c in CATS_ORDER], fontsize=9.5, color=INK)
ax.set_xlim(-100, 155)
ax.set_xticks([-100, -50, 0, 50, 100])
ax.set_xticklabels(['100%', '50%', '0', '50%', '100%'], fontsize=8.5, color=MUTED)
ax.set_xlabel('% of the 30 models  (down ←       → up)', fontsize=9, color=MUTED)
for s in ('top', 'right', 'left'):
    ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color(GRID)
ax.tick_params(axis='y', length=0)
ax.grid(axis='x', color=GRID, lw=0.5, zorder=0)
ax.axhline(len(CATS_ORDER) - 3.5, color=GRID, lw=0.9, zorder=1)

from matplotlib.patches import Patch
handles = [
    Patch(facecolor=UP_COLOR_SIG, label='significant up (p<0.05)'),
    Patch(facecolor=UP_COLOR_NS, label='trend up (n.s.)'),
    Patch(facecolor=DOWN_COLOR_NS, label='trend down (n.s.)'),
    Patch(facecolor=DOWN_COLOR_SIG, label='significant down (p<0.05)'),
]
ax.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.42, -0.34), ncol=4, frameon=False, fontsize=8)

fig.text(0.015, 0.975, 'Direction and significance of change across all 30 senescence models', fontsize=13.5,
          fontweight='bold', color=INK, ha='left', va='top')
subtitle = ('2 reference models (MRC5 CAM/IMT proteomics, Payea 2024 IMR90 etoposide) + 14 SenCat cell lines x '
            'CTIS/IRIS. Competitive gene-resampling permutation test, 2000 resamples per model x category.')
y_sub = 0.925
for ln in textwrap.wrap(subtitle, width=100):
    fig.text(0.015, y_sub, ln, fontsize=8.5, color=MUTED, ha='left', va='top')
    y_sub -= 0.038

fig.subplots_adjust(left=0.34, right=0.80, top=0.80, bottom=0.22)
fig.savefig('Generality_summary.png', dpi=300, facecolor='white')
fig.savefig('Generality_summary.svg', facecolor='white')
print('\nsaved Generality_summary.png/svg, generality_summary_table.csv to:', HERE)