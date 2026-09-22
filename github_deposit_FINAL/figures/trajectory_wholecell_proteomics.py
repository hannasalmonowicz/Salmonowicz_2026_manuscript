#!/usr/bin/env python3
# =====================================================================
# Fig 1E — OXPHOS vs metabolic pathways, time-resolved whole-cell trajectories
# LFQ-MS Time-Resolved Proteomics, whole cell
#
# Reads EV_Table_1_Time_resolved_LFQ_MS.xlsx (single sheet, header on row 1):
#   'Gene_names'  -> 'Gene names'
#   'FC_day1_Pro' -> 'Log2FC Sen D1 vs Prolif'
#   'FC_day3_Pro' -> 'Log2FC Sen D3 vs Prolif'
#   'FC_day9_Pro' -> 'Log2FC Sen D9 vs Prolif'
#
# PROTEIN_LISTS_FILE ('Figure_1_protein_lists_Whole Cell.xlsx') supplies the
# curated gene-category lists this script plots, via dedicated per-pathway
# sheets with a 'Gene' column (Verified_mitoproteome_WC=500,
# OXPHOS_subunits_WC=70, Fatty_acid_oxidation_WC=61, BCAA_metabolism_WC=19,
# Sulfur_metabolism=6, Nucleotide_metabolism=25, 1C_metabolism=4). It is not
# part of EV_Table_1 (which has only the quantitative Log2FC/LFQ data) and
# needs to be submitted as its own supplementary/EV file.
#
# Note: METHIONINE_GENES below is a hardcoded 12-gene set, not read from this
# file. This file's own 'Methionine_metabolism' sheet has only 5 genes; the
# 12-gene set is the correct one (matches the published figure).
# =====================================================================
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ── Custom Matplotlib Settings for Journal Publication ────────────────────────
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ── File Paths ────────────────────────────────────────────────────────────────
DATA_DIR = 'data'
# 'EV_Table_1_*' (with the trailing underscore), not 'EV_Table_1*' -- the deposit's
# data/ folder also holds EV_Table_10_meta_analysis_30_TIS_models.xlsx, which an
# unanchored 'EV_Table_1*' pattern also matches (it starts with 'EV_Table_1'), and
# glob() gives no ordering guarantee between the two files.
ev1_matches = glob.glob(f'{DATA_DIR}/EV_Table_1_*.xlsx') or glob.glob('EV_Table_1_*.xlsx')
if not ev1_matches:
    raise FileNotFoundError(
        "No 'EV_Table_1_*.xlsx' found. Put the submitted EV Table 1 "
        "(Time-resolved LFQ-MS) in ./data/ or next to this script."
    )
ANNOTATIONS_FILE = ev1_matches[0]

list_matches = (glob.glob(f'{DATA_DIR}/Figure_1_protein_lists_Whole*.xlsx')
                or glob.glob('Figure_1_protein_lists_Whole*.xlsx'))
if not list_matches:
    raise FileNotFoundError(
        "No 'Figure_1_protein_lists_Whole*.xlsx' found. Put it in ./data/ "
        "or next to this script — it supplies the curated gene-category "
        "lists (MitoCarta subsets etc.) and is not part of EV_Table_1."
    )
PROTEIN_LISTS_FILE = list_matches[0]
OUT_DIR = '.'

# ── Load Data ─────────────────────────────────────────────────────────────────
# Read by position, not by name — this file only ever has one data sheet, so
# grabbing it by position (0) avoids hardcoding a sheet name that can drift
# between exported copies.
df = pd.read_excel(ANNOTATIONS_FILE, sheet_name=0)
# Strip stray leading/trailing whitespace from column names — this file can have
# columns like ' Log2FC Sen D1 vs Prolif' with a leading space, which breaks an
# exact-name lookup. Stripping here means it works whether or not a given export
# has that space.
df.columns = df.columns.str.strip()
xls_lists = pd.ExcelFile(PROTEIN_LISTS_FILE)

def get_genes(sheet_name):
    """Extract gene list directly from dedicated pathway sheet in Excel."""
    if sheet_name in xls_lists.sheet_names:
        return set(pd.read_excel(xls_lists, sheet_name=sheet_name)['Gene'].dropna())
    print(f"WARNING: Sheet '{sheet_name}' not found in {PROTEIN_LISTS_FILE}")
    return set()

# ── Trajectory Helper ─────────────────────────────────────────────────────────
X = np.array([0, 1, 3, 9])
FC_MAP = {
    'Day 1': 'Log2FC Sen D1 vs Prolif',
    'Day 3': 'Log2FC Sen D3 vs Prolif',
    'Day 9': 'Log2FC Sen D9 vs Prolif',
}

def traj(genes):
    sub = df[df['Gene names'].isin(genes)]
    means = np.array([sub[FC_MAP[tp]].dropna().mean() for tp in ['Day 1', 'Day 3', 'Day 9']])
    sems = np.array([sub[FC_MAP[tp]].dropna().std(ddof=1) / np.sqrt(sub[FC_MAP[tp]].dropna().count())
                      for tp in ['Day 1', 'Day 3', 'Day 9']])
    return np.concatenate([[0.0], means]), np.concatenate([[0.0], sems]), len(sub)

# ── Load Gene Categories directly from Excel Sheets ───────────────────────────
ALL_PROTEINS = set(df['Gene names'].dropna())
MITO_REF = get_genes('Verified_mitoproteome_WC')        # n=500
OXPHOS_GENES = get_genes('OXPHOS_subunits_WC')          # n=70
NEAA = get_genes('NEAA_metabolism')                     # n=14 matched
FAO = get_genes('Fatty_acid_oxidation_WC')               # n=61
BCAA = get_genes('BCAA_metabolism_WC')                    # n=19
SULFUR = get_genes('Sulfur_metabolism')                  # n=6
MITO_NUC_1C = get_genes('Nucleotide_metabolism') | get_genes('1C_metabolism')  # n=22

# GO:BP methionine metabolism
METHIONINE_GENES = {
    'MAT2A', 'MAT2B',                           # SAM synthesis
    'AHCY', 'AHCYL1', 'AHCYL2',                 # SAH hydrolysis
    'MTAP', 'ADI1', 'ENOPH1', 'APIP', 'MRI1',   # methionine salvage (5'-MTA pathway)
    'MSRA',                                     # methionine sulfoxide reductase
    'SMS',                                      # polyamine synthesis (SAM-dependent)
}

# Cytosolic nucleotide & 1C
CYTO_NUC_1C = {
    'ADSL', 'ALDH1L1', 'APRT', 'ATIC', 'CAD', 'CTPS1', 'CTPS2', 'DHFR',
    'DTYMK', 'GART', 'GMPS', 'HINT1', 'HPRT1', 'IMPDH1', 'IMPDH2',
    'MTHFD1', 'MTHFR', 'NT5C', 'NT5C2', 'NT5DC2', 'NUDT5', 'PAICS',
    'PFAS', 'PNP', 'PPAT', 'SHMT1', 'TYMS', 'UCK2', 'UMPS',
}

# ── Panel Configuration ───────────────────────────────────────────────────────
CATEGORY_SPECS = [
    ('All proteins', ALL_PROTEINS, '--', None, 'ref'),
    ('MitoCarta 3.0: all mitochondria', MITO_REF, '--', None, 'ref'),
    ('MitoCarta 3.0: BCAA metabolism', BCAA, '-', 'D', 'main'),
    ('MitoCarta 3.0: fatty acid oxidation', FAO, '-', 'D', 'main'),
    ('MitoCarta 3.0: sulfur metabolism', SULFUR, '-', 'v', 'main'),
    ('MitoCarta 3.0: OXPHOS subunits', OXPHOS_GENES, '-', 'o', 'main'),
    ('MitoCarta 3.0: mito nucleotide & 1C', MITO_NUC_1C, '--', 's', 'main'),
    ('Curated: NEAA metabolism', NEAA, '-', '^', 'de-emphasized'),
    ('GO:BP methionine metabolism', METHIONINE_GENES, '-', 'v', 'main'),
    ('GO:BP & curated: cytosolic nucleotide & 1C', CYTO_NUC_1C, '-', 's', 'main'),
]

REF_COLORS = {'All proteins': '#DADAD7', 'MitoCarta 3.0: all mitochondria': '#B0B0AB'}

PALETTES = {
    'Palette1_OriginalStyle': dict(
        title='Palette 1: Original-style',
        colors={
            'MitoCarta 3.0: BCAA metabolism': '#9CCC65',
            'MitoCarta 3.0: fatty acid oxidation': '#1B5E20',
            'MitoCarta 3.0: sulfur metabolism': '#00ACC1',
            'MitoCarta 3.0: OXPHOS subunits': '#C0392B',
            'MitoCarta 3.0: mito nucleotide & 1C': '#5E35B1',
            'Curated: NEAA metabolism': '#B39DDB',              # light purple
            'GO:BP methionine metabolism': '#34495E',           # bluish, very dark grey
            'GO:BP & curated: cytosolic nucleotide & 1C': '#1565C0',  # rich blue, not purple
        }),
}

# ── Gene manifest -- print the actual curated gene list behind every line ──
print('=' * 78)
print('CURATED GENE LISTS BEHIND EACH LINE')
print('=' * 78)
for label, genes, *_ in CATEGORY_SPECS:
    detected = sorted(df[df['Gene names'].isin(genes)]['Gene names'].unique())
    if label == 'All proteins':
        print(f'\n{label}: n={len(detected)} (entire detected proteome, not listed)')
        continue
    print(f'\n{label}  (n={len(detected)}):')
    print('  ' + ', '.join(detected))
print('=' * 78)

# ── Drawing helper -- produces ribbon figure for Palette 1 ───────────────────
def draw_panel(palette_key):
    palette = PALETTES[palette_key]
    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=(8.5, 7.0))
    ax.axhline(0, color='#DDDDDD', lw=0.9, zorder=0)
    handles = []
    for label, genes, ls, mk, ctype in CATEGORY_SPECS:
        is_ref = ctype == 'ref'
        is_deemph = ctype == 'de-emphasized'
        col = REF_COLORS[label] if is_ref else palette['colors'][label]
        m, s, n = traj(genes)
        if is_ref:
            lw = 2.2
        elif is_deemph:
            lw = 1.6
        else:
            lw = 2.5

        if mk is None:
            ms = 0
        elif is_deemph:
            ms = 4.0
        else:
            ms = 6.0

        if is_ref:
            alpha = 0.6
        elif is_deemph:
            alpha = 0.85
        else:
            alpha = 1.0

        (ln,) = ax.plot(
            X, m, color=col, lw=lw, ls=ls, marker=mk, markersize=ms,
            markerfacecolor=col, markeredgecolor='white', markeredgewidth=0.6,
            label=f'{label} (n={n})', alpha=alpha, zorder=2 if is_ref else 5,
        )
        if not is_ref:
            ax.fill_between(X, m - s, m + s, color=col, alpha=0.10, zorder=1)
        handles.append(ln)

    ax.set_xticks([0, 1, 3, 9])
    ax.set_xticklabels(['Prolif', 'Day 1', 'Day 3', 'Day 9'], fontsize=11, fontweight='bold')
    ax.set_ylabel('Mean log$_2$ FC (senescent vs proliferative)', fontsize=11, fontweight='bold')
    ax.text(0.02, 0.95, 'OXPHOS vs metabolic pathways', transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='top')
    ax.text(0.02, 0.90, f'LFQ-MS Time-Resolved Proteomics, whole cell -- {palette["title"]}',
            transform=ax.transAxes, fontsize=10, style='italic', color='#444444', va='top')
    ax.legend(handles=handles, fontsize=8, loc='upper left', bbox_to_anchor=(0.0, -0.15),
              ncol=2, frameon=False)
    sns.despine(ax=ax)
    plt.tight_layout()
    return fig

# ── Render Palette 1 only, PNG + SVG ──────────────────────────────────────────
palette_key = 'Palette1_OriginalStyle'
fig = draw_panel(palette_key)
fig.savefig(f'{OUT_DIR}/PanelE_{palette_key}.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT_DIR}/PanelE_{palette_key}.svg', bbox_inches='tight')
plt.show()

print(f'Done. Saved PanelE_{palette_key}.png and PanelE_{palette_key}.svg')
