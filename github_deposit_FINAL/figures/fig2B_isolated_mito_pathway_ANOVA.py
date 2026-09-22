"""
Fig 2B -- isolated-mito abundance (left panel) + aggregation (right panel), Ordinary
ANOVA + Bonferroni. This is the exact same figure as before (both panels, one figure);
"abundance" and "aggregation" are the two halves of Figure2B_OrdinaryANOVA_Bonferroni,
not two separate figures.

Inputs -- two kinds:

(1) EV table:
      EV_Table_2*.xlsx -- 'Gene names', per-replicate Corrected Sol/Insol Pro/Sen,
      and 'log2 FC TOTAL Sen_Pro'.

(2) Not part of any EV table -- built by
    verification/build_isolated_mito_annotation_from_curated_lists.py from
    Figure_1_protein_lists_ISOLATED_MITOCHONDRIA.xlsx (a manually curated
    categorization, the isolated-mito analog of Figure_1_protein_lists_Whole_Cell.xlsx
    already used by the trajectory script), written to figures/data/:
      - MitoCarta_pathway_annotation.xlsx -- 'ComplexI'..'ComplexV',
        'Translation', 'mtDNA_maintenance', 'mtRNA metabolism',
        'Lipid_metabolism', 'MitoCarta3.0_MitoPathways', by 'Gene names'.
      - mitoproteome_isolated_mito.txt -- plain-text reference gene list
        (735 genes, 'Verified_mitoproteome_ISO' in the curated file).

This script's own hardcoded SULFUR / NEAA / ONE_C_CORE / NUCLEOTIDE_CORE gene sets below
match the curated file's Sulfur_metabolism / NEAA_metabolism / 1C_metabolism /
Nucleotide_metabolism sheets gene-for-gene. A full run against EV_Table_2 and the curated
file reproduces the manuscript's reported n=4 for 1C metabolism and n=6 for Sulfur
metabolism exactly.
"""
import glob
import warnings
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns

warnings.filterwarnings('ignore')

# ── Custom Matplotlib Settings for Publishing ────────────────────────────────
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# ── File paths ───────────────────────────────────────────────────────────────
DATA_DIR = 'data'
OUT_DIR = '.'


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


# ---- EV table ----
ISO_MITO_DATA_PATH = find_file(
    ['EV_Table_2*.xlsx'],
    'isolated-mito TMT solubility/aggregation data')

# ---- not EV tables -- built by verification/build_isolated_mito_annotation_from_curated_lists.py ----
ISO_MITO_ANNOT_PATH = find_file(
    ['MitoCarta_pathway_annotation*.xlsx', 'Data_processed_HS_Soluble_Insoluble*.xlsx',
     'Annotations_MitoCarta*.xlsx', 'Isolated_mito_MitoCarta_annotations*.xlsx'],
    "MitoCarta pathway/complex annotation for the isolated-mito gene set (needs "
    "columns: 'Gene names', 'ComplexI'..'ComplexV', 'Translation', 'mtDNA_maintenance', "
    "'mtRNA metabolism', 'Lipid_metabolism', 'MitoCarta3.0_MitoPathways') "
    "-- not an EV table, needs to be added as a standalone deposited file")
MITOPROTEOME_ISO_FILE = find_file(
    ['mitoproteome_isolated_mito.txt'],
    'isolated-mito reference gene list (plain text) '
    '-- not an EV table, needs to be added as a standalone deposited file')

print('Required inputs found:')
for _label, _path in [('ISO_MITO_DATA_PATH', ISO_MITO_DATA_PATH),
                       ('ISO_MITO_ANNOT_PATH', ISO_MITO_ANNOT_PATH),
                       ('MITOPROTEOME_ISO_FILE', MITOPROTEOME_ISO_FILE)]:
    print(f'  {_label}: {_path}')
print()

# ── Load data ─────────────────────────────────────────────────────────────────
# sheet_name=0 (by position) + column-whitespace strip: prior EV tables in this
# project have shown up with renamed sheets / stray leading spaces in different
# exported copies of the "same" file.
data_df = pd.read_excel(ISO_MITO_DATA_PATH, sheet_name=0)
data_df.columns = data_df.columns.str.strip()

mcm = 'MitoCarta3.0_MitoPathways'
FCm = 'log2 FC TOTAL Sen_Pro'
ANNOT_COLS = ['ComplexI', 'ComplexII', 'ComplexIII', 'ComplexIV', 'ComplexV',
              'Translation', 'mtDNA_maintenance', 'mtRNA metabolism',
              'Lipid_metabolism', mcm]

# ISO_MITO_ANNOT_PATH can have multiple sheets, and the one we need isn't necessarily
# sheet 0. Rather than assume a position or a hardcoded name, scan every sheet and use
# the first one that actually has 'Gene names' plus all the required annotation columns.
_annot_xl = pd.ExcelFile(ISO_MITO_ANNOT_PATH)
annot_df = None
_annot_sheet_used = None
for _sheet in _annot_xl.sheet_names:
    _d = pd.read_excel(ISO_MITO_ANNOT_PATH, sheet_name=_sheet)
    # str(c) guards against sheets with no real header (integer column names), which
    # a junk/summary sheet elsewhere in the same workbook can have.
    _d.columns = [str(c).strip() for c in _d.columns]
    if 'Gene names' in _d.columns and all(c in _d.columns for c in ANNOT_COLS):
        annot_df = _d
        _annot_sheet_used = _sheet
        break

if annot_df is None:
    raise KeyError(
        f"None of the sheets in {ISO_MITO_ANNOT_PATH} ({_annot_xl.sheet_names}) have "
        f"'Gene names' plus all of {ANNOT_COLS}. Check the filename/sheet -- this needs "
        "to be the sheet with the MitoCarta pathway/complex annotation columns, "
        "not the raw solubility/aggregation data sheet."
    )
print(f'Using annotation sheet: {_annot_sheet_used!r} from {ISO_MITO_ANNOT_PATH}')

# Merge the EV-table data with the (still non-EV-table) MitoCarta annotation, by gene.
dfm = data_df.merge(annot_df[['Gene names'] + ANNOT_COLS], on='Gene names', how='left')

with open(MITOPROTEOME_ISO_FILE) as f:
  MITO_ISO = set(f.read().split())

detected_iso = set(dfm['Gene names'].dropna())
ref_genes = MITO_ISO & detected_iso
ref_fc_iso = dfm[dfm['Gene names'].isin(ref_genes)][FCm].dropna().values
ref_median = np.median(ref_fc_iso)
dfm = dfm.copy()
dfm['FC_rel'] = dfm[FCm] - ref_median


# ── Gene sets ─────────────────────────────────────────────────────────────────
def relg(genes):
  return dfm[dfm['Gene names'].isin(genes)]['FC_rel'].dropna().values


def relcol(col):
  return dfm[dfm[col] == 1]['FC_rel'].dropna().values


SULFUR = {'ETHE1', 'GOT2', 'MPST', 'MSRA', 'SQRDL', 'TST'}
NEAA = {
    'GOT1',
    'GOT2',
    'GPT2',
    'PSAT1',
    'PSPH',
    'PHGDH',
    'SHMT2',
    'PYCR1',
    'PYCR2',
    'ALDH18A1',
    'ASS1',
    'ASL',
    'GLS',
    'GLUD1',
    'ASNS',
    'BCAT2',
}
ONE_C_CORE = {'SHMT2', 'MTHFD2', 'MTHFD1L', 'ALDH1L2'}
NUCLEOTIDE_CORE = {
    'AK2',
    'AK3',
    'AK4',
    'DHODH',
    'DTYMK',
    'DUT',
    'FAM173A',
    'HINT1',
    'NME3',
    'NME4',
    'NT5DC2',
    'NUDT5',
    'NUDT9',
    'PAICS',
    'PPA2',
    'SLC25A23',
    'SLC25A24',
    'SLC25A25',
    'SLC25A4',
    'SLC25A5',
    'SLC25A6',
    'SUCLA2',
    'SUCLG1',
    'SUCLG2',
    'TK2',
}
mito_ge_m = (
    set(dfm[dfm['Translation'] == 1]['Gene names'])
    | set(dfm[dfm['mtDNA_maintenance'] == 1]['Gene names'])
    | set(dfm[dfm['mtRNA metabolism'] == 1]['Gene names'])
)
bcaa_m = set(
    dfm[dfm[mcm].str.contains('Branched-chain', na=False)]['Gene names']
)
fao_m = set(dfm[dfm['Lipid_metabolism'] == 1]['Gene names'])
sulf_m = SULFUR & detected_iso
neaa_m = NEAA & detected_iso
onec_m = ONE_C_CORE & detected_iso
nuc_m = NUCLEOTIDE_CORE & detected_iso

SUB_COLOR = {
    'CI': '#fbb4ae',
    'CII': '#b3cde3',
    'CIII': '#ccebc5',
    'CIV': '#decbe4',
    'CV': '#fed9a6',
}

RAW_CATS = [
    (
        'Complex I',
        relcol('ComplexI'),
        set(dfm[dfm['ComplexI'] == 1]['Gene names']),
        SUB_COLOR['CI'],
    ),
    (
        'Complex II',
        relcol('ComplexII'),
        set(dfm[dfm['ComplexII'] == 1]['Gene names']),
        SUB_COLOR['CII'],
    ),
    (
        'Complex III',
        relcol('ComplexIII'),
        set(dfm[dfm['ComplexIII'] == 1]['Gene names']),
        SUB_COLOR['CIII'],
    ),
    (
        'Complex IV',
        relcol('ComplexIV'),
        set(dfm[dfm['ComplexIV'] == 1]['Gene names']),
        SUB_COLOR['CIV'],
    ),
    (
        'Complex V',
        relcol('ComplexV'),
        set(dfm[dfm['ComplexV'] == 1]['Gene names']),
        SUB_COLOR['CV'],
    ),
    (
        'Mito gene expression',
        relg(mito_ge_m),
        mito_ge_m,
        '#b3e2cd',
    ),
    ('NEAA metabolism', relg(neaa_m), neaa_m, '#cbd5e8'),
    ('Fatty acid oxidation', relg(fao_m), fao_m, '#fdcdac'),
    ('BCAA metabolism', relg(bcaa_m), bcaa_m, '#eccc68'),
    ('1C metabolism', relg(onec_m), onec_m, '#f4a261'),
    ('Nucleotide metabolism', relg(nuc_m), nuc_m, '#e76f51'),
    ('Sulfur metabolism', relg(sulf_m), sulf_m, '#f1948a'),
]

background_full = ref_fc_iso - ref_median

# ── 1. COMPUTE ORDINARY ANOVA + BONFERRONI POST-HOC ─────────────────────────
active_groups = [fc_cat for _, fc_cat, _, _ in RAW_CATS if len(fc_cat) >= 5]
active_groups.append(background_full)
f_stat_p1, anova_p_p1 = stats.f_oneway(*active_groups)

num_valid_comparisons = len([fc for _, fc, _, _ in RAW_CATS if len(fc) >= 5])

# ── 2. TRUE BIOLOGICALLY ABSOLUTE AGGREGATION VALUE TRACKS & ANOVA ───────────
SOL_P = ['Corrected Sol Pro 1', 'Corrected Sol Pro 2', 'Corrected Sol Pro 3']
INS_P = ['Corrected Insol Pro 1', 'Corrected Insol Pro 2', 'Corrected Insol Pro 3']
SOL_S = ['Corrected Sol Sen 1', 'Corrected Sol Sen 2', 'Corrected Sol Sen 3']
INS_S = ['Corrected Insol Sen 1', 'Corrected Insol Sen 2', 'Corrected Insol Sen 3']


def ai_logfc_true(genes_set):
  sub = dfm[dfm['Gene names'].isin(genes_set)]
  logfc = []
  for i in range(3):
    sp = sub[SOL_P[i]].sum()
    ip = sub[INS_P[i]].sum()
    ss = sub[SOL_S[i]].sum()
    is_ = sub[INS_S[i]].sum()
    if (ss + is_) > 0 and (sp + ip) > 0 and is_ > 0 and ip > 0:
      logfc.append(np.log2((is_ / (ss + is_)) / (ip / (sp + ip))))
    else:
      logfc.append(np.nan)
  return np.array(logfc)


ag_values_true = {}
for label, _, genes_set, _ in RAW_CATS:
  ag_values_true[label] = ai_logfc_true(genes_set)

# Mitoproteome baseline reference aggregation values
ag_values_true['Mitoproteome'] = ai_logfc_true(ref_genes)
mitoproteome_ag = ag_values_true['Mitoproteome']

# --- ONE-WAY ANOVA + BONFERRONI FOR AGGREGATION ---
ag_active_groups = [
    ag_values_true[label]
    for label, _, genes_set, _ in RAW_CATS
    if len(ag_values_true[label][~np.isnan(ag_values_true[label])]) >= 2
]
ag_active_groups.append(mitoproteome_ag)
_, anova_p_ag = stats.f_oneway(*ag_active_groups)

ag_sig_labels = {}
num_valid_ag = len(ag_active_groups) - 1

for label, _, genes_set, _ in RAW_CATS:
  triplet = ag_values_true[label]
  clean_triplet = triplet[~np.isnan(triplet)]

  if len(clean_triplet) < 2:
    ag_sig_labels[label] = 'nq'
  else:
    # Pairwise t-test of category aggregation vs Mitoproteome background aggregation
    _, p_raw_ag = stats.ttest_ind(
        clean_triplet, mitoproteome_ag, equal_var=True
    )
    p_bonf_ag = min(p_raw_ag * num_valid_ag, 1.0)

    if p_bonf_ag < 0.001:
      ag_sig_labels[label] = '***'
    elif p_bonf_ag < 0.01:
      ag_sig_labels[label] = '**'
    elif p_bonf_ag < 0.05:
      ag_sig_labels[label] = '*'
    else:
      ag_sig_labels[label] = 'ns'

# ── 3. MAP AND CONSOLIDATE ARRAYS FOR CANVAS GRAPHICS ────────────────────────
CATS_FORNASIERO = []

for label, fc_cat, genes, col in RAW_CATS:
  if len(fc_cat) < 5:
    sig_p1 = 'n too small'
  else:
    _, p_raw_p1 = stats.ttest_ind(fc_cat, background_full, equal_var=True)
    p_bonf = min(p_raw_p1 * num_valid_comparisons, 1.0)
    sig_p1 = (
        '***'
        if p_bonf < 0.001
        else '**'
        if p_bonf < 0.01
        else '*'
        if p_bonf < 0.05
        else 'ns'
    )

  CATS_FORNASIERO.append((label, fc_cat, genes, col, sig_p1))

CATS_FORNASIERO.append(
    ('Mitoproteome', background_full, ref_genes, '#9E9E9E', '—')
)
n_cats = len(CATS_FORNASIERO)

all_means = np.array([v.mean() for v in ag_values_true.values()])
all_ses = np.array(
    [v.std(ddof=1) / np.sqrt(3) for v in ag_values_true.values()]
)
AG_XMIN, AG_XMAX = np.nanmin(all_means - all_ses) - 0.28, max(
    np.nanmax(all_means + all_ses) + 0.28, 0.15
)


# ── 4. VISUALIZATION CANVAS GENERATOR RENDER ENGINE ──────────────────────────
def render_canvas(categories_list, title_string):
  ROW_H = 0.42
  sns.set_style('ticks')
  fig = plt.figure(figsize=(10.6, 0.54 * n_cats + 1.05))
  gs = gridspec.GridSpec(
      1,
      2,
      width_ratios=[2.6, 1.0],
      wspace=0.06,
      left=0.27,
      right=0.96,
      top=0.91,
      bottom=0.11,
      figure=fig,
  )
  ax_ab = fig.add_subplot(gs[0])
  ax_ag = fig.add_subplot(gs[1])

  for ax in (ax_ab, ax_ag):
    for yi in range(n_cats):
      if yi % 2 == 0:
        ax.axhspan(yi - 0.5, yi + 0.5, color='#F7F7F7', zorder=0)

  for yi, (label, fc, genes, col, sig) in enumerate(categories_list):
    y = n_cats - 1 - yi
    if len(fc) == 0:
      ax_ab.text(
          0,
          y,
          'ns',
          va='center',
          ha='left',
          fontsize=8.5,
          color='grey',
          fontweight='bold',
      )
      continue
    q1, med, q3 = np.percentile(fc, [25, 50, 75])
    iqr = q3 - q1
    wl = max(fc.min(), q1 - 1.5 * iqr)
    wh = min(fc.max(), q3 + 1.5 * iqr)
    outliers = fc[(fc < wl) | (fc > wh)]

    ax_ab.add_patch(
        FancyBboxPatch(
            (q1, y - ROW_H / 2),
            max(q3 - q1, 0.01),
            ROW_H,
            boxstyle='round,pad=0.008',
            facecolor=col,
            edgecolor='white',
            lw=0.8,
            alpha=0.88,
            zorder=3,
        )
    )
    ax_ab.plot(
        [med, med],
        [y - ROW_H / 2, y + ROW_H / 2],
        color='white',
        lw=1.8,
        zorder=5,
    )
    ax_ab.plot([wl, q1], [y, y], color='black', lw=1.1, zorder=2, alpha=0.75)
    ax_ab.plot([q3, wh], [y, y], color='black', lw=1.1, zorder=2, alpha=0.75)
    for wx in [wl, wh]:
      ax_ab.plot(
          [wx, wx],
          [y - 0.09, y + 0.09],
          color='black',
          lw=1.1,
          zorder=2,
          alpha=0.75,
      )
    if len(outliers):
      ax_ab.scatter(
          outliers,
          [y] * len(outliers),
          color=col,
          s=7,
          alpha=0.5,
          edgecolors='none',
          zorder=4,
      )

    ax_ab.text(
        0.87,
        y,
        sig,
        va='center',
        ha='left',
        fontsize=8.5,
        color=col,
        fontweight='bold',
        transform=ax_ab.get_yaxis_transform(),
    )
    ax_ab.text(
        -1.18, y, f'n={len(fc)}', va='center', fontsize=6.8, color='#999999'
    )

  ax_ab.axvline(0, color='#888888', ls=':', lw=0.9, alpha=0.85)
  ax_ab.set_xlim(-1.25, 1.05)
  ax_ab.set_ylim(-0.6, n_cats - 0.4)
  ax_ab.set_xticks([-1, -0.5, 0, 0.5])
  ax_ab.set_xlabel(
      r'$\log_2$ FC relative to mitoproteome median'
      + '\n(senescent vs proliferative)',
      fontsize=10,
  )
  ax_ab.xaxis.grid(True, ls='--', alpha=0.2, color='#CCCCCC', zorder=0)
  ax_ab.set_title(
      title_string
      + '\n'
      + r'$\it{TMT\ proteomics\ on\ isolated\ mitochondria}$',
      fontsize=11.5,
      fontweight='bold',
      pad=10,
  )
  ax_ab.set_yticks(range(n_cats))
  ax_ab.set_yticklabels(
      [categories_list[n_cats - 1 - i][0] for i in range(n_cats)],
      fontsize=10,
      fontweight='bold',
  )
  sns.despine(ax=ax_ab, top=True, right=True, left=True)
  ax_ab.tick_params(axis='y', length=0, pad=4)
  ax_ab.tick_params(axis='x', labelsize=9.5)

  # Render Aggregation Index Right Panel (ANOVA + Bonferroni)
  for yi, (label, _, _, col, _) in enumerate(categories_list):
    y = n_cats - 1 - yi
    fc_ag = ag_values_true[label]
    fc_clean = fc_ag[~np.isnan(fc_ag)]
    if len(fc_clean) == 0:
      ax_ag.text(
          0,
          y,
          'nq',
          ha='left',
          va='center',
          fontsize=8,
          color='grey',
          fontweight='bold',
      )
      continue
    mean_fc = np.nanmean(fc_ag)
    se_fc = stats.sem(fc_ag, nan_policy='omit') if len(fc_clean) > 1 else 0
    ax_ag.barh(
        y,
        mean_fc,
        height=ROW_H,
        color=col,
        alpha=0.85,
        edgecolor='white',
        lw=0.6,
        zorder=3,
    )
    if se_fc > 0:
      ax_ag.errorbar(
          mean_fc,
          y,
          xerr=se_fc,
          fmt='none',
          color='black',
          elinewidth=1.0,
          capsize=2.2,
          zorder=4,
          alpha=0.85,
      )

    sig_ag = '—' if label == 'Mitoproteome' else ag_sig_labels.get(label, 'ns')
    err_edge = mean_fc + se_fc if mean_fc >= 0 else mean_fc - se_fc
    ax_ag.text(
        err_edge + (0.05 if mean_fc >= 0 else -0.05),
        y,
        sig_ag,
        ha='left' if mean_fc >= 0 else 'right',
        va='center',
        fontsize=8,
        color=col,
        fontweight='bold',
        zorder=6,
    )

  ax_ag.axvline(0, color='#888888', lw=1.0, zorder=2)
  ax_ag.set_xlim(AG_XMIN, AG_XMAX)
  ax_ag.set_ylim(-0.6, n_cats - 0.4)
  ax_ag.set_yticks([])
  ax_ag.set_xlabel(
      r'Aggregation index $\log_2$ FC'
      + '\n(Sen vs Prolif)\n$n=3$ biological replicates',
      fontsize=9.5,
  )
  ax_ag.xaxis.grid(True, ls='--', alpha=0.2, color='#CCCCCC', zorder=0)
  ax_ag.set_title('Aggregation', fontsize=11.5, fontweight='bold', pad=10)
  sns.despine(ax=ax_ag, top=True, right=True, left=True)
  ax_ag.tick_params(axis='x', labelsize=9)
  return fig


# ── 5. EXECUTE GRAPHICS PIPELINE & EXPORT ─────────────────────────────────────
fig_fornasiero = render_canvas(
    CATS_FORNASIERO,
    'OXPHOS Pathways vs Mitoproteome (Ordinary ANOVA + Bonferroni)',
)

# Export SVG, PNG, and PDF formats
fig_fornasiero.savefig(
    f'{OUT_DIR}/Figure2B_OrdinaryANOVA_Bonferroni.svg',
    bbox_inches='tight',
)
fig_fornasiero.savefig(
    f'{OUT_DIR}/Figure2B_OrdinaryANOVA_Bonferroni.png',
    dpi=300,
    bbox_inches='tight',
)
fig_fornasiero.savefig(
    f'{OUT_DIR}/Figure2B_OrdinaryANOVA_Bonferroni.pdf',
    bbox_inches='tight',
)

plt.close('all')
print(
    'Done. Generated single Ordinary ANOVA + Bonferroni figure'
    ' (Figure2B_OrdinaryANOVA_Bonferroni in SVG/PNG/PDF).'
)