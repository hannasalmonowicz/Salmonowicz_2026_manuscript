#!/usr/bin/env python3
"""
Fig 1G -- 30-model dot-matrix, MRC5 column from EV_Table_8 (Sen CTRL vs
Prolif CTRL; CAM and IMT are separate processing arms in that file, CTRL is
the arm plotted here, as the direct analog of "senescent vs proliferating"
used for every other model).

Inputs:

EV tables:
      EV_Table_8*.xlsx   -- MRC5 proteomics
      EV_Table_10*.xlsx  -- all 14 SenCat lines
      EV_Table_12*.xlsx  -- Payea et al. 2024 IMR90 etoposide dataset

Pathway gene lists (OXPHOS subunits, fatty acid oxidation, BCAA, NEAA, sulfur
metabolism, 1C/folate) and MITO_GE_GENES (gene symbols annotated to the MitoCarta3.0
"Translation / mtDNA maintenance / mtRNA metabolism" pathways, used for the
Mito_gene_expression category) are fixed reference sets, hardcoded directly below
rather than read from a separate file.
"""
import glob
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
N_PERM = 2000
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
EV_TABLE_8_PATH = find_file(['EV_Table_8*.xlsx'], 'EV Table 8')
ANERILLAS_PATH = find_file(['EV_Table_10*.xlsx'], '30-model TIS meta-analysis (Anerillas)')
PAYEA_PATH = find_file(
    ['EV_Table_12*.xlsx'],
    'EV Table 12 / Payea et al. 2024 IMR90 etoposide dataset')

print('All required inputs found:')
for _label, _path in [('EV_TABLE_8_PATH', EV_TABLE_8_PATH), ('ANERILLAS_PATH', ANERILLAS_PATH),
                       ('PAYEA_PATH', PAYEA_PATH)]:
    print(f'  {_label}: {_path}')
print()

# Hardcoded (curated, cross-checked against MitoCarta3.0; see
# no external file needed.
MITO_GE_GENES = {
    'AARS2', 'ALKBH1', 'ANGEL2', 'APEX1', 'ATAD3A', 'ATAD3B', 'AURKAIP1', 'CARS2', 'CDK5RAP1',
    'CHCHD1', 'COA3', 'COX14', 'DAP3', 'DARS2', 'DDX28', 'DHX30', 'DNA2', 'DUS2', 'EARS2',
    'ELAC2', 'ENDOG', 'ERAL1', 'EXD2', 'EXOG', 'FARS2', 'FASTK', 'FASTKD1', 'FASTKD2',
    'FASTKD3', 'FASTKD5', 'GADD45GIP1', 'GARS1', 'GATB', 'GATC', 'GFM1', 'GFM2', 'GRSF1',
    'GTPBP10', 'GTPBP3', 'GUF1', 'HARS2', 'HEMK1', 'HSD17B10', 'IARS2', 'KARS1', 'KGD4',
    'LACTB2', 'LARS2', 'LIG3', 'LRPPRC', 'MALSU1', 'MARS2', 'METAP1D', 'METTL15', 'METTL17',
    'METTL5', 'METTL8', 'MGME1', 'MIEF1', 'MPV17L2', 'MRM1', 'MRM2', 'MRM3', 'MRPL1', 'MRPL10',
    'MRPL11', 'MRPL12', 'MRPL13', 'MRPL14', 'MRPL15', 'MRPL16', 'MRPL17', 'MRPL18', 'MRPL19',
    'MRPL2', 'MRPL20', 'MRPL21', 'MRPL22', 'MRPL23', 'MRPL24', 'MRPL27', 'MRPL28', 'MRPL3',
    'MRPL30', 'MRPL32', 'MRPL33', 'MRPL34', 'MRPL35', 'MRPL36', 'MRPL37', 'MRPL38', 'MRPL39',
    'MRPL4', 'MRPL40', 'MRPL41', 'MRPL42', 'MRPL43', 'MRPL44', 'MRPL46', 'MRPL47', 'MRPL48',
    'MRPL49', 'MRPL50', 'MRPL51', 'MRPL52', 'MRPL53', 'MRPL54', 'MRPL55', 'MRPL57', 'MRPL58',
    'MRPL9', 'MRPS10', 'MRPS11', 'MRPS12', 'MRPS14', 'MRPS15', 'MRPS16', 'MRPS17', 'MRPS18A',
    'MRPS18B', 'MRPS18C', 'MRPS2', 'MRPS21', 'MRPS22', 'MRPS23', 'MRPS24', 'MRPS25', 'MRPS26',
    'MRPS27', 'MRPS28', 'MRPS30', 'MRPS31', 'MRPS33', 'MRPS34', 'MRPS35', 'MRPS5', 'MRPS6',
    'MRPS7', 'MRPS9', 'MRRF', 'MTERF3', 'MTERF4', 'MTFMT', 'MTG1', 'MTG2', 'MTIF2', 'MTIF3',
    'MTO1', 'MTPAP', 'MTRES1', 'MTRF1', 'MTRF1L', 'MUTYH', 'NARS2', 'NGRN', 'NOA1', 'NSUN2',
    'NSUN4', 'OGG1', 'OSGEPL1', 'OXA1L', 'PARS2', 'PDE12', 'PDF', 'PIF1', 'PNPT1', 'POLB',
    'POLDIP2', 'POLG', 'POLG2', 'POLQ', 'POLRMT', 'PPA2', 'PRORP', 'PTCD1', 'PTCD2', 'PTCD3',
    'PUS1', 'PUSL1', 'QRSL1', 'QTRT1', 'RARS2', 'RBFA', 'RCC1L', 'RECQL4', 'REXO2', 'RMND1',
    'RNASEH1', 'RPUSD3', 'RPUSD4', 'SARS2', 'SLIRP', 'SSBP1', 'SUPV3L1', 'TACO1', 'TARS2',
    'TBRG4', 'TEFM', 'TFAM', 'TFB1M', 'TFB2M', 'THG1L', 'TIMM21', 'TOP3A', 'TRIT1', 'TRMT1',
    'TRMT10C', 'TRMT2B', 'TRMT5', 'TRMT61B', 'TRMU', 'TRNT1', 'TRUB2', 'TSFM', 'TUFM', 'TWNK',
    'UNG', 'VARS2', 'WARS2', 'YARS2', 'YBEY', 'YRDC'
}  # n=222

# Pathway gene lists (fixed reference sets, not derived from any input file).
gene_lists = {
    'OXPHOS_subunits': [
        'ATP5F1A', 'ATP5F1B', 'ATP5F1C', 'ATP5F1D', 'ATP5PB', 'ATP5PD', 'ATP5ME', 'ATP5PF',
        'ATP5MF', 'ATP5MG', 'ATP5PO', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C',
        'COX7A2', 'COX7A2L', 'COX7C', 'CYC1', 'CYCS', 'HCCS', 'MT-ATP6', 'MT-ATP8', 'MT-CO1',
        'MT-CO2', 'MT-ND1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3',
        'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1',
        'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7',
        'NDUFB8', 'NDUFB9', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS7',
        'NDUFS8', 'NDUFV1', 'NDUFV2', 'SDHA', 'SDHB', 'UQCR11', 'UQCRB', 'UQCRC1', 'UQCRC2',
        'UQCRFS1', 'UQCRH', 'UQCRQ',
    ],
    'Fatty_acid_oxidation': [
        'ACAA1', 'ACAA2', 'ACACA', 'ACAD10', 'ACAD11', 'ACADM', 'ACADSB', 'ACADVL', 'ACAT1',
        'ACOT13', 'ACOT7', 'ACOT9', 'ACSF2', 'ACSF3', 'ACSL1', 'AGK', 'CPT1A', 'CPT2', 'CRAT',
        'CROT', 'CYB5R3', 'DBI', 'DECR1', 'DHRS1', 'ECH1', 'ECHDC1', 'ECHS1', 'ECI1', 'ECI2',
        'ETFA', 'ETFB', 'ETFDH', 'FASN', 'FDPS', 'FDXR', 'GCSH', 'HADH', 'HADHA', 'HADHB',
        'HINT2', 'HSD17B10', 'HSD17B4', 'IDI1', 'LACTB', 'LYPLA1', 'MCEE', 'MGST3', 'MUT',
        'NDUFAB1', 'OSBPL1A', 'PCCB', 'PLSCR3', 'PRDX6', 'PTGES2', 'PTPMT1', 'SCP2', 'SLC25A1',
        'SLC25A20', 'SPTLC2', 'TAMM41', 'TSPO',
    ],
    'BCAA_metabolism': [
        'ACADSB', 'ACAT1', 'ALDH6A1', 'BCAT2', 'BCKDHA', 'DBT', 'DLD', 'ECHS1', 'ETFA', 'ETFB',
        'ETFDH', 'HADHA', 'HIBADH', 'HIBCH', 'HMGCL', 'HSD17B10', 'IVD', 'MCCC1', 'MCCC2',
    ],
    'NEAA_metabolism': [
        'ALDH18A1', 'ASL', 'ASNS', 'ASS1', 'BCAT2', 'GLS', 'GOT1', 'GOT2', 'PHGDH', 'PSAT1',
        'PSPH', 'PYCR1', 'PYCR2', 'SHMT2',
    ],
    'Sulfur_metabolism': ['ETHE1', 'GOT2', 'MPST', 'MSRA', 'SQOR', 'TST'],
    'OneC_folate_metabolism': [
        'ALDH1L1', 'MTHFD1', 'MTHFR', 'SHMT1', 'DHFR', 'TYMS', 'MTHFD1L', 'MTHFD2', 'ALDH1L2',
        'SHMT2', 'GART', 'ATIC',
    ],
}
print('1C-folate category has', len(gene_lists['OneC_folate_metabolism']), 'genes')


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

# ---------------------------------------------------------- MRC5 (EV_Table_8)
# sheet_name=0 (by position) + column-whitespace strip: prior EV tables in this project
# have shown up with renamed sheets / stray leading spaces in different exported copies.
cam = pd.read_excel(EV_TABLE_8_PATH, sheet_name=0)
cam.columns = cam.columns.str.strip()
cam = cam.dropna(subset=['Gene names']).drop_duplicates(subset=['Gene names'], keep='first')
cam['Gene names'] = cam['Gene names'].map(rn)
cam = cam.set_index('Gene names')
cam_fc = cam['Log2FC Sen CTRL vs Prolif CTRL']
MRC5_LABEL = 'MRC5'
recs, bg = records_for_series(cam_fc, MRC5_LABEL, 'Reference')
model_background[MRC5_LABEL] = bg
cam_bg_mean = cam_fc.dropna().mean()
cam_corrected = cam_fc.dropna() - cam_bg_mean
ge_genes_cam = set(g for g in cam_fc.dropna().index if g in MITO_GE_GENES)
add_mito_ge(recs, cam_corrected, MRC5_LABEL, 'Reference', ge_genes_cam)
all_records += recs
print(MRC5_LABEL, 'full-proteome background size:', len(bg), '| Mito_gene_expression genes:', len(ge_genes_cam))

# ---------------------------------------------------------- IMR90 (Payea 2024)
payea = pd.read_excel(PAYEA_PATH, sheet_name=0)
payea.columns = payea.columns.str.strip()
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
ge_genes_payea = set(g for g in payea_fc.index if g in MITO_GE_GENES)
add_mito_ge(recs, payea_corrected, PAYEA_LABEL, 'Reference', ge_genes_payea)
all_records += recs

# ---------------------------------------------------------- Anerillas SenCat
CELLTYPE_SHEETS = {'WI38': 'WI38', 'BJ': 'BJ', 'HSAEC': 'HSAEC', 'HEKn': 'HEKn', 'HCAEC': 'HCAEC',
                    'HUVEC': 'HUVEC', 'BMMSC': 'BMMSC', 'HVSMC': 'HVSMC', 'HSKM': 'HSKM', 'PBMC': 'PBMC',
                    'PreAdipo': 'PreAdipo', 'NHO': 'NHO', 'NHA': 'NHA', 'HEMn': 'HEM'}

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
        ge_genes = set(g for g in diff.dropna().index if g in MITO_GE_GENES)
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

# ================================================================== plotting
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


def build_figure(mrc5_label, mrc5_tag, out_prefix):
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

    fig.text(0.01, 0.99, 'Direction, magnitude, and significance across 30 models', fontsize=14,
              fontweight='bold', color=INK, ha='left', va='top')

    fig.subplots_adjust(left=0.125, right=0.90, top=0.865, bottom=0.24)
    fig.savefig(f'{out_prefix}.png', dpi=300, facecolor=SURFACE)
    fig.savefig(f'{out_prefix}.svg', facecolor=SURFACE)
    plt.close(fig)
    print('saved', f'{out_prefix}.png/svg')


build_figure(MRC5_LABEL, 'DOX', 'Generality_matrix')
