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
# Pathway gene lists (curated whole-cell MitoCarta subsets, NEAA, nucleotide/1C,
# methionine metabolism, etc.) are fixed reference sets and are written directly
# below rather than read from a separate file.
# =====================================================================
import glob
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

# ── Pathway gene lists (fixed reference sets, not derived from any input file) ─
ALL_PROTEINS = set(df['Gene names'].dropna())

MITO_REF = {
    'AARS2', 'AASS', 'ABAT', 'ABCD3', 'ABHD10', 'ABHD11', 'ACAA1', 'ACAA2', 'ACACA', 'ACAD10',
    'ACAD11', 'ACAD9', 'ACADM', 'ACADSB', 'ACADVL', 'ACAT1', 'ACLY', 'ACO2', 'ACOT13', 'ACOT7',
    'ACOT9', 'ACSF2', 'ACSF3', 'ACSL1', 'AFG3L2', 'AGK', 'AHCYL1', 'AIFM1', 'AK2', 'AK3', 'AKAP1',
    'AKR7A2', 'ALDH18A1', 'ALDH1B1', 'ALDH1L2', 'ALDH2', 'ALDH4A1', 'ALDH6A1', 'ALDH7A1',
    'ALDH9A1', 'APEX1', 'APOA1BP', 'APOO', 'APOOL', 'ARL2', 'ARMCX1', 'ARMCX3', 'ATAD1', 'ATAD3A',
    'ATAD3B', 'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5F1', 'ATP5H', 'ATP5I', 'ATP5J', 'ATP5J2',
    'ATP5L', 'ATP5O', 'ATPAF1', 'AURKAIP1', 'BAX', 'BCAT2', 'BCKDHA', 'BCL2L1', 'BCL2L13', 'BCS1L',
    'BID', 'BPHL', 'C1QBP', 'C20orf24', 'C6orf203', 'CARKD', 'CARS2', 'CASP3', 'CASP8', 'CAT',
    'CBR3', 'CCBL2', 'CCDC51', 'CHCHD2;CHCHD2P9', 'CHCHD3', 'CHCHD6', 'CISD1', 'CLPP', 'CLPX',
    'COA3', 'COA4', 'COA6', 'COASY', 'COMT', 'COQ9', 'COX20', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1',
    'COX6B1', 'COX6C', 'COX7A2', 'COX7A2L', 'COX7C', 'CPOX', 'CPT1A', 'CPT2', 'CRAT', 'CROT', 'CS',
    'CYB5B', 'CYB5R3', 'CYC1', 'CYCS', 'DAP3', 'DARS2', 'DBI', 'DBT', 'DCXR', 'DDX28', 'DECR1',
    'DHRS1', 'DHX30', 'DIABLO', 'DLAT', 'DLD', 'DLST', 'DNAJA3', 'DNAJC11', 'DNAJC19', 'DNM1L',
    'DTYMK', 'DUS2', 'DUT', 'EARS2', 'ECH1', 'ECHDC1', 'ECHS1', 'ECI1', 'ECI2', 'ECSIT', 'ELAC2',
    'ENDOG', 'ETFA', 'ETFB', 'ETFDH', 'ETHE1', 'EXOG', 'FAHD1', 'FAM213A', 'FASN', 'FASTKD2',
    'FDPS', 'FDX1L', 'FDXR', 'FECH', 'FH', 'FIS1', 'FKBP10', 'FKBP8', 'FLAD1', 'FOXRED1', 'FTH1',
    'FXN', 'GADD45GIP1', 'GARS', 'GBAS', 'GCDH', 'GCSH', 'GFM1', 'GLRX5', 'GLS', 'GLUD1;GLUD2',
    'GOT2', 'GPD2', 'GPT2;GPT', 'GPX1', 'GRHPR', 'GRPEL1', 'GRSF1', 'GSR', 'GSTK1', 'GSTZ1', 'GUK1',
    'HADH', 'HADHA', 'HADHB', 'HAGH', 'HARS2', 'HCCS', 'HIBADH', 'HIBCH', 'HIGD2A', 'HINT1',
    'HINT2', 'HMGCL', 'HSD17B10', 'HSD17B4', 'HSPA9', 'HSPD1', 'HSPE1', 'HTRA2', 'IARS2', 'IDE',
    'IDH2', 'IDH3A', 'IDH3B', 'IDH3G', 'IDI1', 'IMMT', 'ISCA2', 'ISCU', 'IVD', 'KARS', 'LACTB',
    'LACTB2', 'LAP3', 'LARS2', 'LDHB', 'LETM1', 'LIG3', 'LONP1', 'LRPPRC', 'LYPLA1', 'LYRM7',
    'MAVS', 'MCCC1', 'MCCC2', 'MCEE', 'MCU', 'MDH2', 'ME2', 'MFF', 'MFN1', 'MFN2', 'MGST1',
    'MGST3', 'MICU1', 'MICU2', 'MIPEP', 'MMAB', 'MPST', 'MPV17', 'MRPL1', 'MRPL10', 'MRPL11',
    'MRPL12', 'MRPL13', 'MRPL14', 'MRPL15', 'MRPL16', 'MRPL17', 'MRPL18', 'MRPL19', 'MRPL2',
    'MRPL21', 'MRPL22', 'MRPL23', 'MRPL24', 'MRPL27', 'MRPL28', 'MRPL3', 'MRPL37', 'MRPL38',
    'MRPL39', 'MRPL4', 'MRPL40', 'MRPL41', 'MRPL43', 'MRPL44', 'MRPL46', 'MRPL47', 'MRPL48',
    'MRPL49', 'MRPL50', 'MRPL52', 'MRPL54', 'MRPL57', 'MRPL9', 'MRPS11', 'MRPS12', 'MRPS14',
    'MRPS16', 'MRPS17', 'MRPS18A', 'MRPS18B', 'MRPS2', 'MRPS22', 'MRPS23', 'MRPS25', 'MRPS26',
    'MRPS27', 'MRPS30', 'MRPS31', 'MRPS34', 'MRPS35', 'MRPS36', 'MRPS5', 'MRPS7', 'MRPS9', 'MRRF',
    'MSRA', 'MSRB3', 'MT-ATP6', 'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-ND1', 'MTCH1', 'MTCH2',
    'MTHFD1L', 'MTHFD2', 'MTX1', 'MTX2', 'MUT', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13',
    'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1',
    'NDUFAF2', 'NDUFAF3', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3', 'NDUFB4', 'NDUFB5',
    'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5',
    'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NFS1', 'NFU1', 'NIPSNAP1', 'NIT1', 'NIT2', 'NLN',
    'NME3', 'NNT', 'NRD1', 'NSUN2', 'NT5DC2', 'NUDT19', 'NUDT2', 'NUDT5', 'NUDT9', 'OAT', 'OGDH',
    'OPA1', 'OSBPL1A', 'OXCT1', 'OXR1', 'PAICS', 'PAM16', 'PARK7', 'PC', 'PCCB', 'PCK2', 'PDE12',
    'PDHA1', 'PDHB', 'PDHX', 'PDPR', 'PGAM5', 'PHB', 'PHB2', 'PITRM1', 'PLSCR3', 'PMPCA', 'PMPCB',
    'PNPO', 'PNPT1', 'POLDIP2', 'PPA2', 'PPIF', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6',
    'PRKACA', 'PROSC', 'PTCD3', 'PTGES2', 'PTPMT1', 'PUS1', 'PYCR1', 'PYCR2', 'QDPR', 'QIL1',
    'QTRT1', 'RARS2', 'RDH13', 'REXO2', 'RFK', 'RHOT1', 'RHOT2', 'RMDN3', 'ROMO1', 'RPIA',
    'SAMM50', 'SCO1', 'SCP2', 'SDHA', 'SDHB', 'SDSL', 'SFXN1', 'SFXN3', 'SHMT2', 'SLC25A1',
    'SLC25A11', 'SLC25A12', 'SLC25A13', 'SLC25A20', 'SLC25A22', 'SLC25A24', 'SLC25A29', 'SLC25A3',
    'SLC25A4', 'SLC25A5', 'SLC25A6', 'SLIRP', 'SMIM20', 'SNAP29', 'SOD1', 'SOD2', 'SPG7', 'SPR',
    'SPTLC2', 'SQRDL', 'SSBP1', 'STOML2', 'SUCLA2', 'SUCLG1', 'SUCLG2', 'SYNJ2BP', 'TACO1',
    'TAMM41', 'TBRG4', 'TFAM', 'TFB2M', 'TIMM13', 'TIMM17A', 'TIMM44', 'TIMM50', 'TIMM8A', 'TIMM9',
    'TMEM11', 'TMEM126A', 'TOMM20', 'TOMM22', 'TOMM34', 'TOMM40', 'TOMM5', 'TOMM6', 'TOMM70A',
    'TRAP1', 'TRMT1', 'TRMT10C', 'TRNT1', 'TSFM', 'TSPO', 'TST', 'TUFM', 'TXN2', 'TXNRD1',
    'TXNRD2', 'UQCR11', 'UQCRB', 'UQCRC1', 'UQCRC2', 'UQCRFS1;UQCRFS1P1', 'UQCRH', 'UQCRQ',
    'VDAC1', 'VDAC2', 'VDAC3', 'XPNPEP3', 'YARS2', 'YME1L1', 'YRDC',
}  # n=500

OXPHOS_GENES = {
    'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5F1', 'ATP5H', 'ATP5I', 'ATP5J', 'ATP5J2', 'ATP5L',
    'ATP5O', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A2', 'COX7A2L',
    'COX7C', 'CYC1', 'CYCS', 'HCCS', 'MT-ATP6', 'MT-ATP8', 'MT-CO1', 'MT-CO2', 'MT-ND1',
    'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA4', 'NDUFA5', 'NDUFA6',
    'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB2', 'NDUFB3',
    'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFS1', 'NDUFS2', 'NDUFS3',
    'NDUFS4', 'NDUFS5', 'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'SDHA', 'SDHB', 'UQCR11', 'UQCRB',
    'UQCRC1', 'UQCRC2', 'UQCRFS1;UQCRFS1P1', 'UQCRH', 'UQCRQ',
}  # n=70

NEAA = {
    'ALDH18A1', 'ASL', 'ASNS', 'ASS1', 'BCAT2', 'GLS', 'GLUD1', 'GOT1', 'GOT2', 'GPT2', 'PHGDH',
    'PSAT1', 'PSPH', 'PYCR1', 'PYCR2', 'SHMT2',
}  # n=16 in the curated list; n=14 matched to the detected whole-cell proteome

FAO = {
    'ACAA1', 'ACAA2', 'ACACA', 'ACAD10', 'ACAD11', 'ACADM', 'ACADSB', 'ACADVL', 'ACAT1', 'ACOT13',
    'ACOT7', 'ACOT9', 'ACSF2', 'ACSF3', 'ACSL1', 'AGK', 'CPT1A', 'CPT2', 'CRAT', 'CROT', 'CYB5R3',
    'DBI', 'DECR1', 'DHRS1', 'ECH1', 'ECHDC1', 'ECHS1', 'ECI1', 'ECI2', 'ETFA', 'ETFB', 'ETFDH',
    'FASN', 'FDPS', 'FDXR', 'GCSH', 'HADH', 'HADHA', 'HADHB', 'HINT2', 'HSD17B10', 'HSD17B4',
    'IDI1', 'LACTB', 'LYPLA1', 'MCEE', 'MGST3', 'MUT', 'NDUFAB1', 'OSBPL1A', 'PCCB', 'PLSCR3',
    'PRDX6', 'PTGES2', 'PTPMT1', 'SCP2', 'SLC25A1', 'SLC25A20', 'SPTLC2', 'TAMM41', 'TSPO',
}  # n=61

BCAA = {
    'ACADSB', 'ACAT1', 'ALDH6A1', 'BCAT2', 'BCKDHA', 'DBT', 'DLD', 'ECHS1', 'ETFA', 'ETFB',
    'ETFDH', 'HADHA', 'HIBADH', 'HIBCH', 'HMGCL', 'HSD17B10', 'IVD', 'MCCC1', 'MCCC2',
}  # n=19

SULFUR = {'ETHE1', 'GOT2', 'MPST', 'MSRA', 'SQRDL', 'TST'}  # n=6

NUCLEOTIDE_METABOLISM = {
    'AK2', 'AK3', 'AK4', 'DHODH', 'DTYMK', 'DUT', 'FAM173A', 'HINT1', 'NME3', 'NME4', 'NT5DC2',
    'NUDT5', 'NUDT9', 'PAICS', 'PPA2', 'SLC25A23', 'SLC25A24', 'SLC25A25', 'SLC25A4', 'SLC25A5',
    'SLC25A6', 'SUCLA2', 'SUCLG1', 'SUCLG2', 'TK2',
}  # n=25
ONE_C_METABOLISM = {'ALDH1L2', 'MTHFD1L', 'MTHFD2', 'SHMT2'}  # n=4
MITO_NUC_1C = NUCLEOTIDE_METABOLISM | ONE_C_METABOLISM  # n=22 (Nucleotide_metabolism | 1C_metabolism)

# GO:BP methionine metabolism -- hardcoded, matches the published figure. (The
# curated file's own 'Methionine_metabolism' sheet has only 5 genes; this 12-gene
# set is the correct/verified one.)
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
