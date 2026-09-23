"""
Fig 2B -- isolated-mito abundance (left panel) + aggregation (right panel), Ordinary
ANOVA + Bonferroni. This is the exact same figure as before (both panels, one figure);
"abundance" and "aggregation" are the two halves of Figure2B_OrdinaryANOVA_Bonferroni,
not two separate figures.

Input -- one EV table only:
      EV_Table_2*.xlsx -- 'Gene names', per-replicate Corrected Sol/Insol Pro/Sen,
      and 'log2 FC TOTAL Sen_Pro'.

Every gene-category definition below (Complex I-V, mito gene expression, fatty acid
oxidation, BCAA/sulfur/NEAA/1C/nucleotide metabolism, and the 735-gene mitoproteome
reference background) is hardcoded as a Python set literal, extracted from a manually
curated categorization (the isolated-mito analog of the one already used by the
trajectory script) and cross-checked against the MitoCarta3.0 pathway annotation. A full run against
EV_Table_2 reproduces the manuscript's reported n=4 for 1C metabolism and n=6 for
Sulfur metabolism exactly.
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


# ── Hardcoded gene sets (curated, cross-checked against MitoCarta3.0) --
# no external file needed ──────────────────────────────────────────────
MITOPROTEOME_ISO = {
    'AADAT', 'AARS2', 'AASS', 'ABAT', 'ABCB10', 'ABCB6', 'ABCB7', 'ABCB8', 'ABCD1', 'ABCD3',
    'ABHD10', 'ABHD11', 'ACAA1', 'ACAA2', 'ACACA', 'ACAD10', 'ACAD11', 'ACAD8', 'ACAD9',
    'ACADM', 'ACADS', 'ACADSB', 'ACADVL', 'ACAT1', 'ACLY', 'ACO2', 'ACOT13', 'ACOT7', 'ACOT9',
    'ACP6', 'ACSF2', 'ACSF3', 'ACSL1', 'ACSS3', 'AFG3L2', 'AGK', 'AGPAT5', 'AHCYL1', 'AIFM1',
    'AIFM2', 'AK2', 'AK3', 'AK4', 'AKAP1', 'AKAP10', 'AKR7A2', 'ALDH18A1', 'ALDH1B1',
    'ALDH1L1', 'ALDH1L2', 'ALDH2', 'ALDH3A2', 'ALDH4A1', 'ALDH5A1', 'ALDH6A1', 'ALDH7A1',
    'ALDH9A1', 'ALKBH1', 'APEX1', 'APOA1BP', 'APOO', 'APOOL', 'ARL2;ARL2-SNX15', 'ARMC10',
    'ARMCX1', 'ARMCX3', 'ATAD1', 'ATAD3A', 'ATAD3B', 'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D',
    'ATP5E;ATP5EP2', 'ATP5F1', 'ATP5G1;ATP5G3;ATP5G2', 'ATP5H', 'ATP5I', 'ATP5J', 'ATP5J2',
    'ATP5L', 'ATP5O', 'ATP5S', 'ATP5SL', 'ATPAF1', 'ATPAF2', 'ATPIF1', 'AURKAIP1', 'BAD',
    'BAK1', 'BAX', 'BCAT2', 'BCKDHA', 'BCKDHB', 'BCKDK', 'BCL2L1', 'BCL2L13', 'BCS1L', 'BID',
    'BLOC1S1', 'BNIP3', 'BNIP3L', 'BOLA3', 'BPHL', 'C19orf52', 'C1QBP', 'C20orf24', 'C2orf47',
    'C6orf203', 'CARKD', 'CARS2', 'CAT', 'CBR4', 'CCBL2', 'CCDC109B', 'CCDC51', 'CDK5RAP1',
    'CHCHD1', 'CHCHD3', 'CHCHD6', 'CISD1', 'CISD3', 'CLPB', 'CLPP', 'CLPX', 'CMC1', 'COA1',
    'COA3', 'COA4', 'COA6', 'COA7', 'COASY', 'COMT', 'COQ10B', 'COQ3', 'COQ5', 'COQ6', 'COQ7',
    'COQ9', 'COX11', 'COX15', 'COX18', 'COX19', 'COX20', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1',
    'COX6B1', 'COX6C', 'COX7A1', 'COX7A2', 'COX7A2L', 'COX7B', 'COX7C', 'CPOX', 'CPT1A',
    'CPT2', 'CRAT', 'CROT', 'CS', 'CYB5B', 'CYB5R3', 'CYC1', 'CYCS', 'CYP27A1', 'D2HGDH',
    'DAP3', 'DARS2', 'DBI', 'DBT', 'DCXR', 'DDX28', 'DECR1', 'DHODH', 'DHRS1', 'DHRS4',
    'DHTKD1', 'DHX30', 'DIABLO', 'DLAT', 'DLD', 'DLST', 'DNAJA3', 'DNAJC11', 'DNAJC15', 'DNLZ',
    'DNM1L', 'DTYMK', 'DUS2', 'DUT', 'EARS2', 'ECH1', 'ECHDC1', 'ECHS1', 'ECI1;DCI', 'ECI2',
    'ECSIT', 'EHHADH', 'ELAC2', 'ENDOG', 'ERAL1', 'ETFA', 'ETFB', 'ETFDH', 'ETHE1', 'EXD2',
    'EXOG', 'FAHD1', 'FAM173A', 'FAM210B', 'FAM213A', 'FARS2', 'FASN', 'FASTKD1', 'FASTKD2',
    'FASTKD5', 'FDPS', 'FDX1', 'FDX1L', 'FDXR', 'FECH', 'FH', 'FIS1', 'FKBP10', 'FKBP8',
    'FLAD1', 'FOXRED1', 'FTH1', 'FUNDC1', 'FUNDC2', 'FXN', 'GADD45GIP1', 'GARS', 'GATB;PET112',
    'GBAS', 'GCAT', 'GCDH', 'GCSH', 'GFER', 'GFM1', 'GFM2', 'GHITM', 'GLRX5', 'GLS',
    'GLUD1;GLUD2', 'GOT2', 'GPAM', 'GPD2', 'GPT2', 'GPX1', 'GPX4', 'GRHPR', 'GRPEL1', 'GRSF1',
    'GSR', 'GSTK1', 'GTPBP10', 'GTPBP3', 'GUF1', 'HADH', 'HADHA', 'HADHB', 'HAGH', 'HARS2',
    'HCCS', 'HEMK1', 'HIBADH', 'HIBCH', 'HIGD1A', 'HIGD2A', 'HINT1', 'HINT2', 'HMGCL',
    'HSD17B10', 'HSD17B4', 'HSPA9', 'HSPD1', 'HSPE1', 'HTRA2', 'IARS2', 'IBA57', 'ICT1', 'IDE',
    'IDH2', 'IDH3A', 'IDH3B', 'IDH3G', 'IDI1', 'IMMT', 'ISCA2', 'ISCU', 'IVD', 'KARS',
    'KIAA0391', 'L2HGDH', 'LACE1', 'LACTB', 'LACTB2', 'LAP3', 'LARS2', 'LDHB', 'LETM1',
    'LETMD1', 'LIG3', 'LRPPRC', 'LYPLA1', 'LYPLAL1', 'LYRM4', 'LYRM7', 'MACROD1', 'MALSU1',
    'MAOA', 'MAOB', 'MARC2', 'MARCH5', 'MARS2', 'MAVS', 'MCAT', 'MCCC1', 'MCCC2', 'MCEE',
    'MCU', 'MDH2', 'ME2', 'ME3', 'MECR', 'METTL15', 'METTL17', 'MFF', 'MFN1', 'MFN2', 'MGARP',
    'MGME1', 'MGST1', 'MGST3', 'MICU1', 'MICU2', 'MICU3', 'MIEF1', 'MIPEP', 'MLYCD', 'MMAA',
    'MMAB', 'MOCS1', 'MP68;C14orf2', 'MPC2', 'MPST', 'MPV17', 'MPV17L2', 'MRM1', 'MRPL1',
    'MRPL10', 'MRPL11', 'MRPL12', 'MRPL13', 'MRPL14', 'MRPL15', 'MRPL16', 'MRPL17', 'MRPL18',
    'MRPL19', 'MRPL2', 'MRPL20', 'MRPL21', 'MRPL22', 'MRPL23', 'MRPL24', 'MRPL27', 'MRPL28',
    'MRPL3', 'MRPL30', 'MRPL32', 'MRPL33', 'MRPL34', 'MRPL35', 'MRPL37', 'MRPL38', 'MRPL39',
    'MRPL4', 'MRPL40', 'MRPL41', 'MRPL42', 'MRPL43', 'MRPL44', 'MRPL46', 'MRPL47', 'MRPL48',
    'MRPL49', 'MRPL50', 'MRPL51', 'MRPL53', 'MRPL54', 'MRPL55', 'MRPL57', 'MRPL9', 'MRPS10',
    'MRPS11', 'MRPS12', 'MRPS14', 'MRPS15', 'MRPS16', 'MRPS17;hCG_1984214', 'MRPS18A',
    'MRPS18B', 'MRPS18C', 'MRPS2', 'MRPS21', 'MRPS22', 'MRPS23', 'MRPS24', 'MRPS25', 'MRPS26',
    'MRPS27', 'MRPS28', 'MRPS30', 'MRPS31', 'MRPS33', 'MRPS34', 'MRPS35', 'MRPS36', 'MRPS5',
    'MRPS6', 'MRPS7', 'MRPS9', 'MRRF', 'MSRA', 'MSRB2', 'MSRB3', 'MT-ATP6', 'MT-CO1', 'MT-CO2',
    'MT-CO3', 'MT-CYB', 'MT-ND1', 'MT-ND2', 'MT-ND4', 'MT-ND5', 'MT-ND6', 'MTCH1', 'MTCH2',
    'MTERF1', 'MTERF3', 'MTERF4', 'MTFMT', 'MTG1', 'MTHFD2', 'MTIF2', 'MTIF3', 'MTO1', 'MTPAP',
    'MTRF1L', 'MTX1', 'MTX2', 'MTX3', 'MUL1', 'MUT', 'MYO19', 'NADK2', 'NARS2', 'NBR1',
    'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2', 'NDUFA3', 'NDUFA4',
    'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFAF1', 'NDUFAF2',
    'NDUFAF3', 'NDUFAF4', 'NDUFAF5', 'NDUFAF7', 'NDUFB1', 'NDUFB10', 'NDUFB11', 'NDUFB3',
    'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9', 'NDUFC2;KCTD14;NDUFC2-KCTD14',
    'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS6', 'NDUFS7', 'NDUFS8', 'NDUFV1',
    'NDUFV2', 'NDUFV3', 'NFS1', 'NFU1', 'NIPSNAP1', 'NIT2', 'NLN', 'NLRX1', 'NME3', 'NME4',
    'NNT', 'NOA1', 'NRD1', 'NSUN2', 'NSUN4', 'NT5DC2', 'NUBPL', 'NUDT19', 'NUDT5', 'NUDT9',
    'OAT', 'OCIAD2', 'OGDH', 'OPA1', 'OSBPL1A', 'OSGEPL1', 'OXA1L', 'OXCT1', 'OXSM', 'PAICS',
    'PAM16;CORO7-PAM16', 'PARK7', 'PARL', 'PARS2', 'PC', 'PCBD2', 'PCCA', 'PCCB', 'PCK2',
    'PDE12', 'PDF', 'PDHA1', 'PDHB', 'PDHX', 'PDK1', 'PDK2', 'PDP1', 'PDPR', 'PDSS2', 'PEO1',
    'PET100', 'PGAM5', 'PGS1', 'PHB', 'PHB2', 'PISD', 'PITRM1', 'PLSCR3;TMEM256-PLSCR3',
    'PMPCA', 'PMPCB', 'PNKD', 'PNPLA8', 'PNPO', 'PNPT1', 'POLDIP2', 'POLG', 'POLG2', 'POLRMT',
    'PPA2', 'PPIF', 'PPOX', 'PRDX2', 'PRDX3', 'PRDX4', 'PRDX5', 'PRDX6', 'PRKACA;KIN27',
    'PROSC', 'PTCD3', 'PTGES2', 'PTPMT1', 'PUS1', 'PUSL1', 'PYCR1', 'PYCR2', 'QDPR',
    'QIL1;C19orf70', 'QRSL1', 'RAB24', 'RARS2', 'RBFA', 'RDH13', 'RDH14', 'REXO2', 'RHOT1',
    'RHOT2', 'RMDN3', 'RMND1', 'RNASEH1', 'RNMTL1', 'RPUSD4', 'RSAD1', 'SAMM50', 'SARS2',
    'SCO1', 'SCO2', 'SCP2', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SELO', 'SFXN1', 'SFXN3',
    'SHMT2', 'SIRT5', 'SLC25A1', 'SLC25A10', 'SLC25A11', 'SLC25A12', 'SLC25A13', 'SLC25A15',
    'SLC25A16', 'SLC25A19', 'SLC25A20', 'SLC25A21', 'SLC25A22', 'SLC25A23', 'SLC25A24',
    'SLC25A25', 'SLC25A26', 'SLC25A29', 'SLC25A3', 'SLC25A30', 'SLC25A32', 'SLC25A4',
    'SLC25A42', 'SLC25A46', 'SLC25A5', 'SLC25A51', 'SLC25A6', 'SLC30A9', 'SLIRP', 'SMDT1',
    'SNAP29', 'SOD1', 'SOD2', 'SPG7', 'SPIRE1', 'SPR', 'SPTLC2', 'SQRDL', 'SSBP1', 'STARD7',
    'STOML2', 'SUCLA2', 'SUCLG1', 'SUCLG2', 'SUGCT', 'SUOX', 'SUPV3L1', 'SURF1', 'SYNJ2BP',
    'TACO1', 'TAMM41', 'TARS2', 'TAZ', 'TBRG4', 'TEFM', 'TFAM', 'TFB1M', 'TFB2M', 'THEM4',
    'TIMM10', 'TIMM13', 'TIMM17A', 'TIMM17B', 'TIMM21', 'TIMM22', 'TIMM23', 'TIMM44', 'TIMM50',
    'TIMM8A', 'TIMM9', 'TIMMDC1', 'TK2', 'TMEM11', 'TMEM126A', 'TMEM70', 'TMLHE', 'TOMM20',
    'TOMM22', 'TOMM34', 'TOMM40', 'TOMM40L', 'TOMM7', 'TOMM70A', 'TRAP1', 'TRMT10C', 'TRMT2B',
    'TRMT5', 'TRMT61B', 'TRMU', 'TRNT1', 'TSFM', 'TST', 'TTC19', 'TUFM', 'TXN2', 'TXNRD1',
    'TXNRD2', 'UNG', 'UQCC1', 'UQCC2', 'UQCR10', 'UQCR11', 'UQCRB', 'UQCRC1', 'UQCRC2',
    'UQCRFS1;UQCRFS1P1', 'UQCRH;UQCRHL', 'UQCRQ', 'USMG5', 'VARS2', 'VDAC1', 'VDAC2', 'VDAC3',
    'WARS2', 'WBSCR16', 'XPNPEP3', 'YARS2', 'YME1L1', 'ZADH2'
}  # n=735

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

COMPLEXI_GENES = {
    'ACAD9', 'AIFM1', 'ATP5SL', 'COA1', 'ECSIT', 'FOXRED1', 'MT-ND1', 'MT-ND2', 'MT-ND4',
    'MT-ND5', 'MT-ND6', 'NDUFA1', 'NDUFA10', 'NDUFA11', 'NDUFA12', 'NDUFA13', 'NDUFA2',
    'NDUFA3', 'NDUFA5', 'NDUFA6', 'NDUFA7', 'NDUFA8', 'NDUFA9', 'NDUFAB1', 'NDUFAF1',
    'NDUFAF2', 'NDUFAF3', 'NDUFAF4', 'NDUFAF5', 'NDUFAF7', 'NDUFB1', 'NDUFB10', 'NDUFB11',
    'NDUFB3', 'NDUFB4', 'NDUFB5', 'NDUFB6', 'NDUFB7', 'NDUFB8', 'NDUFB9',
    'NDUFC2;KCTD14;NDUFC2-KCTD14', 'NDUFS1', 'NDUFS2', 'NDUFS3', 'NDUFS4', 'NDUFS5', 'NDUFS6',
    'NDUFS7', 'NDUFS8', 'NDUFV1', 'NDUFV2', 'NDUFV3', 'NUBPL', 'TIMMDC1', 'TMEM126A', 'TMEM70'
}  # n=56

COMPLEXII_GENES = {
    'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD'
}  # n=5

COMPLEXIII_GENES = {
    'BCS1L', 'CYC1', 'LYRM7', 'MT-CYB', 'TTC19', 'UQCC1', 'UQCC2', 'UQCR10', 'UQCR11', 'UQCRB',
    'UQCRC1', 'UQCRC2', 'UQCRFS1;UQCRFS1P1', 'UQCRH;UQCRHL', 'UQCRQ'
}  # n=15

COMPLEXIV_GENES = {
    'CMC1', 'COA1', 'COA3', 'COA4', 'COA6', 'COA7', 'COX11', 'COX15', 'COX18', 'COX19',
    'COX20', 'COX4I1', 'COX5A', 'COX5B', 'COX6A1', 'COX6B1', 'COX6C', 'COX7A1', 'COX7A2',
    'COX7A2L', 'COX7B', 'COX7C', 'HIGD1A', 'MT-CO1', 'MT-CO2', 'MT-CO3', 'NDUFA4', 'PET100',
    'SCO1', 'SCO2', 'SURF1', 'TACO1', 'TIMM21'
}  # n=33

COMPLEXV_GENES = {
    'ATP5A1', 'ATP5B', 'ATP5C1', 'ATP5D', 'ATP5E;ATP5EP2', 'ATP5F1', 'ATP5G1;ATP5G3;ATP5G2',
    'ATP5H', 'ATP5I', 'ATP5J', 'ATP5J2', 'ATP5L', 'ATP5O', 'ATP5S', 'ATPAF1', 'ATPAF2',
    'ATPIF1', 'MP68;C14orf2', 'MT-ATP6', 'TMEM70', 'USMG5'
}  # n=21

FAO_ISO_GENES = {
    'ACAA1', 'ACAA2', 'ACACA', 'ACAD10', 'ACAD11', 'ACADM', 'ACADS', 'ACADSB', 'ACADVL',
    'ACAT1', 'ACOT13', 'ACOT7', 'ACOT9', 'ACP6', 'ACSF2', 'ACSF3', 'ACSL1', 'ACSS3', 'AGK',
    'AGPAT5', 'CBR4', 'CPT1A', 'CPT2', 'CRAT', 'CROT', 'CYB5R3', 'CYP27A1', 'DBI', 'DECR1',
    'DHRS1', 'ECH1', 'ECHDC1', 'ECHS1', 'ECI1;DCI', 'ECI2', 'EHHADH', 'ETFA', 'ETFB', 'ETFDH',
    'FASN', 'FDPS', 'FDX1', 'FDXR', 'GCSH', 'GPAM', 'HADH', 'HADHA', 'HADHB', 'HINT2',
    'HSD17B10', 'HSD17B4', 'IDI1', 'LACTB', 'LYPLA1', 'LYPLAL1', 'MCAT', 'MCEE', 'MECR',
    'MGST3', 'MLYCD', 'MUT', 'NDUFAB1', 'OSBPL1A', 'OXSM', 'PCCA', 'PCCB', 'PGS1', 'PISD',
    'PLSCR3;TMEM256-PLSCR3', 'PNPLA8', 'PRDX6', 'PTGES2', 'PTPMT1', 'SCP2', 'SLC25A1',
    'SLC25A20', 'SPTLC2', 'STARD7', 'TAMM41', 'TAZ', 'THEM4', 'ZADH2'
}  # n=82

MITO_GE_ISO_GENES = {
    'AARS2', 'ALKBH1', 'APEX1', 'ATAD3A', 'ATAD3B', 'AURKAIP1', 'C6orf203', 'CARS2',
    'CDK5RAP1', 'CHCHD1', 'COA3', 'DAP3', 'DARS2', 'DDX28', 'DHX30', 'DUS2', 'EARS2', 'ELAC2',
    'ENDOG', 'ERAL1', 'EXD2', 'EXOG', 'FARS2', 'FASTKD1', 'FASTKD2', 'FASTKD5', 'GADD45GIP1',
    'GARS', 'GATB;PET112', 'GFM1', 'GFM2', 'GRSF1', 'GTPBP10', 'GTPBP3', 'GUF1', 'HARS2',
    'HEMK1', 'HSD17B10', 'IARS2', 'ICT1', 'KARS', 'KIAA0391', 'LACTB2', 'LARS2', 'LIG3',
    'LRPPRC', 'MALSU1', 'MARS2', 'METTL15', 'METTL17', 'MGME1', 'MIEF1', 'MPV17L2', 'MRM1',
    'MRPL1', 'MRPL10', 'MRPL11', 'MRPL12', 'MRPL13', 'MRPL14', 'MRPL15', 'MRPL16', 'MRPL17',
    'MRPL18', 'MRPL19', 'MRPL2', 'MRPL20', 'MRPL21', 'MRPL22', 'MRPL23', 'MRPL24', 'MRPL27',
    'MRPL28', 'MRPL3', 'MRPL30', 'MRPL32', 'MRPL33', 'MRPL34', 'MRPL35', 'MRPL37', 'MRPL38',
    'MRPL39', 'MRPL4', 'MRPL40', 'MRPL41', 'MRPL42', 'MRPL43', 'MRPL44', 'MRPL46', 'MRPL47',
    'MRPL48', 'MRPL49', 'MRPL50', 'MRPL51', 'MRPL53', 'MRPL54', 'MRPL55', 'MRPL57', 'MRPL9',
    'MRPS10', 'MRPS11', 'MRPS12', 'MRPS14', 'MRPS15', 'MRPS16', 'MRPS17;hCG_1984214',
    'MRPS18A', 'MRPS18B', 'MRPS18C', 'MRPS2', 'MRPS21', 'MRPS22', 'MRPS23', 'MRPS24', 'MRPS25',
    'MRPS26', 'MRPS27', 'MRPS28', 'MRPS30', 'MRPS31', 'MRPS33', 'MRPS34', 'MRPS35', 'MRPS36',
    'MRPS5', 'MRPS6', 'MRPS7', 'MRPS9', 'MRRF', 'MTERF1', 'MTERF3', 'MTERF4', 'MTFMT', 'MTG1',
    'MTIF2', 'MTIF3', 'MTO1', 'MTPAP', 'MTRF1L', 'NARS2', 'NOA1', 'NSUN2', 'NSUN4', 'OSGEPL1',
    'OXA1L', 'PARS2', 'PDE12', 'PDF', 'PEO1', 'PNPT1', 'POLDIP2', 'POLG', 'POLG2', 'POLRMT',
    'PPA2', 'PTCD3', 'PUS1', 'PUSL1', 'QRSL1', 'RARS2', 'RBFA', 'REXO2', 'RMND1', 'RNASEH1',
    'RNMTL1', 'RPUSD4', 'SARS2', 'SLIRP', 'SSBP1', 'SUPV3L1', 'TACO1', 'TARS2', 'TBRG4',
    'TEFM', 'TFAM', 'TFB1M', 'TFB2M', 'TIMM21', 'TRMT10C', 'TRMT2B', 'TRMT5', 'TRMT61B',
    'TRMU', 'TRNT1', 'TSFM', 'TUFM', 'UNG', 'VARS2', 'WARS2', 'WBSCR16', 'YARS2'
}  # n=191

BCAA_ISO_GENES = {
    'ACAD8', 'ACADSB', 'ACAT1', 'ALDH6A1', 'BCAT2', 'BCKDHA', 'BCKDHB', 'BCKDK', 'DBT', 'DLD',
    'ECHS1', 'ETFA', 'ETFB', 'ETFDH', 'HADHA', 'HIBADH', 'HIBCH', 'HMGCL', 'HSD17B10', 'IVD',
    'MCCC1', 'MCCC2'
}  # n=22


# ---- EV table ----
ISO_MITO_DATA_PATH = find_file(
    ['EV_Table_2*.xlsx'],
    'isolated-mito TMT solubility/aggregation data')

print('Required inputs found:')
print(f'  ISO_MITO_DATA_PATH: {ISO_MITO_DATA_PATH}')
print()

# ── Load data ─────────────────────────────────────────────────────────────────
# sheet_name=0 (by position) + column-whitespace strip: prior EV tables in this
# project have shown up with renamed sheets / stray leading spaces in different
# exported copies of the "same" file.
data_df = pd.read_excel(ISO_MITO_DATA_PATH, sheet_name=0)
data_df.columns = data_df.columns.str.strip()

FCm = 'log2 FC TOTAL Sen_Pro'
dfm = data_df.copy()

MITO_ISO = MITOPROTEOME_ISO

detected_iso = set(dfm['Gene names'].dropna())
ref_genes = MITO_ISO & detected_iso
ref_fc_iso = dfm[dfm['Gene names'].isin(ref_genes)][FCm].dropna().values
ref_median = np.median(ref_fc_iso)
dfm = dfm.copy()
dfm['FC_rel'] = dfm[FCm] - ref_median


# ── Gene sets ─────────────────────────────────────────────────────────────────
def relg(genes):
  return dfm[dfm['Gene names'].isin(genes)]['FC_rel'].dropna().values


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
ci_m = COMPLEXI_GENES & detected_iso
cii_m = COMPLEXII_GENES & detected_iso
ciii_m = COMPLEXIII_GENES & detected_iso
civ_m = COMPLEXIV_GENES & detected_iso
cv_m = COMPLEXV_GENES & detected_iso
mito_ge_m = MITO_GE_ISO_GENES & detected_iso
bcaa_m = BCAA_ISO_GENES & detected_iso
fao_m = FAO_ISO_GENES & detected_iso
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
    ('Complex I', relg(ci_m), ci_m, SUB_COLOR['CI']),
    ('Complex II', relg(cii_m), cii_m, SUB_COLOR['CII']),
    ('Complex III', relg(ciii_m), ciii_m, SUB_COLOR['CIII']),
    ('Complex IV', relg(civ_m), civ_m, SUB_COLOR['CIV']),
    ('Complex V', relg(cv_m), cv_m, SUB_COLOR['CV']),
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