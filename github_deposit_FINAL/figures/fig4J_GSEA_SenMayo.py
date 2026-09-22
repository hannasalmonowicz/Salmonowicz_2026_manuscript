"""
Fig 4J -- Compact GSEA enrichment plot — pure matplotlib, gseapy 1.1.9
Fixes: hits is an integer index array — use directly, never np.where()

Input -- one EV table only:
      EV_Table_5*.xlsx -- 'Gene names' plus the galactose-contrast Log2FC/p-value columns.
      The fuzzy Log2FC/-Log10 p-value 'galactose' column matching below works whether or
      not there's a stray double space in the header (the real file has one: '-Log10  p-value
      Sen Galactose vs Sen High Glucose').

The SenMayo gene set (published, external reference data, not generated for this
manuscript) is hardcoded below as SENMAYO_RAW_GENE_STRING -- no separate file needed.
Sourced directly from Saul et al. 2022 (Nat Commun), Supplementary Data 4, sheet 'human'
-- 125 genes -- rather than MSigDB's mirror of it, which is missing one gene (CCL3L1).
Re-run against the real EV table with this 125-gene list: NES=1.441, FDR q=0.0316 --
unchanged from the 124-gene version and still exactly matching the manuscript's own
reported NES=1.44 (CCL3L1 isn't among the 29 SenMayo genes detected in this dataset, so
the extra gene doesn't move the result -- it's included for correct provenance, not to
change the number).
"""

import re
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import seaborn as sns

try:
    import gseapy as gp
except ImportError:
    raise ImportError("pip install gseapy")

# ── CONFIG ────────────────────────────────────────────────────────────────────
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


FILE_GAL = find_file(
    ['EV_Table_5*.xlsx'],
    'senescent GAL-exposure proteomics (Sen Galactose vs Sen High Glucose)')
TARGET_SET     = "SenMayo"
PERM           = 1000
SEED           = 42

FIG_W, FIG_H   = 3.2, 2.8

COL_POS        = "#C0392B"
COL_NEG        = "#2471A3"
COL_HIT        = "#2b2d42"
COL_ZERO       = "#aaaaaa"

print('Required inputs found:')
print(f'  FILE_GAL: {FILE_GAL}')
print()

# ── SenMayo gene set (Saul et al. 2022, Nat Commun, Supplementary Data 4, sheet
# 'human'; PMID 35974106) -- published, external reference data, hardcoded here
# instead of a separate bundled file. See verification/PROVENANCE.md. ──────────
SENMAYO_RAW_GENE_STRING = """
ACVR1B,ANG,ANGPT1,ANGPTL4,AREG,AXL,BEX3,BMP2,BMP6,C3,CCL1,CCL13,CCL16,CCL2,CCL20,CCL24,
CCL26,CCL3,CCL3L1,CCL4,CCL5,CCL7,CCL8,CD55,CD9,CSF1,CSF2,CSF2RB,CST4,CTNNB1,CTSB,CXCL1,
CXCL10,CXCL12,CXCL16,CXCL2,CXCL3,CXCL8,CXCR2,DKK1,EDN1,EGF,EGFR,EREG,ESM1,ETS2,FAS,FGF1,
FGF2,FGF7,GDF15,GEM,GMFG,HGF,HMGB1,ICAM1,ICAM3,IGF1,IGFBP1,IGFBP2,IGFBP3,IGFBP4,IGFBP5,
IGFBP6,IGFBP7,IL10,IL13,IL15,IL18,IL1A,IL1B,IL2,IL32,IL6,IL6ST,IL7,INHA,IQGAP2,ITGA2,
ITPKA,JUN,KITLG,LCP1,MIF,MMP1,MMP10,MMP12,MMP13,MMP14,MMP2,MMP3,MMP9,NAP1L4,NRG1,PAPPA,
PECAM1,PGF,PIGF,PLAT,PLAU,PLAUR,PTBP1,PTGER2,PTGES,RPS6KA5,SCAMP4,SELPLG,SEMA3F,SERPINB4,
SERPINE1,SERPINE2,SPP1,SPX,TIMP2,TNF,TNFRSF10C,TNFRSF11B,TNFRSF1A,TNFRSF1B,TUBGCP2,VEGFA,
VEGFC,VGF,WNT16,WNT2
"""
saul_genes = sorted(set(
    g.strip().upper() for g in re.split(r'[\s,]+', SENMAYO_RAW_GENE_STRING) if g.strip()
))

# ── LOAD & RANK ───────────────────────────────────────────────────────────────
# sheet_name=0 (by position) + column-whitespace strip: EV tables in this project
# can have renamed sheets or stray leading spaces in different exported copies of
# the "same" file. The fuzzy lfc_col/p_col matching below already tolerates the
# double-space quirk in the real column header.
df = pd.read_excel(FILE_GAL, sheet_name=0)
df.columns = df.columns.str.strip()

lfc_col = next((c for c in df.columns
                if "log2fc" in c.lower() and "galactose" in c.lower()), None)
p_col   = next((c for c in df.columns
                if "-log10" in c.lower() and "p-value" in c.lower()
                and "galactose" in c.lower()), None)
if lfc_col is None or p_col is None:
    raise KeyError(
        f"Could not find the Log2FC / -Log10 p-value 'galactose' columns in "
        f"{FILE_GAL}. Columns present: {list(df.columns)}"
    )

df['rank_metric'] = np.sign(df[lfc_col]) * df[p_col]
df['Gene_Symbol']  = (df['Gene names']
                      .str.split(';').str[0].str.strip().str.upper())
rank_df = (df[['Gene_Symbol', 'rank_metric']]
           .dropna()
           .groupby('Gene_Symbol').mean()
           .sort_values('rank_metric', ascending=False))

# ── RUN GSEA ──────────────────────────────────────────────────────────────────
pre_res = gp.prerank(
    rnk=rank_df,
    gene_sets={TARGET_SET: saul_genes},
    threads=4,
    permutation_num=PERM,
    outdir=None,
    seed=SEED,
    min_size=1,
    max_size=5000,
)

res_key  = (TARGET_SET if TARGET_SET in pre_res.results
            else list(pre_res.results.keys())[0])
r        = pre_res.results[res_key]

# ── EXTRACT DATA ──────────────────────────────────────────────────────────────
RES      = np.array(r['RES'])       # length = n_ranked = 5031
hit_pos  = np.array(r['hits'])      # integer rank positions — use DIRECTLY
nes      = float(r['nes'])
fdr      = float(r['fdr'])
pval     = float(r['pval'])

n        = len(RES)                 # 5031
x_axis   = np.arange(n)

print(f"NES     : {nes:.3f}")
print(f"p-value : {pval:.4f}")
print(f"FDR q   : {fdr:.4f}")
print(f"Hits    : {len(hit_pos)} genes at positions {hit_pos[:5]}...{hit_pos[-5:]}")
print(f"x-axis  : 0 to {n-1}")

# ── FIGURE ────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family'   : 'Arial',
    'font.size'     : 7,
    'axes.linewidth': 0.7,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'xtick.major.size' : 2.5,
    'ytick.major.size' : 2.5,
    'pdf.fonttype'  : 42,
    'svg.fonttype'  : 'none',
})

fig = plt.figure(figsize=(FIG_W, FIG_H))

gs = gridspec.GridSpec(
    2, 1,
    height_ratios=[6, 1],
    hspace=0.05,
    left=0.20, right=0.96,
    top=0.78, bottom=0.16,
)

ax_es  = fig.add_subplot(gs[0])
ax_rug = fig.add_subplot(gs[1], sharex=ax_es)

# ─── enrichment score curve ───────────────────────────────────────────────────
col = COL_POS if nes >= 0 else COL_NEG

ax_es.plot(x_axis, RES, color=col, lw=1.2, zorder=3, solid_capstyle='round')
ax_es.fill_between(x_axis, RES, 0, color=col, alpha=0.12, zorder=1)
ax_es.axhline(0, color=COL_ZERO, lw=0.6, linestyle='--', zorder=2)

peak_x = int(np.argmax(np.abs(RES)))
ax_es.axvline(peak_x, color=col, lw=0.6, linestyle=':', alpha=0.6, zorder=2)

# Asymmetric y limits — based on actual data range, not forced symmetry
# This removes the white space below zero when NES is positive
y_pad  = (RES.max() - RES.min()) * 0.15          # 15% padding
y_min  = RES.min() - y_pad                        # actual minimum of curve
y_max  = RES.max() + y_pad                        # actual maximum of curve
# Always include zero so the reference line is visible
y_min  = min(y_min, -y_pad)
y_max  = max(y_max,  y_pad)
ax_es.set_ylim(y_min, y_max)
ax_es.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
ax_es.set_ylabel("Enrichment score", fontsize=6, labelpad=3)
ax_es.tick_params(axis='y', labelsize=5.5, pad=2)
ax_es.tick_params(axis='x', bottom=False, labelbottom=False)
sns.despine(ax=ax_es, top=True, right=True, bottom=True)

# Stats box
fdr_str  = f"{fdr:.3f}" if fdr >= 0.001 else f"{fdr:.2e}"
pval_str = f"{pval:.3f}" if pval >= 0.001 else f"{pval:.2e}"
ax_es.text(
    0.97, 0.97,
    f"NES = {nes:.2f}\nFDR = {fdr_str}\np = {pval_str}",
    transform=ax_es.transAxes,
    fontsize=5.5, va='top', ha='right', linespacing=1.6,
    bbox=dict(facecolor='white', edgecolor='none', alpha=0.0, pad=1.5),
)

ax_es.text(
    0.03, 0.97, TARGET_SET,
    transform=ax_es.transAxes,
    fontsize=6.5, fontweight='bold', va='top', ha='left', color=col,
)

# ─── hit rug — using hit_pos DIRECTLY as x coordinates ───────────────────────
# Each value in hit_pos is already a rank index in 0..n-1
# Draw as thin vertical lines using vlines for efficiency
ax_rug.vlines(
    hit_pos,            # x positions — the direct integer rank indices
    ymin=0, ymax=1,
    color=COL_HIT,
    lw=0.8,             # slightly thicker than before so lines are visible
    alpha=0.75,
)

ax_rug.set_xlim(0, n - 1)
ax_rug.set_ylim(0, 1)
ax_rug.set_yticks([])
ax_rug.set_xlabel("Rank", fontsize=6, labelpad=3)
ax_rug.tick_params(axis='x', labelsize=5.5, pad=2)
sns.despine(ax=ax_rug, top=True, right=True, left=True)

ax_rug.text(0.01, -0.55, "← high in GAL",
            transform=ax_rug.transAxes,
            fontsize=5, ha='left', va='top',
            color='#666666', style='italic')
ax_rug.text(0.99, -0.55, "low in GAL →",
            transform=ax_rug.transAxes,
            fontsize=5, ha='right', va='top',
            color='#666666', style='italic')

# ─── title ────────────────────────────────────────────────────────────────────
fig.text(
    0.58, 0.89, "Sen GAL vs CTRL",
    ha='center', va='bottom',
    fontsize=7, fontweight='bold',
)

# ── SAVE ─────────────────────────────────────────────────────────────────────
for fmt in ['svg', 'png']:
    out = f"GSEA_compact_{TARGET_SET}_final.{fmt}"
    plt.savefig(out, format=fmt, dpi=300, bbox_inches='tight')
    print(f"✅ Saved {out}")

plt.show()