"""
Fig 5E -- Compact GSEA enrichment plot — pure matplotlib, gseapy 1.1.9
Comparison    : Sen CAM vs Sen CTRL (EXACT MATCH)
Gene Set      : Custom ISR/ATF4 Target Gene List

Input -- EV table:
  EV_Table_8*.xlsx -- has the exact 'Log2FC Sen CAM vs Sen CTRL' / '-Log10 p-value Sen
  CAM vs Sen CTRL' column names this script requires verbatim (EXACT_LFC_COL/
  EXACT_PVAL_COL below), so no fuzzy matching is needed.
"""

import glob
import re
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


FILE_PATH = find_file(
    ['EV_Table_8*.xlsx'],
    'WC CAM IMT proteomics (Sen CAM vs Sen CTRL, CAM-exposed dataset)')
print(f'Using FILE_PATH: {FILE_PATH}\n')
TARGET_SET   = "ISR_ATF4_Signature"
PERM         = 1000
SEED         = 42

FIG_W, FIG_H = 3.2, 2.8

COL_POS      = "#C0392B"
COL_NEG      = "#2471A3"
COL_HIT      = "#2b2d42"
COL_ZERO     = "#aaaaaa"

# ── EXACT TARGET COLUMNS ──────────────────────────────────────────────────────
EXACT_LFC_COL = "Log2FC Sen CAM vs Sen CTRL"
EXACT_PVAL_COL = "-Log10 p-value Sen CAM vs Sen CTRL"

# ── CUSTOM GENE LIST ──────────────────────────────────────────────────────────
RAW_GENE_STRING = """
Mthfd2, Slc7a5, Aldh18a1, Abcc4, Cth, Atf5, Eprs1, Lonp1, Atf4, Atf6, Ddr2, Slc1a4, Iars1, Tsc22d3, Asns, Nfil3, Slc25a33, Pakap, Hspa5, Gtpbp2, Shmt2, Zc3h11a, Snai2, Nfe2l1, Ncoa7, Ypel5, Gfpt1, Nfu1, Cebpb
Trib3, Phgdh, Arhgef2, Psat1, Ddit3, Pycr1, Got1, Rhbdd1, Cars1, Slc20a1, Slc6a9, Lars1, Yars1, Paqr3, Atf3, Ddit4, Ubald2, Tars1, Gars1, Gpt2, Cebpg, Xbp1, Inpp5b, Stc2, Sesn2, Tmem11, Ube2g2, Mxd1, Clcn3, Eif4ebp1
"""

# Clean, normalize to uppercase, and deduplicate genes
custom_genes = sorted(list(set(
    g.strip().upper() for g in re.split(r'[\s,]+', RAW_GENE_STRING) if g.strip()
)))

print(f"📋 Loaded {len(custom_genes)} unique genes into '{TARGET_SET}' gene set.")

# ── LOAD & RANK ───────────────────────────────────────────────────────────────
# sheet_name=0 (by position) + column-whitespace strip: prior EV tables in this
# project have shown up with renamed sheets / stray leading spaces in different
# exported copies of the "same" file.
df = pd.read_excel(FILE_PATH, sheet_name=0)
df.columns = df.columns.str.strip()

# Verify exact column match
if EXACT_LFC_COL not in df.columns or EXACT_PVAL_COL not in df.columns:
    raise KeyError(f"Could not find exact columns '{EXACT_LFC_COL}' or '{EXACT_PVAL_COL}' in Excel file.")

lfc_col = EXACT_LFC_COL
p_col   = EXACT_PVAL_COL

print(f"✅ EXACT MATCH Log2FC Column : {lfc_col}")
print(f"✅ EXACT MATCH P-value Column: {p_col}")

# Ensure numeric types
lfc_numeric = pd.to_numeric(df[lfc_col], errors='coerce')
p_numeric   = pd.to_numeric(df[p_col], errors='coerce')

# Calculate rank metric (-log10(p) * sign(log2FC))
df['rank_metric'] = np.sign(lfc_numeric) * p_numeric

# Target 'Gene names' column from the EV table
df['Gene_Symbol']  = (df['Gene names']
                      .astype(str)
                      .str.split(';').str[0].str.strip().str.upper())

# Rank dataframe
rank_df = (df[['Gene_Symbol', 'rank_metric']]
           .dropna()
           .groupby('Gene_Symbol').mean()
           .sort_values('rank_metric', ascending=False))

# ── RUN GSEA ──────────────────────────────────────────────────────────────────
pre_res = gp.prerank(
    rnk=rank_df,
    gene_sets={TARGET_SET: custom_genes},
    threads=4,
    permutation_num=PERM,
    outdir=None,
    seed=SEED,
    min_size=1,
    max_size=5000,
)

res_key = (TARGET_SET if TARGET_SET in pre_res.results
           else list(pre_res.results.keys())[0])
r       = pre_res.results[res_key]

# ── EXTRACT DATA ──────────────────────────────────────────────────────────────
RES      = np.array(r['RES'])
hit_pos  = np.array(r['hits'])
nes      = float(r['nes'])
fdr      = float(r['fdr'])
pval     = float(r['pval'])

n        = len(RES)
x_axis   = np.arange(n)

# Extract leading edge genes (genes contributing to the ES peak)
lead_genes = set(r.get('matched_genes', []))

# ── INSPECT DETECTED VS MISSING GENES ─────────────────────────────────────────
ranked_gene_list = list(rank_df.index)
hit_genes = [ranked_gene_list[pos] for pos in hit_pos]
missing_genes = sorted(list(set(custom_genes) - set(hit_genes)))

# Create detailed dataframe for detected genes
detected_details = []
for pos in hit_pos:
    gene = ranked_gene_list[pos]
    
    # Retrieve original values from df
    sub = df[df['Gene_Symbol'] == gene]
    log2fc_val = sub[lfc_col].mean() if not sub.empty else np.nan
    pval_val   = sub[p_col].mean() if not sub.empty else np.nan
    
    detected_details.append({
        'Rank_Position': pos,
        'Gene_Symbol': gene,
        'Rank_Metric': rank_df.loc[gene, 'rank_metric'],
        'Log2FC': log2fc_val,
        '-Log10_PValue': pval_val,
        'Is_Leading_Edge': gene in lead_genes
    })

hits_df = pd.DataFrame(detected_details).sort_values('Rank_Position')

print("\n" + "="*70)
print(f"📊 GENE DETECTABILITY SUMMARY FOR '{TARGET_SET}' (Sen CAM vs Sen CTRL)")
print("="*70)
print(f"Total input genes       : {len(custom_genes)}")
print(f"Detected in proteomics  : {len(hits_df)} ({len(hits_df)/len(custom_genes)*100:.1f}%)")
print(f"Missing from proteomics : {len(missing_genes)} ({len(missing_genes)/len(custom_genes)*100:.1f}%)")
print(f"Leading Edge genes      : {len(lead_genes)}")

print("\n✅ DETECTED GENES / PROTEINS (Sorted by Rank):")
print(hits_df.to_string(index=False))

if missing_genes:
    print("\n❌ MISSING GENES / PROTEINS (Not detected in dataset):")
    print(", ".join(missing_genes))

# Save detailed report to CSV
report_filename = f"{TARGET_SET}_detected_genes_EXACT_Sen_CAM_vs_CTRL.csv"
hits_df.to_csv(report_filename, index=False)
print(f"\n💾 Saved detailed gene table to: {report_filename}")
print("="*70 + "\n")

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

es_lim = max(abs(RES)) * 1.25
ax_es.set_ylim(-es_lim, es_lim)
ax_es.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, symmetric=True))
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
ax_rug.vlines(
    hit_pos,
    ymin=0, ymax=1,
    color=COL_HIT,
    lw=0.8,
    alpha=0.75,
)

ax_rug.set_xlim(0, n - 1)
ax_rug.set_ylim(0, 1)
ax_rug.set_yticks([])
ax_rug.set_xlabel("Rank", fontsize=6, labelpad=3)
ax_rug.tick_params(axis='x', labelsize=5.5, pad=2)
sns.despine(ax=ax_rug, top=True, right=True, left=True)

# Annotations for Sen CAM vs Sen CTRL
ax_rug.text(0.01, -0.55, "← high in Sen CAM",
            transform=ax_rug.transAxes,
            fontsize=5, ha='left', va='top',
            color='#666666', style='italic')
ax_rug.text(0.99, -0.55, "high in Sen CTRL →",
            transform=ax_rug.transAxes,
            fontsize=5, ha='right', va='top',
            color='#666666', style='italic')

# ─── title ────────────────────────────────────────────────────────────────────
fig.text(
    0.58, 0.89, "Sen CAM vs Sen CTRL",
    ha='center', va='bottom',
    fontsize=7, fontweight='bold',
)

# ── SAVE ─────────────────────────────────────────────────────────────────────
for fmt in ['svg', 'png']:
    out = f"GSEA_compact_{TARGET_SET}_EXACT_Sen_CAM_vs_CTRL.{fmt}"
    plt.savefig(out, format=fmt, dpi=300, bbox_inches='tight')
    print(f"✅ Saved {out}")

plt.show()