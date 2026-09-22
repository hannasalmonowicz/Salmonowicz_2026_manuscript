"""
Fig 5G -- ISR/ATF4 focused metabolic heatmap, 4 contrasts.

Inputs -- EV tables:
  EV_Table_8*.xlsx -- Prolif CAM, Senescence, and Sen CAM contrasts
  EV_Table_5*.xlsx -- Carbon Source Gal/Glu contrast (sheet 'KSz10_ProteinGroups',
                      with the 'Sen Galactose vs Sen High Glucose' Log2FC/p-value
                      columns this script fuzzy-matches on)
"""
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

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


MAIN_PATH = find_file(['EV_Table_8*.xlsx'], 'WC CAM IMT proteomics (Prolif CAM/Senescence/Sen CAM contrasts)')
GAL_PATH = find_file(['EV_Table_5*.xlsx'], 'senescent GAL-exposure proteomics (Carbon Source Gal/Glu contrast)')
print(f'Using MAIN_PATH: {MAIN_PATH}')
print(f'Using GAL_PATH: {GAL_PATH}\n')

# -------------------------------
# 1️⃣ Load & Clean Data
# -------------------------------
# sheet_name=0 (by position): prior EV tables in this project have shown up with
# renamed sheets in different exported copies of the "same" file.
df_main = pd.read_excel(MAIN_PATH, sheet_name=0)
df_gal = pd.read_excel(GAL_PATH, sheet_name=0)

def clean_gene_names(df):
    df = df.copy()
    df.columns = df.columns.str.strip()
    gene_col = next((col for col in df.columns if "gene names" in col.lower()), "Gene names")
    df["Gene names"] = df[gene_col].astype(str).str.split(";").str[0].str.strip()
    return df.groupby("Gene names", as_index=False).mean(numeric_only=True).set_index("Gene names")

df_main_clean = clean_gene_names(df_main)
df_gal_clean = clean_gene_names(df_gal)

# -------------------------------
# 2️⃣ Focused ISR Metabolic Gene Set — trimmed to the 10 genes actually
# shown in the published Fig 5G panel (SHMT2, MTHFD2, MTHFD1L, ALDH1L2,
# PHGDH, PSPH, PYCR1, PYCR2, GLS, SLC7A11), in that exact row order.
# "Redox & Nitrogen" previously had 6 genes; only GLS is in the published
# panel. SLC7A11 is new — added here under "Redox & Nitrogen" as the
# cystine/glutamate antiporter feeding glutathione synthesis. Flag if the
# real panel grouped/colored it differently.
# -------------------------------
gene_groups = {
    "1-Carbon Metabolism": ["SHMT2", "MTHFD2", "MTHFD1L", "ALDH1L2"],
    "Serine Biosynthesis": ["PHGDH", "PSPH"],
    "Proline Biosynthesis": ["PYCR1", "PYCR2"],
    "Redox & Nitrogen":    ["GLS", "SLC7A11"],
}
target_genes = [gene for sublist in gene_groups.values() for gene in sublist]
gene_to_group = {gene: group for group, genes in gene_groups.items() for gene in genes}

# -------------------------------
# 3️⃣ Build Matrices (ORDER: Prolif CAM, Senescence, Sen CAM, Gal)
# -------------------------------
valid_idx = df_main_clean.index.union(df_gal_clean.index).intersection(target_genes)
valid_idx = [g for g in target_genes if g in valid_idx]

def get_col_data(df, search_term, genes):
    fc_col = next((c for c in df.columns if search_term.lower() in c.lower() and "log2fc" in c.lower()), None)
    p_col = next((c for c in df.columns if search_term.lower() in c.lower() and "p-value" in c.lower()), None)
    if fc_col:
        data = df[fc_col].reindex(genes).to_frame().clip(-1, 1)
        sig = ["*" if val >= 1.3 else "" for val in df[p_col].reindex(genes).fillna(0)]
        return data, pd.DataFrame(sig, index=genes, columns=[fc_col])
    return pd.DataFrame(0.0, index=genes, columns=[search_term]), pd.DataFrame("", index=genes, columns=[search_term])

# Dataset Mapping
data_1, sig_1 = get_col_data(df_main_clean, "Prolif CAM vs Prolif CTRL", valid_idx)
data_2, sig_2 = get_col_data(df_main_clean, "Sen CTRL vs Prolif CTRL", valid_idx)
data_3, sig_3 = get_col_data(df_main_clean, "Sen CAM vs Sen CTRL", valid_idx)
data_4, sig_4 = get_col_data(df_gal_clean, "Sen Galactose vs Sen High Glucose", valid_idx)

data_1.columns = ["Prolif CAM\n(CAM/Prolif)"]
data_2.columns = ["Senescence\n(Sen/Prolif)"]
data_3.columns = ["Sen CAM\n(CAM/Sen)"]
data_4.columns = ["Carbon Source\n(Gal/Glu)"]

# -------------------------------
# 4️⃣ Plotting
# -------------------------------
fig = plt.figure(figsize=(12, 9))
gs = fig.add_gridspec(1, 6, width_ratios=[0.2, 1, 1, 1, 1, 2], wspace=0.6)

ax_cbar = fig.add_subplot(gs[0])
axes = [fig.add_subplot(gs[i+1]) for i in range(4)]
datasets, sigs = [data_1, data_2, data_3, data_4], [sig_1, sig_2, sig_3, sig_4]

def draw_borders(ax, shape):
    for i in range(shape[0] + 1): ax.axhline(i, color='black', lw=1.5)
    for j in range(shape[1] + 1): ax.axvline(j, color='black', lw=1.5)

heat_params = dict(cmap="coolwarm", center=0, vmin=-1.0, vmax=1.0, square=True, fmt="")

for i, ax in enumerate(axes):
    sns.heatmap(datasets[i], ax=ax, annot=sigs[i], cbar=(i==0), cbar_ax=(ax_cbar if i==0 else None), **heat_params)
    draw_borders(ax, datasets[i].shape)
    plt.setp(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticks([]); ax.set_xlabel(""); ax.set_ylabel("")

# Final axis formatting
axes[3].set_yticks(np.arange(len(valid_idx)) + 0.5)
axes[3].set_yticklabels(valid_idx, rotation=0, fontsize=10)
axes[3].yaxis.tick_right()
axes[3].tick_params(axis='y', pad=15)

group_colors = {
    "1-Carbon Metabolism": "#333333", "Serine Biosynthesis": "#E69F00",
    "Proline Biosynthesis": "#56B4E9", "Redox & Nitrogen": "#009E73"
}

for i, gene in enumerate(valid_idx):
    color = group_colors[gene_to_group[gene]]
    axes[0].add_patch(plt.Rectangle((-0.6, i), 0.2, 1, color=color, transform=axes[0].get_yaxis_transform(), clip_on=False))

ax_cbar.set_title("Log2FC", fontsize=9, pad=10)
fig.suptitle(f"Mitochondrial ISR Focused Panel (n={len(valid_idx)})", fontsize=15, y=0.98)

handles = [plt.Rectangle((0,0),1,1, color=color) for color in group_colors.values()]
axes[0].legend(handles, group_colors.keys(), title="Pathway", bbox_to_anchor=(-1.8, 1.02), loc='upper right', frameon=False)

plt.savefig("ISR_Focused_Panel.svg", bbox_inches="tight")
# JPG has no transparency channel, so force a white background rather than
# letting matplotlib's default transparent one turn black on export.
plt.savefig("ISR_Focused_Panel.jpg", dpi=300, bbox_inches="tight", facecolor="white")
plt.show()
