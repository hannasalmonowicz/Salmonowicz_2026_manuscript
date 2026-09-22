#!/usr/bin/env python3
# =====================================================================
# GO Biological Process enrichment bubble plots -- EV Fig. 7B
# Three scenarios: Prolif CAM vs Prolif CTRL, Sen CAM vs Sen CTRL, and
# Sen CAM vs Sen CTRL excluding Prolif CAM effects.
#
# Quantitative data: EV_Table_8 (read by sheet position, not name) --
# 'Gene names', 'Log2FC Prolif CAM vs Prolif CTRL', 'q-value Prolif CAM vs
# Prolif CTRL', 'Log2FC Sen CAM vs Sen CTRL', 'q-value Sen CAM vs Sen CTRL',
# and 'GO_Names_P'.
#
# GO annotation: EV_Table_8's own 'GO_Names_P' column, not a separate file.
# It holds, for each row, the semicolon-joined set of GO Biological Process
# term names for that row's gene(s) -- unioned across every constituent
# symbol when 'Gene names' is a semicolon-joined protein group (e.g.
# 'UQCRFS1;UQCRFS1P1'). It was generated once from the original, unmodified
# QuickGO Biological Process annotation export (GO release 2025-06-01,
# exported 2025-06-29), keeping only human (TAXON ID 9606) Biological
# Process (GO ASPECT 'P') rows -- see verification/PROVENANCE.md for the
# exact derivation and match-rate numbers.
# =====================================================================
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import Normalize
import matplotlib.cm as cm

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


# ---- quantitative data + GO annotation, both from EV_Table_8 ----
EV_TABLE_8_PATH = find_file(['EV_Table_8*.xlsx'], 'EV Table 8')
df = pd.read_excel(EV_TABLE_8_PATH, sheet_name=0)

# GO term column definitions and FDR cut-off
gobp_col = "GO_Names_P"
fdr_cut = 0.05

# Define CAM-only scenarios and their corresponding columns
scenarios = [
    {
        "title": "Prolif CAM vs Prolif CTRL",
        "q_value": "q-value Prolif CAM vs Prolif CTRL",
        "log2fc": "Log2FC Prolif CAM vs Prolif CTRL",
        "flag": "prolif_cam_vs_ctrl"
    },
    {
        "title": "Sen CAM vs Sen CTRL",
        "q_value": "q-value Sen CAM vs Sen CTRL",
        "log2fc": "Log2FC Sen CAM vs Sen CTRL",
        "flag": "sen_cam_vs_ctrl"
    },
    {
        "title": "Sen CAM vs Sen CTRL (Excluding Prolif CAM effects)",
        "q_value": "q-value Sen CAM vs Sen CTRL",
        "log2fc": "Log2FC Sen CAM vs Sen CTRL",
        "flag": "sen_cam_vs_ctrl_excl_prolif_cam"
    }
]

def get_gobp_terms_explode(df_in):
    """Splits GO terms on semicolons and explodes them into individual rows."""
    return df_in[gobp_col].dropna().str.split(';').explode().str.strip()

def get_gobp_enrichment_with_fdr(df_filtered, q_value_col, topn=20):
    """Return top N GO terms with counts and median/min FDR for a given q-value column."""
    if len(df_filtered) == 0:
        return pd.DataFrame(columns=['go_term', 'gene_count', 'median_fdr', 'min_fdr'])

    df_temp = df_filtered[[gobp_col, q_value_col]].dropna(subset=[gobp_col]).copy()
    df_temp[gobp_col] = df_temp[gobp_col].str.split(';')
    df_exploded = df_temp.explode(gobp_col)
    df_exploded[gobp_col] = df_exploded[gobp_col].str.strip()
    df_exploded = df_exploded[(df_exploded[gobp_col].notnull()) & (df_exploded[gobp_col] != '')]

    go_stats = df_exploded.groupby(gobp_col).agg({
        q_value_col: ['count', 'median', 'min']
    }).reset_index()
    go_stats.columns = ['go_term', 'gene_count', 'median_fdr', 'min_fdr']

    go_stats = go_stats.sort_values('gene_count', ascending=False).head(topn)
    return go_stats

def filter_data(df, direction, scenario):
    """Filters the DataFrame based on significance and fold-change direction for CAM scenarios."""

    # Base filter for the current scenario's q_value and log2fc
    q_val_col = scenario["q_value"]
    log2fc_col = scenario["log2fc"]

    sig_mask = (df[q_val_col] < fdr_cut)
    fc_mask = (df[log2fc_col] > 0) if direction == "up" else (df[log2fc_col] < 0)
    filtered_df = df[sig_mask & fc_mask]

    # Apply additional filter for excluding Prolif CAM effects from Sen CAM
    if scenario["flag"] == "sen_cam_vs_ctrl_excl_prolif_cam":
        q_val_prolif_col = "q-value Prolif CAM vs Prolif CTRL"
        log2fc_prolif_col = "Log2FC Prolif CAM vs Prolif CTRL"

        sig_prolif_mask = (df.loc[filtered_df.index, q_val_prolif_col] < fdr_cut)
        fc_prolif_mask = (df.loc[filtered_df.index, log2fc_prolif_col] > 0) if direction == "up" else (df.loc[filtered_df.index, log2fc_prolif_col] < 0)

        filtered_df = filtered_df[~(sig_prolif_mask & fc_prolif_mask)]

    return filtered_df

def plot_bubble_chart(df_down, df_up, scenario_title, q_value_col, filename_flag):
    """
    Generates and saves a bubble plot from the filtered data.
    """
    bubble_down = get_gobp_enrichment_with_fdr(df_down, q_value_col, topn=20)
    bubble_up = get_gobp_enrichment_with_fdr(df_up, q_value_col, topn=20)

    bubble_down['direction'] = "DOWN"
    bubble_up['direction'] = "UP"

    bubble_all = pd.concat([bubble_down, bubble_up], ignore_index=True)

    if bubble_all.empty:
        print(f"No data to plot for {scenario_title}")
        return

    # Order and base positions, increase spacing for better visibility and remove jitter
    bubble_all = bubble_all.sort_values(['direction', 'gene_count'], ascending=[True, False]).reset_index(drop=True)
    bubble_all['y_pos'] = bubble_all.index * 1.5

    # Column layout
    center_x = 0.5
    dodge = 0.13
    x_map = {'DOWN': center_x - dodge, 'UP': center_x + dodge}
    bubble_all['x_pos'] = bubble_all['direction'].map(x_map)

    # -log10 median FDR
    bubble_all['neg_log10_fdr'] = -np.log10(bubble_all['median_fdr'].clip(lower=1e-300))

    # Cap color scale at 90th percentile to compress range
    clean_vals = bubble_all['neg_log10_fdr'].replace([np.inf, -np.inf], np.nan).dropna()
    vmax_cap = np.percentile(clean_vals, 90) if len(clean_vals) > 0 else 1.0
    vmin = clean_vals.min() if len(clean_vals) > 0 else 0.0
    vmax = max(vmax_cap, vmin + 1e-6)

    # Plot
    sns.set_style("whitegrid")
    fig_height = max(6, len(bubble_all) * 0.33)
    fig, ax = plt.subplots(1, 1, figsize=(8.5, fig_height))

    # Define color maps for each direction
    cmap_down = cm.Blues
    cmap_up = cm.Reds
    size_scale = 28

    # Plot per direction with separate color maps
    for direction in ['DOWN', 'UP']:
        subset = bubble_all[bubble_all['direction'] == direction]

        # Use the correct color map based on direction
        cmap = cmap_down if direction == 'DOWN' else cmap_up

        ax.scatter(
            subset['x_pos'],
            subset['y_pos'],
            s=subset['gene_count'] * size_scale,
            c=subset['neg_log10_fdr'],
            cmap=cmap,
            norm=Normalize(vmin=vmin, vmax=vmax),
            edgecolors='black',
            linewidth=0.4,
            alpha=0.85,
            zorder=3
        )

    # Y labels
    ax.set_yticks(bubble_all['y_pos'])
    ax.set_yticklabels(bubble_all['go_term'], fontsize=10)
    ax.invert_yaxis()

    # X axis formatting
    ax.set_xticks([center_x - dodge, center_x + dodge])
    ax.set_xticklabels(['Down', 'Up'], fontsize=13, fontweight='bold')
    ax.set_xlim(center_x - 0.5, center_x + 0.5)
    ax.set_xlabel('', fontsize=12)

    # Spines & grid
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.grid(axis='x', color='lightgray', linestyle='-', linewidth=0.5, alpha=0.4)
    ax.set_axisbelow(True)

    # Create two color bars on the right
    cbar1 = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap_down),
                         ax=ax, shrink=0.55, pad=0.08, aspect=20, location='right')
    cbar1.set_label('-log10(FDR) Down', fontsize=11, fontweight='bold')
    cbar1.ax.tick_params(labelsize=9)


    cbar2 = fig.colorbar(cm.ScalarMappable(norm=Normalize(vmin=vmin, vmax=vmax), cmap=cmap_up),
                         ax=ax, shrink=0.55, pad=0.02, aspect=20, location='right')
    cbar2.set_label('-log10(FDR) Up', fontsize=11, fontweight='bold')
    cbar2.ax.tick_params(labelsize=9)


    # Size legend
    unique_counts = np.array(sorted(bubble_all['gene_count'].astype(float).unique()))
    if unique_counts.size == 0:
        leg_counts = []
    elif unique_counts.size == 1:
        leg_counts = [int(unique_counts[0])]
    elif unique_counts.size == 2:
        leg_counts = [int(unique_counts[0]), int(unique_counts[1])]
    else:
        median_val = int(np.median(unique_counts))
        leg_counts = [int(unique_counts[0]), median_val, int(unique_counts[-1])]

    legend_elems = []
    for c in leg_counts:
        legend_elems.append(plt.scatter([], [], s=c * size_scale, c='lightgray', edgecolors='black', linewidth=0.4))
    legend_labels = [str(int(c)) for c in leg_counts]
    if legend_elems:
        legend = ax.legend(legend_elems, legend_labels, title='Count', loc='lower right', bbox_to_anchor=(0.98, 0.02), frameon=True)
        legend.get_title().set_fontweight('bold')

    plt.tight_layout()

    # Set the main title
    plt.suptitle(f"GO Biological Process Enrichment: {scenario_title}\n(CAM Data)", y=1.01, fontsize=18)

    # Save as SVG (vector) and JPG
    fname_base = f"Bubble plot_{filename_flag}_CAM_data"
    plt.savefig(f"{fname_base}.svg", format='svg', bbox_inches='tight')
    plt.savefig(f"{fname_base}.jpg", format='jpg', dpi=300, bbox_inches='tight')

    plt.show()

    print(f"Saved SVG: {fname_base}.svg")

# --- Main loop over CAM scenarios ---
for scenario in scenarios:
    filtered_down = filter_data(df, "down", scenario)
    filtered_up = filter_data(df, "up", scenario)
    plot_bubble_chart(filtered_down, filtered_up, scenario["title"], scenario["q_value"], scenario["flag"])
