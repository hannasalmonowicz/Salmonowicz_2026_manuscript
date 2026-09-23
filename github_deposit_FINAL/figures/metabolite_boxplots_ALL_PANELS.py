"""
ALL MANUSCRIPT METABOLITE BOXPLOTS (Fig 2G/H, 4C-G, 5C, 6C/D, EV6G/H,
EV10F/G-H -- 21 panels).

Inputs -- both are EV tables, nothing else needed:

(1) EV_Table_3*.xlsx -- sheets 'Neg_mode_Normalized-BMIS' / 'Pos_mode_Normalized-BMIS',
    standard MS-DIAL layout (a 'Class' label cell marking where the per-sample columns
    start, a 'Metabolite name' column, sample columns grouped by Class into Prolif_CAM,
    Prolif_CTRL, Prolif_GAL, Sen_CAM, Sen_CTRL, Sen_GAL, QC). The parser below locates
    both by searching the sheet rather than assuming a fixed row/column.

(2) EV_Table_11_Metabolomics_Statistical_Results*.xlsx -- two sheets,
    'Negative_mode_statistics' / 'Positive_mode_statistics', each with a 'Contrast'
    column ("<A> vs <B>") and an 'adjusted p-value (fdr)' column -- carries the FDR
    used for the significance bracket on every single-metabolite panel (ratio panels
    keep a directly computed Welch's t-test, unchanged -- see STATISTICS below).
    Optional: if this file isn't present, single-metabolite panels still render from
    the real EV_Table_3 data, with the bracket labelled 'FDR pending' instead of a
    value.

STATISTICS: single-metabolite comparisons come from the real facility FDR, never
recomputed. Ratio panels (e.g. "NAD / NADH") are the exception and keep a Welch's
t-test computed directly on the log2 ratio.

COLORS: CAM_PAL/CAM_PT are the final published colors. GAL_PAL/GAL_PT and
BASE_PAL/BASE_PT/SEN_GAL_PAL/SEN_GAL_PT are still placeholders -- swap in the real
GAL-track/baseline hex codes at the palette definitions below.
"""

import glob
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

# ---------- 1. CONFIGURATION ----------
DATA_DIR = "data"
ENC = "latin-1"
OUT = "panels"
os.makedirs(OUT, exist_ok=True)


def find_file(patterns, label):
    """Find a required input by glob pattern(s), checking ./data/ first, then
    the current folder, and raises FileNotFoundError naming the missing
    file."""
    for pat in patterns:
        matches = glob.glob(f"{DATA_DIR}/{pat}") or glob.glob(pat)
        if matches:
            return matches[0]
    raise FileNotFoundError(
        f"Missing required input for {label}: none of {patterns} found in "
        f"./{DATA_DIR}/ or next to this script. Check the filename/location and rerun."
    )


EV_TABLE_3_PATH = find_file(
    ["EV_Table_3*.xlsx"],
    "untargeted metabolomics BMIS-normalized quantification matrix "
    "(Neg_mode_Normalized-BMIS / Pos_mode_Normalized-BMIS sheets)",
)
print(f"Using EV_TABLE_3_PATH: {EV_TABLE_3_PATH}\n")

# Optional: EV_Table_11 (significance brackets). Not required to run the script --
# if absent, brackets just read "FDR pending". Matched by a number-agnostic glob,
# same convention as every other multi-digit EV table in this repo, so it's found
# regardless of what number it ends up as (or a placeholder like "EV_Table_X_...").
_ev11_matches = glob.glob(f"{DATA_DIR}/EV_Table_*_Metabolomics_Statistical_Results*.xlsx") \
    or glob.glob("EV_Table_*_Metabolomics_Statistical_Results*.xlsx")
EV_STATS_PATH = _ev11_matches[0] if _ev11_matches else None
print(f"Using EV_STATS_PATH: {EV_STATS_PATH}"
      + ("" if EV_STATS_PATH else "  (not found -- brackets will read 'FDR pending')") + "\n")

# ---------- 2. BMIS MATRIX PARSER (EV_Table_3, xlsx) ----------
def parse_bmis_ev_table(path, sheet_name):
    """Parses one BMIS-normalized quantification sheet of the EV table.
    Standard MS-DIAL layout (a 'Class' label cell marking where per-sample
    columns start, a 'Metabolite name' column labelling each compound row),
    located by searching the sheet rather than assuming a fixed row/column."""
    raw = pd.read_excel(path, sheet_name=sheet_name, header=None, dtype=str)

    class_hits = [
        (r, c) for r in range(min(10, raw.shape[0])) for c in range(raw.shape[1])
        if str(raw.iat[r, c]).strip() == "Class"
    ]
    if not class_hits:
        raise KeyError(
            f"Could not find a 'Class' label in the first 10 rows of "
            f"{path} [{sheet_name}]. Sheet layout may have changed -- check by hand."
        )
    class_row, lblcol = class_hits[0]

    name_hits = [
        (r, c) for r in range(min(10, raw.shape[0])) for c in range(raw.shape[1])
        if str(raw.iat[r, c]).strip() == "Metabolite name"
    ]
    if not name_hits:
        raise KeyError(
            f"Could not find a 'Metabolite name' column in the first 10 rows of "
            f"{path} [{sheet_name}]. Sheet layout may have changed -- check by hand."
        )
    hdr_row, namecol = name_hits[0]

    scols = list(range(lblcol + 1, raw.shape[1]))
    classes = raw.iloc[class_row, scols].tolist()
    # MS-DIAL/BMIS marks an ambiguous ID with a trailing '*' -- strip it so
    # the name matches what TARGETS / NAME_ALIASES below actually type.
    names = raw.iloc[hdr_row + 1:, namecol].map(
        lambda s: str(s).strip().rstrip("*").strip()
    ).values
    vals = raw.iloc[hdr_row + 1:, scols].apply(pd.to_numeric, errors="coerce")

    vals.index = names
    vals.columns = [f"S{i}" for i in range(len(scols))]
    return vals, np.array(classes)


neg_v, neg_cls = parse_bmis_ev_table(EV_TABLE_3_PATH, "Neg_mode_Normalized-BMIS")
pos_v, pos_cls = parse_bmis_ev_table(EV_TABLE_3_PATH, "Pos_mode_Normalized-BMIS")


def qc_cv(vals, cls):
    qc = vals.loc[:, cls == "QC"]
    return (qc.std(axis=1, ddof=1) / qc.mean(axis=1)) * 100


neg_cv, pos_cv = qc_cv(neg_v, neg_cls), qc_cv(pos_v, pos_cls)

# ---------- 3. NAME RESOLUTION (unchanged) ----------
NAME_ALIASES = {
    "2-Oxoglutarate": ["Oxoglutaric acid", "2-Oxoglutaric acid",
                        "2-Oxoglutarate", "Alpha-ketoglutaric acid",
                        "alpha-Ketoglutarate"],
    "Succinate": ["Succinic acid", "Succinate"],
    "Fumarate": ["Fumaric acid", "Fumarate"],
    "Oxoglutaric acid": ["2-Oxoglutaric acid", "alpha-Ketoglutarate",
                          "Alpha-ketoglutaric acid"],
    "Lactic acid": ["Lactate"],
    "Aspartic acid": ["Aspartate"],
    "Malic acid": ["Malate"],
    "Glutamic acid": ["Glutamate"],
}


def resolve_name(name, df):
    """Finds exact string match in dataframe index using aliases."""
    if name in df.index:
        return name
    for alias in NAME_ALIASES.get(name, []):
        if alias in df.index:
            return alias
    return None


def pick_mode(compound_name):
    neg_key = resolve_name(compound_name, neg_v)
    pos_key = resolve_name(compound_name, pos_v)
    if neg_key is None and pos_key is None:
        raise KeyError(
            f"Compound '{compound_name}' not found in either polarity. "
            f"Check spelling / add an alias in NAME_ALIASES."
        )
    n = neg_cv.get(neg_key, np.inf) if neg_key else np.inf
    p = pos_cv.get(pos_key, np.inf) if pos_key else np.inf
    return "neg" if n <= p else "pos"


def get_series(compound_name):
    """Returns log2 BMIS values and class labels. 'A / B' -> ratio via
    log2(A/B) = log2(A) - log2(B), same as before."""
    if " / " in compound_name:
        num_name, den_name = compound_name.split(" / ")
        s_num, cls = get_series(num_name)
        s_den, _ = get_series(den_name)
        return (s_num - s_den), cls

    mode = pick_mode(compound_name)
    v, c = (neg_v, neg_cls) if mode == "neg" else (pos_v, pos_cls)
    key = resolve_name(compound_name, v)
    if key is None:
        raise KeyError(f"Compound '{compound_name}' not found in {mode} mode BMIS.")
    return np.log2(v.loc[key]), c


# ---------- 4. STATISTICAL ENGINE ----------
# Ratio panels: Welch's t-test computed directly on the log2 ratio, as
# before -- this is the confirmed exception to the "stats come from the
# facility file" rule.
def compute_welch_p(data_dict, cond1, cond2):
    v1, v2 = data_dict[cond1], data_dict[cond2]
    if len(v1) < 2 or len(v2) < 2:
        return np.nan
    _, p = ttest_ind(v1, v2, equal_var=False, nan_policy="omit")
    return p


def fmt_p(p):
    if pd.isna(p):
        return "n.d."
    if p < 0.001:
        return "p<0.001"
    return f"p={p:.3f}"


# Single-metabolite panels: FDR looked up from EV_Table_11, for the exact
# contrast being drawn, in whichever polarity pick_mode() already selected for
# that compound -- never recomputed. STATS_MISSING collects every (contrast,
# polarity) this run actually needed but could not find (file absent, or that
# contrast/compound not in it), so the end-of-run summary says precisely what's
# still pending rather than a blanket "supply everything" ask.
_EV_STATS_SHEET_CACHE = {}
STATS_MISSING = set()


def _load_ev_stats_sheet(polarity):
    """Loads (and caches) one sheet of EV_Table_11. Returns None if the file
    isn't present -- callers treat that as 'not yet supplied', same as before."""
    if EV_STATS_PATH is None:
        return None
    if polarity in _EV_STATS_SHEET_CACHE:
        return _EV_STATS_SHEET_CACHE[polarity]
    sheet = "Negative_mode_statistics" if polarity == "negative" else "Positive_mode_statistics"
    df = pd.read_excel(EV_STATS_PATH, sheet_name=sheet)
    df.columns = [str(c).strip() for c in df.columns]
    df["_clean_id"] = (
        df["Identification"].astype(str).str.strip().str.rstrip("*").str.strip()
    )
    df = df[~df["_clean_id"].str.startswith("IS:")]
    _EV_STATS_SHEET_CACHE[polarity] = df
    return df


def _load_stats_file(cond_a, cond_b, polarity):
    df = _load_ev_stats_sheet(polarity)
    if df is None:
        return None
    for contrast in (f"{cond_a} vs {cond_b}", f"{cond_b} vs {cond_a}"):
        sub = df[df["Contrast"] == contrast]
        if not sub.empty:
            return sub
    return None


def fdr_lookup(compound_name, cond_a, cond_b):
    """Looks up 'adjusted p-value (fdr)' for compound_name in EV_Table_11,
    for the cond_a/cond_b contrast, in the polarity pick_mode() already
    selected for this compound. Returns (value, found) -- value is NaN and
    found is False when the file/contrast isn't available, which the caller
    renders as 'FDR pending'."""
    mode = pick_mode(compound_name)
    stats_df = _load_stats_file(cond_a, cond_b, "negative" if mode == "neg" else "positive")
    if stats_df is None:
        STATS_MISSING.add((cond_a, cond_b, "negative" if mode == "neg" else "positive"))
        return np.nan, False
    key = resolve_name(compound_name, neg_v if mode == "neg" else pos_v)
    row = stats_df[stats_df["_clean_id"] == key]
    if row.empty:
        return np.nan, False
    return float(row["adjusted p-value (fdr)"].iloc[0]), True


def fmt_fdr(fdr, found):
    if not found:
        return "FDR pending"
    if pd.isna(fdr):
        return "n.d."
    if fdr < 0.001:
        return "FDR<0.001"
    return f"FDR={fdr:.3f}"


# ---------- 5. STYLE ----------
plt.rcdefaults()
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "svg.fonttype": "none",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.0,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "axes.grid": False,
})

# CAM_PAL/CAM_PT: final published colors.
CAM_PAL = {"Prolif_CTRL": "#f3b98a", "Prolif_CAM": "#f3b98a",
           "Sen_CTRL": "#9ba4b0", "Sen_CAM": "#e2453b"}
CAM_PT = {"Prolif_CTRL": "#e08a3c", "Prolif_CAM": "#e08a3c",
          "Sen_CTRL": "#6f6f6f", "Sen_CAM": "#b02a22"}

# PLACEHOLDER -- same style logic (Prolif = orange family, Sen_CTRL = grey,
# treatment = a second color), not the final GAL-track hex codes.
GAL_PAL = {"Prolif_CTRL": "#f3b98a", "Prolif_GAL": "#f3b98a",
           "Sen_CTRL": "#9ba4b0", "Sen_GAL": "#7fb3d5"}
GAL_PT = {"Prolif_CTRL": "#e08a3c", "Prolif_GAL": "#e08a3c",
          "Sen_CTRL": "#6f6f6f", "Sen_GAL": "#3d7ab0"}

# PLACEHOLDER -- 2-condition version of the same family.
BASE_PAL = {"Prolif_CTRL": "#f3b98a", "Sen_CTRL": "#9ba4b0"}
BASE_PT = {"Prolif_CTRL": "#e08a3c", "Sen_CTRL": "#6f6f6f"}

SEN_GAL_PAL = {"Sen_CTRL": "#9ba4b0", "Sen_GAL": "#7fb3d5"}
SEN_GAL_PT = {"Sen_CTRL": "#6f6f6f", "Sen_GAL": "#3d7ab0"}


# ---------- 6. GENERALIZED PLOTTING ENGINE ----------
def plot_box(ax, compound, title, conditions, pal, pt, pairs, base_contrast=None):
    S, cls = get_series(compound)
    XS = np.arange(len(conditions))
    data_dict = {c: S.values[cls == c] for c in conditions}
    data_list = [data_dict[c] for c in conditions]

    is_ratio = " / " in compound

    def stat_label(c_left, c_right):
        if is_ratio:
            val = compute_welch_p(data_dict, c_left, c_right)
            return val, fmt_p(val)
        val, found = fdr_lookup(compound, c_left, c_right)
        return val, fmt_fdr(val, found)

    bp = ax.boxplot(
        data_list, positions=XS, widths=0.45, patch_artist=True,
        medianprops=dict(color="k", lw=1.5),
        whiskerprops=dict(color="k", lw=0.9),
        capprops=dict(color="k", lw=0.9),
        boxprops=dict(lw=0.0), showfliers=False,
    )
    for patch, c in zip(bp["boxes"], conditions):
        patch.set_facecolor(pal[c])

    rng = np.random.RandomState(42)
    for i, c in enumerate(conditions):
        y = data_dict[c]
        x_pos = np.full(len(y), i) + rng.uniform(-0.10, 0.10, len(y))
        ax.scatter(x_pos, y, s=16, c=pt[c], edgecolor="k", lw=0.4, zorder=5)

    ymax = max(np.max(d) for d in data_list if len(d) > 0)
    ymin = min(np.min(d) for d in data_list if len(d) > 0)
    span = ymax - ymin + 1e-9

    y0 = ymax + span * 0.08
    for (i0, i1), (c_left, c_right) in pairs:
        _, label = stat_label(c_left, c_right)
        ax.plot([i0, i0, i1, i1], [y0, y0 + span * 0.03, y0 + span * 0.03, y0],
                lw=0.9, c="k")
        ax.text((i0 + i1) / 2, y0 + span * 0.045, label,
                ha="center", va="bottom", fontsize=7.5)

    top_y = y0
    if base_contrast is not None:
        c_left, c_right = base_contrast
        i0, i1 = conditions.index(c_left), conditions.index(c_right)
        val, label = stat_label(c_left, c_right)
        yb = ymax + span * 0.30
        is_ns = (not pd.isna(val)) and val >= 0.05
        lab = "ns" if is_ns else label
        ax.plot([i0, i0, i1, i1], [yb, yb + span * 0.03, yb + span * 0.03, yb],
                lw=0.9, c="0.55", ls=(0, (4, 2)))
        ax.text((i0 + i1) / 2, yb + span * 0.045, lab, ha="center", va="bottom",
                fontsize=7.5, color=("0.4" if is_ns else "k"),
                style=("italic" if is_ns else "normal"))
        top_y = yb

    ax.set_ylim(ymin - span * 0.10, top_y + span * 0.22)
    ax.set_xticks(XS)
    ax.set_xticklabels([c.replace("_", "\n") for c in conditions], fontsize=7.5)
    ax.set_title(title, fontsize=10.0, fontstyle="italic", fontweight="bold", pad=5)
    ax.set_ylabel("Log2 Ratio" if " / " in compound else "Log2 BMIS", fontsize=8.0)

    primary = compound if " / " not in compound else compound.split(" / ")[0]
    mode_tag = "positive mode" if pick_mode(primary) == "pos" else "negative mode"
    ax.text(0.02, 0.98, mode_tag, transform=ax.transAxes, fontsize=7.0,
            color="#3b5aa8", fontstyle="italic", va="top", ha="left")
    ax.tick_params(labelsize=7.5)
    ax.grid(False)


# ---------- 7. SHAPE PRESETS ----------
def plot_cam_track(ax, compound, title):
    conditions = ["Prolif_CTRL", "Prolif_CAM", "Sen_CTRL", "Sen_CAM"]
    pairs = [((0, 1), ("Prolif_CTRL", "Prolif_CAM")),
             ((2, 3), ("Sen_CTRL", "Sen_CAM"))]
    plot_box(ax, compound, title, conditions, CAM_PAL, CAM_PT, pairs,
             base_contrast=("Prolif_CTRL", "Sen_CTRL"))


def plot_gal_track(ax, compound, title):
    conditions = ["Prolif_CTRL", "Prolif_GAL", "Sen_CTRL", "Sen_GAL"]
    pairs = [((0, 1), ("Prolif_CTRL", "Prolif_GAL")),
             ((2, 3), ("Sen_CTRL", "Sen_GAL"))]
    plot_box(ax, compound, title, conditions, GAL_PAL, GAL_PT, pairs,
             base_contrast=("Prolif_CTRL", "Sen_CTRL"))


def plot_baseline(ax, compound, title):
    conditions = ["Prolif_CTRL", "Sen_CTRL"]
    pairs = [((0, 1), ("Prolif_CTRL", "Sen_CTRL"))]
    plot_box(ax, compound, title, conditions, BASE_PAL, BASE_PT, pairs,
             base_contrast=None)


def plot_sen_only_gal(ax, compound, title):
    conditions = ["Sen_CTRL", "Sen_GAL"]
    pairs = [((0, 1), ("Sen_CTRL", "Sen_GAL"))]
    plot_box(ax, compound, title, conditions, SEN_GAL_PAL, SEN_GAL_PT, pairs,
             base_contrast=None)


SHAPE_FUNCS = {
    "baseline": plot_baseline,
    "sen_only_gal": plot_sen_only_gal,
    "gal_track": plot_gal_track,
    "cam_track": plot_cam_track,
}

# ---------- 8. ALL 21 MANUSCRIPT PANELS ----------
# (fig, shape, compound-as-in-BMIS-file[/ "A / B" for a ratio], display title)
TARGETS = [
    ("Fig 2H", "baseline", "Hypoxanthine", "Hypoxanthine"),
    ("Fig 2G", "baseline", "NAD / NADH", "NAD+/NADH"),
    ("Fig 4C", "sen_only_gal", "Glycerol 3-phosphate", "Glycerol 3-phosphate (G3-P)"),
    ("Fig 4D", "sen_only_gal", "Oxoglutaric acid", "alpha-Ketoglutarate (aKG)"),
    ("Fig 4E", "gal_track", "NAD / NADH", "NAD+/NADH (GAL track)"),
    ("Fig 4F", "gal_track", "Lactic acid", "Lactate"),
    ("Fig 4G", "gal_track", "Aspartic acid / Malic acid", "Aspartate/Malate"),
    ("Fig 5C", "cam_track", "NAD / NADH", "NAD+/NADH (CAM track)"),
    ("Fig 6C", "cam_track", "Glutamic acid", "Glutamate"),
    ("Fig 6C", "cam_track", "Aspartic acid", "Aspartate"),
    ("Fig 6C", "cam_track", "Glycine", "Glycine"),
    ("Fig 6C", "cam_track", "Serine", "Serine"),
    ("Fig 6D", "cam_track", "Cytidine", "Cytidine"),
    ("Fig 6D", "cam_track", "Uracil", "Uracil"),
    ("EV6G", "baseline", "Glycerol 3-phosphate", "Glycerol 3-phosphate (G3-P)"),
    ("EV6H", "baseline", "Oxoglutaric acid", "alpha-Ketoglutarate (aKG)"),
    ("EV10F", "cam_track", "Serine / Glycine", "Serine/Glycine"),
    ("EV10G-H", "cam_track", "Hypoxanthine", "Hypoxanthine"),
    ("EV10G-H", "cam_track", "Inosine", "Inosine"),
    ("EV10G-H", "cam_track", "Guanine", "Guanine"),
    ("EV10G-H", "cam_track", "Guanosine", "Guanosine"),
]

# Not yet added -- legend text doesn't name specific metabolites:
#   EV9B (NEAA), EV9C (BCAA), EV10I (pyrimidines), EV10J (pyrimidine-linked + beta-alanine)

# ---------- 9. RENDER ----------
if __name__ == "__main__":
    rendered, skipped = [], []
    for i, (fig_name, shape, compound, title) in enumerate(TARGETS, 1):
        fn = SHAPE_FUNCS[shape]
        f, ax = plt.subplots(figsize=(3.2, 3.4))
        try:
            fn(ax, compound, f"{fig_name} - {title}")
        except KeyError as e:
            plt.close(f)
            skipped.append((fig_name, title, str(e)))
            continue
        safe = "".join(ch if ch.isalnum() else "_" for ch in title).strip("_")
        out_png = f"{OUT}/{i:02d}_{shape}_{safe}.png"
        f.tight_layout()
        f.savefig(out_png, dpi=300, bbox_inches="tight")
        f.savefig(out_png.replace(".png", ".svg"), bbox_inches="tight")
        plt.close(f)
        rendered.append((fig_name, title, out_png))

    print(f"Rendered {len(rendered)}/{len(TARGETS)} panels to {OUT}/")
    for fig_name, title, path in rendered:
        print(f"  {fig_name:10s} {title:35s} -> {path}")

    if skipped:
        print(f"\n{len(skipped)} panel(s) not rendered -- compound not found in "
              f"EV_Table_3:")
        for fig_name, title, reason in skipped:
            print(f"  {fig_name:10s} {title:35s} {reason}")

    if STATS_MISSING:
        print(
            "\nEV_Table_11 contrast(s) still needed for full transparency "
            "(single-metabolite significance brackets above currently read "
            "'FDR pending' until these are found):"
        )
        for cond_a, cond_b, polarity in sorted(STATS_MISSING):
            print(f"  {cond_a} vs {cond_b}  ({polarity} mode)")
        print(
            "Each panel's boxes and data are already from the real EV_Table_3 "
            "values -- only the significance annotation is pending."
        )
