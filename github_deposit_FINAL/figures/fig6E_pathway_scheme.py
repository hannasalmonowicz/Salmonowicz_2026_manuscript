# =====================================================================
# Fig 6E -- CAM track pathway scheme: NEAA + TCA + purine & pyrimidine (cata)bolism
# Node colour = direction in Sen_CAM vs Sen_CTRL (facility BMIS stats).
#   red  = higher in Sen_CAM (accumulates),  blue = lower (depletes)
#   * = FDR (adjusted p) < 0.05;  full colour = FDR<0.05, pale = trend
#
# For each metabolite, values are reported in whichever ionization mode
# (negative/positive) has the lower %CV in QC -- computed from the pooled QC
# samples, so it doesn't depend on which two groups are being compared, making
# it an outcome-independent measure of technical precision. A mode with %CV in
# QC >= 30% (the facility's own stated reliability cutoff) is excluded from
# consideration. Every node's full neg/pos t / FDR / %CV and which mode was
# used is recorded in AUDIT and printed/saved, so the choice is never silent.
#
# Input: EV_Table_X_Metabolomics_Statistical_Results.xlsx (or whatever EV
# number is ultimately assigned), which combines the facility's 10 separate
# statistical-results exports. Every value is read verbatim from the
# facility's own numbers -- nothing here recomputes a statistic.
#
# Put that file in a folder named "data/" next to this notebook:
#   EV_Table_X_Metabolomics_Statistical_Results.xlsx
#   (two sheets: Negative_mode_statistics / Positive_mode_statistics,
#    each with a 'Contrast' column -- this script filters both sheets to
#    Contrast == "Sen_CAM vs Sen_CTRL", the only contrast this figure uses.)
#
# Requires: pandas, numpy, matplotlib, openpyxl
# =====================================================================
import pandas as pd, numpy as np, glob
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
DATA = "data"          # folder holding the combined EV stats table
CV_MAX = 30.0           # facility's own QC-reliability cutoff (%CV in QC)
CONTRAST = "Sen_CAM vs Sen_CTRL"   # the only contrast this figure draws


def find_file(patterns, label):
    """Find a required input by glob pattern(s), checking ./data/ first, then the
    current folder. Raises loudly and by name if not found -- never silently skips
    a required input."""
    for pat in patterns:
        matches = glob.glob(f"{DATA}/{pat}") or glob.glob(pat)
        if matches:
            return matches[0]
    raise FileNotFoundError(
        f"Missing required input for {label}: none of {patterns} found in "
        f"./{DATA}/ or next to this script. Check the filename/location and rerun."
    )


EV_STATS_PATH = find_file(
    ["EV_Table_X_Metabolomics_Statistical_Results*.xlsx",
     "EV_Table_*_Metabolomics_Statistical_Results*.xlsx"],
    "combined metabolomics statistical-results EV table "
    "(Negative_mode_statistics / Positive_mode_statistics sheets)",
)
print(f"Using EV_STATS_PATH: {EV_STATS_PATH}\n")

# ---------- 1. LOAD THE Sen_CAM vs Sen_CTRL STAT TABLES (neg + pos) ----------
def load_stat_sheet(sheet_name, contrast):
    df = pd.read_excel(EV_STATS_PATH, sheet_name=sheet_name)
    df.columns = [str(c).strip() for c in df.columns]
    if "Contrast" not in df.columns:
        raise KeyError(
            f"Sheet '{sheet_name}' of {EV_STATS_PATH} has no 'Contrast' column -- "
            f"columns present: {list(df.columns)}. Check the file by hand."
        )
    df = df[df["Contrast"] == contrast].copy()
    if df.empty:
        raise KeyError(
            f"No rows with Contrast == '{contrast}' in sheet '{sheet_name}' of "
            f"{EV_STATS_PATH}. Contrasts present: "
            f"{sorted(pd.read_excel(EV_STATS_PATH, sheet_name=sheet_name)['Contrast'].unique())}"
        )
    df["is_IS"] = df["Identification"].astype(str).str.startswith("IS:")
    df["name"]  = df["Identification"].astype(str).str.rstrip("*").str.strip()
    return df[~df["is_IS"]]


STAT = {
    "neg": load_stat_sheet("Negative_mode_statistics", CONTRAST),
    "pos": load_stat_sheet("Positive_mode_statistics", CONTRAST),
}
# ---------- 2. PER-COMPOUND DIRECTION + FDR (FDR-priority, CV-gated) ----------
AUDIT = []  # collects full per-compound record for the printed audit table

def _mode_record(mode, cmpd):
    """Return dict with t/fdr/cv for cmpd in this polarity, or None if absent/duplicated."""
    r = STAT[mode][STAT[mode]["name"] == cmpd]
    if len(r) == 0:
        return None
    if len(r) > 1:
        # Duplicate stripped-name collision in this polarity -- surface it explicitly
        # rather than silently taking row 0. None of the compounds in this scheme
        # currently hit this, but guard against it if the compound list changes.
        return {"dup": True, "n_rows": len(r)}
    cv  = pd.to_numeric(r["%CV in QC"], errors="coerce").values[0]
    t   = float(r["Statistical Value"].values[0])
    fdr = float(r["adjusted p-value (fdr)"].values[0])
    return {"t": t, "fdr": fdr, "cv": cv, "dup": False}

def node_stat(cmpd):
    """Return (t, fdr) for a compound. Polarity is selected ONLY by %CV in QC
    (the original rule) -- this is an outcome-independent, technical-precision
    metric (computed from pooled QC samples, not from the compared groups), so
    it cannot bias toward significance the way an FDR-based tie-break would.
    The only change from the original logic is an explicit floor: a polarity
    with %CV in QC >= 30% (the facility's own stated reliability cutoff) is
    excluded from consideration rather than silently allowed to win on a technicality.
    Also appends a full audit record to AUDIT for later printing/verification."""
    neg = _mode_record("neg", cmpd)
    pos = _mode_record("pos", cmpd)

    candidates = {}
    for mode, rec in (("neg", neg), ("pos", pos)):
        if rec is None or rec.get("dup"):
            continue
        if rec["cv"] is not None and not np.isnan(rec["cv"]) and rec["cv"] >= CV_MAX:
            continue  # fails facility's own QC-reliability gate
        candidates[mode] = rec

    chosen_mode, reason = None, ""
    if not candidates:
        chosen_mode = None
        reason = "no usable polarity (missing, duplicated, or CV>=30%)"
    elif len(candidates) == 1:
        chosen_mode = next(iter(candidates))
        reason = f"only {chosen_mode} mode available"
    else:
        chosen_mode = min(candidates, key=lambda m: candidates[m]["cv"])
        reason = f"{chosen_mode} has lower %CV in QC ({candidates[chosen_mode]['cv']} vs {candidates['pos' if chosen_mode=='neg' else 'neg']['cv']})"

    AUDIT.append({
        "compound": cmpd,
        "neg_t": None if neg is None or neg.get("dup") else neg["t"],
        "neg_fdr": None if neg is None or neg.get("dup") else neg["fdr"],
        "neg_cv": None if neg is None or neg.get("dup") else neg["cv"],
        "pos_t": None if pos is None or pos.get("dup") else pos["t"],
        "pos_fdr": None if pos is None or pos.get("dup") else pos["fdr"],
        "pos_cv": None if pos is None or pos.get("dup") else pos["cv"],
        "mode_used": chosen_mode,
        "reason": reason,
    })

    if chosen_mode is None:
        return None
    rec = candidates[chosen_mode]
    return (rec["t"], rec["fdr"])

def color_for(cmpd):
    """Returns (face_color, edge_color, arrow_label, detail_label).
    detail_label is the exact FDR value and which polarity it came from --
    read straight from the facility file, not derived/rounded for display
    beyond 3 decimals."""
    s = node_stat(cmpd)
    if s is None:
        return "#e8e8e8", "0.5", "n/m", "", "0.4"
    t, fdr = s
    mode = AUDIT[-1]["mode_used"]  # node_stat() just appended this compound's record
    dark = fdr < 0.05   # dark, saturated fill -> needs light text for contrast
    if t > 0:                                            # accumulates -> red
        face = "#c6161d" if dark else ("#f4a582" if fdr < 0.15 else "#fddbc7")
    else:                                                # depletes -> blue
        face = "#2166ac" if dark else ("#92c5de" if fdr < 0.15 else "#d1e5f0")
    lab = ("↑" if t > 0 else "↓") + ("*" if dark else "")
    detail = f"FDR={fdr:.3f} ({mode})"
    detail_color = "#f5d5d3" if dark else "0.25"  # pale tint readable on dark fill, dark gray on pale fill
    return face, "k", lab, detail, detail_color
# ---------- 3. DRAWING PRIMITIVES ----------
fig, ax = plt.subplots(figsize=(15, 11))
ax.set_xlim(0, 100); ax.set_ylim(0, 74); ax.axis("off")
def node(x, y, cmpd, short=None, w=13, h=5.2, fs=8):
    face, ec, lab, detail, detail_color = color_for(cmpd)
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h,
                 boxstyle="round,pad=0.3,rounding_size=0.8",
                 fc=face, ec=ec, lw=1.0, zorder=3))
    text_main = "w" if detail_color.startswith("#f5") else "k"
    ax.text(x, y + 1.1, (short or cmpd), ha="center", va="center",
            fontsize=fs, zorder=4, fontweight="bold", color=text_main)
    if lab != "n/m":
        ax.text(x, y - 0.6, lab, ha="center", va="center",
                fontsize=fs - 1, zorder=4, color=text_main)
        ax.text(x, y - 2.0, detail, ha="center", va="center",
                fontsize=fs - 3.3, zorder=4, color=detail_color)
def arrow(p1, p2, style="-|>", color="0.35", lw=1.4, ls="-", rad=0.0, shr=8):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=13,
                 color=color, lw=lw, ls=ls, shrinkA=shr, shrinkB=shr, zorder=2,
                 connectionstyle=f"arc3,rad={rad}"))
# ---------- 4. TITLE ----------
ax.text(50, 72, "NEAA, TCA cycle and purine / pyrimidine (cata)bolism under chloramphenicol",
        ha="center", fontsize=12.5, fontweight="bold")
ax.text(50, 69.4, "Node colour = direction in senescent cells (Sen CAM vs Sen CTRL); "
        "red ↑ accumulates, blue ↓ depletes, * FDR<0.05",
        ha="center", fontsize=8.5, color="0.3")
ax.text(50, 67.6, "FDR = facility's adjusted p-value, read directly from the combined EV stats table (not recalculated here). "
        "Each node also states which polarity (neg/pos) it came from -- picked by lower %CV in QC.",
        ha="center", fontsize=7.5, color="0.4", style="italic")
# ---------- 5. NON-ESSENTIAL AMINO ACIDS (left column) ----------
ax.text(15, 65, "Non-essential amino acids", ha="center", fontsize=10, fontweight="bold", color="#08519c")
node(15, 60, "Aspartic acid", "Aspartate")
node(15, 54, "Glutamine", "Glutamine")
node(15, 48, "Glycine", "Glycine")
node(15, 42, "Serine", "Serine")
arrow((15, 44.1), (15, 46), style="<|-|>", color="0.5")
# ---------- 6. PENTOSE-PHOSPHATE PATHWAY (left, mid) ----------
ax.text(15, 37, "Pentose-phosphate pathway", ha="center", fontsize=9.5, fontweight="bold", color="#238b45")
node(15, 32, "Ribose 5-phosphate", "Ribose-5-P")
arrow((18, 32), (43, 55.5), color="#238b45", rad=-0.3, lw=1.6)
ax.text(33, 43, "PRPP ribose\nbackbone", fontsize=7, color="#238b45", ha="center", style="italic")
# ---------- 7. PURINE SYNTHESIS (centre top) ----------
ax.text(45, 65, "Purine synthesis", ha="center", fontsize=10, fontweight="bold", color="#a50f15")
node(45, 57, "Inosine Monophosphate", "IMP", w=11)
node(45, 51, "Adenosine monophosphate", "AMP", w=11)
node(58, 51, "Guanosine monophosphate", "GMP", w=11)
arrow((21.5, 48), (39.5, 56.2), color="#08519c", rad=-0.15, lw=1.6)
arrow((21.5, 60), (39.5, 57.8), color="#08519c", rad=0.05, lw=1.6)
arrow((21.5, 54), (39.5, 57.0), color="#08519c", rad=-0.05, lw=1.6)
ax.text(30, 60.5, "Gly, Asp, Gln\n(N & C donors)", fontsize=7, color="#08519c", ha="center", style="italic")
arrow((45, 54.9), (45, 53.1))
arrow((47, 55.5), (56, 53))
# ---------- 8. PURINE CATABOLISM (right top) ----------
ax.text(80, 65, "Purine catabolism", ha="center", fontsize=10, fontweight="bold", color="#a50f15")
node(80, 59, "Inosine", "Inosine", w=11)
node(80, 53, "Hypoxanthine", "Hypoxanthine", w=13)
node(80, 47, "Xanthine", "Xanthine", w=11)
node(80, 41, "Allantoin", "Allantoin", w=11)
node(93, 53, "Guanosine", "Guanosine", w=11)
node(93, 47, "Guanine", "Guanine", w=10)
arrow((50.5, 57), (74.5, 59), color="0.35", rad=-0.1)
arrow((80, 56.9), (80, 55.1)); arrow((80, 50.9), (80, 49.1)); arrow((80, 44.9), (80, 43.1))
arrow((63.5, 51), (88, 53), color="0.35", rad=-0.05)
arrow((93, 50.9), (93, 49.1)); arrow((90, 46), (83, 48), color="0.35", rad=0.1)
# ---------- 9. PYRIMIDINE SYNTHESIS & CATABOLISM (centre band) ----------
ax.text(62, 37, "Pyrimidine synthesis & catabolism", ha="center", fontsize=9.5, fontweight="bold", color="#6a51a3")
node(40, 31, "UMP", "UMP (n.m.)", w=13)
node(56, 31, "Uridine", "Uridine", w=11)
node(72, 31, "Uracil", "Uracil", w=10)
node(88, 31, "Beta-Alanine", "β-Alanine", w=11)
node(56, 24, "Cytidine", "Cytidine", w=11)
arrow((17, 40), (34, 33), color="#08519c", rad=-0.25, lw=1.6)
ax.text(27, 47, "Asp, Gln\n(ring atoms)", fontsize=7, color="#08519c", ha="center", style="italic")
arrow((46.5, 31), (50.5, 31)); arrow((61.5, 31), (67, 31)); arrow((77, 31), (82.5, 31))
arrow((54, 28.9), (54, 26.1), style="-|>")
arrow((58, 26.1), (58, 28.9), style="-|>", color="0.6")
# salvage feedback (purine ribose -> R5P)
arrow((80, 60.8), (18, 33.5), color="#238b45", ls="--", rad=0.34, lw=1.2)
ax.text(55, 49, "ribose salvage → R5P", fontsize=7, color="#238b45", ha="center", style="italic", rotation=9)
# ---------- 10. TCA CYCLE (bottom strip) ----------
ax.text(52, 18.5, "TCA cycle", ha="center", fontsize=10, fontweight="bold", color="#b35806")
tca_y = 12
node(34, tca_y, "Malic acid",       "Malate",        w=11)
node(48, tca_y, "Fumaric acid",     "Fumarate",      w=11)
node(62, tca_y, "Succinic acid",    "Succinate",     w=11)
node(76, tca_y, "Oxoglutaric acid", "α-KG",     w=11)
node(90, tca_y, "cis-Aconitic acid","cis-Aconitate", w=12)
# oxidative-direction chain: aconitate -> a-KG -> succinate -> fumarate -> malate
arrow((84, tca_y), (81.5, tca_y)); arrow((70.5, tca_y), (67.5, tca_y))
arrow((56.5, tca_y), (53.5, tca_y)); arrow((42.5, tca_y), (39.5, tca_y))
# return to OAA (not measured) -> feeds aspartate
arrow((30, tca_y+1.2), (26, tca_y+3.0), color="#b35806", rad=0.2, lw=1.2)
ax.text(24, tca_y+4.2, "→ OAA", fontsize=7, color="#b35806", ha="center", style="italic")
# anaplerosis: NEAA <-> TCA carbon skeletons
arrow((15, 57.9), (32, 14.5), color="#e08214", rad=-0.22, lw=1.6, style="<|-|>")
ax.text(35, 25, "Asp ↔ OAA/malate", fontsize=7, color="#e08214", ha="center", style="italic")
arrow((16, 52), (74, 14.6), color="#e08214", rad=0.30, lw=1.4, style="<|-|>")
ax.text(60, 20.5, "Glu/Gln ↔ α-KG", fontsize=7, color="#e08214", ha="center", style="italic")
# ---------- 11. LEGEND ----------
leg = [
    Line2D([0],[0], marker="s", ls="", mfc="#c6161d", mec="k", ms=12, label="↑ accumulates (FDR<0.05)"),
    Line2D([0],[0], marker="s", ls="", mfc="#fddbc7", mec="k", ms=12, label="↑ trend"),
    Line2D([0],[0], marker="s", ls="", mfc="#2166ac", mec="k", ms=12, label="↓ depletes (FDR<0.05)"),
    Line2D([0],[0], marker="s", ls="", mfc="#d1e5f0", mec="k", ms=12, label="↓ trend"),
    Line2D([0],[0], marker="s", ls="", mfc="#e8e8e8", mec="0.5", ms=12, label="not measured (n.m.)"),
    Line2D([0],[0], color="#08519c", lw=2, label="NEAA nitrogen/carbon donation"),
    Line2D([0],[0], color="#238b45", lw=2, label="PPP ribose contribution"),
    Line2D([0],[0], color="#e08214", lw=2, label="NEAA ↔ TCA (anaplerosis)"),
]
ax.legend(handles=leg, loc="lower center", ncol=4, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.03))
fig.savefig("CAM_pathway_scheme_FIXED.png", dpi=200, bbox_inches="tight")
fig.savefig("CAM_pathway_scheme_FIXED.svg", bbox_inches="tight")
plt.close(fig)
# ---------- 12. FULL AUDIT TABLE (both polarities + which was chosen + why) ----------
allc = (["Aspartic acid","Glutamine","Glycine","Serine","Ribose 5-phosphate"] +
        ["Inosine Monophosphate","Adenosine monophosphate","Guanosine monophosphate"] +
        ["Inosine","Hypoxanthine","Xanthine","Allantoin","Guanosine","Guanine"] +
        ["Uridine","Uracil","Beta-Alanine","Cytidine"] +
        ["Malic acid","Fumaric acid","Succinic acid","Oxoglutaric acid","cis-Aconitic acid"])
AUDIT = []  # reset (node() calls above already populated it once during drawing; redo cleanly)
for c in allc:
    node_stat(c)  # appends the full neg/pos/mode_used/reason record to AUDIT
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 220)
audit_df = pd.DataFrame(AUDIT)
print(audit_df.to_string(index=False))
audit_df.to_csv("CAM_pathway_scheme_audit.csv", index=False)
print("\nSaved: CAM_pathway_scheme_FIXED.png/.svg and CAM_pathway_scheme_audit.csv")