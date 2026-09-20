#!/usr/bin/env python3
"""
Builds figures/data/GOBP_annotations_CAM_IMT.xlsx -- the per-gene GO
Biological Process annotation that GOenrichment_CAM_vs_CTRL_bubbleplots.py
needs (its GO_Names_P column) and that isn't part of any submitted EV table.

Not a figure-producing script -- kept here with the other verification/
build scripts.

Source: a raw QuickGO annotation export (tab-separated, one row per
(gene product, GO term) pair, all three GO aspects mixed together --
Function/Component/Process), exported from QuickGO on 2025-06-29 (GO release
2025-06-01), filename 'QuickGO-annotations-1751219500456-20250629 (3).tsv'.
Not included in this deposit -- at ~90MB it's a raw external-database dump.
Point QUICKGO_PATH below at your own copy to rebuild.

What this script does, nothing more:
  1. Filters the QuickGO export to GO ASPECT == 'P' (Biological Process)
     and TAXON ID == 9606 (human) -- drops Function/Component rows and the
     11 non-human (Neanderthal, taxon 63221) rows in the raw export.
  2. Collapses to one row per (SYMBOL, GO NAME) pair (some genes carry the
     same GO term from multiple evidence codes/references in the raw
     export -- deduplicated here, nothing recomputed).
  3. For every unique 'Gene names' value actually present in EV_Table_8
     (including the 158 semicolon-joined protein-group entries, e.g.
     'UQCRFS1;UQCRFS1P1'), looks up each constituent gene symbol and unions
     their BP term sets, so the output's 'Gene names' values match
     EV_Table_8's exactly and the two can be pd.merge()'d directly on that
     column with no further key-matching logic needed in the plotting
     script.
  4. Writes 'Gene names' + 'GO_Names_P' (terms sorted, joined with '; ')
     to figures/data/GOBP_annotations_CAM_IMT.xlsx. Genes with zero BP
     terms in the export get GO_Names_P = NaN (dropped downstream by the
     plotting script's own .dropna()).

Verified against the real files: 136,264 raw BP rows -> 102,833
unique (symbol, term) pairs across 19,345 symbols; matched 4,485 of
EV_Table_8's 4,811 unique 'Gene names' values to at least one BP term
(93.2%) -- the remaining 6.8% simply have no Biological Process annotation
in this QuickGO export, not a script bug (spot-checked several: obscure/
uncharacterized ORFs and a few outdated gene symbols).
"""
import pandas as pd

QUICKGO_PATH = 'QuickGO-annotations-1751219500456-20250629 (3).tsv'  # point at your copy
EV8_PATH = '../figures/data/EV_Table_8_LFQ_MS_CAM_exposed_prolif_and_sen.xlsx'
OUT_PATH = '../figures/data/GOBP_annotations_CAM_IMT.xlsx'


def main():
    go = pd.read_csv(QUICKGO_PATH, sep='\t')
    go = go[(go['GO ASPECT'] == 'P') & (go['TAXON ID'] == 9606)]
    go = go[['SYMBOL', 'GO NAME']].dropna().drop_duplicates()
    symbol_to_terms = go.groupby('SYMBOL')['GO NAME'].apply(lambda s: sorted(set(s))).to_dict()
    print(f"QuickGO: {len(go)} unique (symbol, BP term) pairs, "
          f"{len(symbol_to_terms)} symbols with >=1 BP term")

    ev8 = pd.read_excel(EV8_PATH, sheet_name='KSz_Jun2025_set1_Proteins')
    gene_name_values = ev8['Gene names'].dropna().unique()
    print(f"EV_Table_8: {len(gene_name_values)} unique 'Gene names' values to annotate")

    rows = []
    n_matched = 0
    for gn in gene_name_values:
        parts = [p.strip() for p in str(gn).split(';') if p.strip()]
        terms = set()
        for p in parts:
            terms.update(symbol_to_terms.get(p, []))
        if terms:
            n_matched += 1
        rows.append({'Gene names': gn, 'GO_Names_P': '; '.join(sorted(terms)) if terms else None})

    out = pd.DataFrame(rows)
    out.to_excel(OUT_PATH, index=False)
    print(f"Matched {n_matched}/{len(gene_name_values)} 'Gene names' values to >=1 BP term "
          f"({100 * n_matched / len(gene_name_values):.1f}%)")
    print(f"Wrote {OUT_PATH}")


if __name__ == '__main__':
    main()
