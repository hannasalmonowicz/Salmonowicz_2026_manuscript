#!/usr/bin/env python3
"""
Builds the two non-EV-table inputs fig2B_isolated_mito_pathway_ANOVA.py needs:

1. figures/data/MitoCarta_pathway_annotation.xlsx -- pathway/complex annotation
   (columns 'Gene names', 'ComplexI'..'ComplexV', 'Translation',
   'mtDNA_maintenance', 'mtRNA metabolism', 'Lipid_metabolism',
   'MitoCarta3.0_MitoPathways').
2. figures/data/mitoproteome_isolated_mito.txt -- the reference background
   gene list fig2B normalizes against (see below).

Not a figure-producing script -- kept here with the other verification/
build scripts.

Source: the official MitoCarta3.0 human gene list (sheet 'A Human
MitoCarta3.0', 1136 genes), 'Human.MitoCarta3.0.xls' -- the standard
Broad Institute download, not generated for this manuscript.

Why build this instead of reusing 5_Annotations.xlsx's MitoCartaPathways
sheet directly: that sheet only covers the ~4142 genes in the whole-cell
(time-resolved) dataset's own universe, and fig2B is about the isolated-
mitochondria dataset (EV_Table_2) -- a different gene set. Rebuilding from
the full 1136-gene official MitoCarta3.0 list instead gives complete
coverage regardless of which genes EV_Table_2 happens to contain.

mitoproteome_isolated_mito.txt's definition follows the manuscript's EV
Methods, "Mitochondrial proteome data normalization and relative
enrichment", directly: "A curated, high-confidence mitochondrial proteome
annotation (MitoCarta3.0) was filtered to the subset of proteins detected in
this experiment to define a reference background population." fig2B's own
code already does exactly that filtering step in-script (intersects this
file's gene set with its detected 'Gene names') -- so this file is simply the
full MitoCarta3.0 gene symbol list, one gene per line, matching the published
method exactly.

How the boolean columns are derived, and how that's verified as correct:
MitoCarta3.0_MitoPathways is a pipe-separated list of '>'-delimited pathway
hierarchy paths per gene, e.g.
'OXPHOS > Complex I > CI subunits | Metabolism > Metals and cofactors > ...'.
A gene gets ComplexI=1 iff 'Complex I' appears as an exact path SEGMENT
(not a substring match, which would wrongly also flag 'Complex II' /
'Complex III' etc. since they share the 'Complex I' prefix). Same exact-
segment rule for the other 8 categories. This rule was reverse-engineered
from, and checked against, 5_Annotations.xlsx's own MitoCartaPathways
sheet (already-confirmed, already used elsewhere in this project) and
reproduces all 4142 of its rows x 9 boolean columns with ZERO mismatches --
i.e. this is the same categorization already trusted in this project,
just applied to the full official gene list instead of one dataset's
subset of it.
"""
import pandas as pd

MITOCARTA_PATH = 'Human.MitoCarta3.0.xls'  # point at your copy
OUT_PATH = '../figures/data/MitoCarta_pathway_annotation.xlsx'
MITOPROTEOME_OUT_PATH = '../figures/data/mitoproteome_isolated_mito.txt'

CATEGORIES = {
    'ComplexI': 'Complex I', 'ComplexII': 'Complex II', 'ComplexIII': 'Complex III',
    'ComplexIV': 'Complex IV', 'ComplexV': 'Complex V', 'Translation': 'Translation',
    'mtDNA_maintenance': 'mtDNA maintenance', 'mtRNA metabolism': 'mtRNA metabolism',
    'Lipid_metabolism': 'Lipid metabolism',
}


def category_flags(pathway_str):
    out = {c: 0 for c in CATEGORIES}
    if pd.isna(pathway_str):
        return out
    for path in str(pathway_str).split('|'):
        segments = [s.strip() for s in path.split('>')]
        for col, segment in CATEGORIES.items():
            if segment in segments:
                out[col] = 1
    return out


def main():
    mc = pd.read_excel(MITOCARTA_PATH, sheet_name='A Human MitoCarta3.0')
    print(f"MitoCarta3.0: {len(mc)} genes")

    flags = mc['MitoCarta3.0_MitoPathways'].apply(category_flags).apply(pd.Series)
    out = pd.concat([mc[['Symbol', 'MitoCarta3.0_MitoPathways']], flags], axis=1)
    out = out.rename(columns={'Symbol': 'Gene names'})
    out = out[['Gene names'] + list(CATEGORIES) + ['MitoCarta3.0_MitoPathways']]

    for col in CATEGORIES:
        print(f"  {col}: {int(out[col].sum())} genes")

    out.to_excel(OUT_PATH, index=False)
    print(f"Wrote {OUT_PATH}")

    # mitoproteome_isolated_mito.txt -- fig2B's reference background population. Per
    # the EV Methods ("Mitochondrial proteome data normalization and relative
    # enrichment"): "A curated, high-confidence mitochondrial proteome annotation
    # (MitoCarta3.0) was filtered to the subset of proteins detected in this
    # experiment to define a reference background population." That's exactly what
    # fig2B's own code does (intersects this file's gene set with its detected genes)
    # -- so this file is just the full MitoCarta3.0 gene symbol list, plain text,
    # one token per gene, matching the manuscript's stated method exactly.
    with open(MITOPROTEOME_OUT_PATH, 'w') as f:
        f.write('\n'.join(sorted(mc['Symbol'].dropna().astype(str))))
    print(f"Wrote {MITOPROTEOME_OUT_PATH} ({mc['Symbol'].notna().sum()} genes)")


if __name__ == '__main__':
    main()
