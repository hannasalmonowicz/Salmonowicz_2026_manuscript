#!/usr/bin/env python3
"""
Builds the two non-EV-table inputs fig2B_isolated_mito_pathway_ANOVA.py needs,
from a manually curated reference file,
'Figure_1_protein_lists_ISOLATED_MITOCHONDRIA.xlsx' -- the isolated-mito
analog of 'Figure_1_protein_lists_Whole_Cell.xlsx' already used by the
trajectory script.

Supersedes verification/build_MitoCarta_annotation_isolated_mito.py: that
script derived the same two files automatically from the raw official
MitoCarta3.0 database, which was a reasonable approximation (cross-checked
to reproduce this project's own established categorization logic with zero
mismatches) but not dataset-specific. This curated file IS dataset-specific
-- it already restricts every category to genes actually detected in the
isolated-mito experiment (EV_Table_2), separates structural subunits from
assembly factors before recombining them, and documents one correction
(Sulfur_metabolism, n=4 -> n=6) and one exclusion (PNKD, a spurious Complex
IV pathway-string match, flagged for removal in its own Flagged_review
sheet) that the automatic derivation could not know about. Cross-checked:
this file's SULFUR (6 genes), NEAA (16), 1C-core (4) and Nucleotide (25)
category definitions are byte-for-byte identical, gene for gene, to the
sets already hardcoded inside fig2B_isolated_mito_pathway_ANOVA.py itself --
i.e. this file and this script already agree; this build script just
supplies the remaining pieces (Complex I-V, Mito gene expression, Fatty
acid oxidation/BCAA tagging, and the mitoproteome reference) the same way.

Outputs (both in figures/data/):
  1. MitoCarta_pathway_annotation.xlsx -- 'Gene names' + ComplexI..ComplexV +
     Translation/mtDNA_maintenance/mtRNA metabolism/Lipid_metabolism (boolean)
     + MitoCarta3.0_MitoPathways (kept only as a carrier for the literal
     substring 'Branched-chain', which is all fig2B's own code actually reads
     from that column -- see its bcaa_m computation).
  2. mitoproteome_isolated_mito.txt -- the reference background gene list,
     now 'Verified_mitoproteome_ISO' (735 genes) instead of the full
     1136-gene MitoCarta3.0 database.

fig2B's Translation / mtDNA_maintenance / mtRNA metabolism columns are only
ever consumed as a union (mito_ge_m = Translation | mtDNA_maintenance |
mtRNA_metabolism) -- so all of Mito_gene_expression_ISO's 191 genes are
written to the Translation column alone (the other two left 0); the union
is identical either way.
"""
import pandas as pd

CURATED_PATH = 'Figure_1_protein_lists_ISOLATED_MITOCHONDRIA.xlsx'  # point at your copy
ANNOT_OUT_PATH = '../figures/data/MitoCarta_pathway_annotation.xlsx'
MITOPROTEOME_OUT_PATH = '../figures/data/mitoproteome_isolated_mito.txt'

COMPLEX_SHEETS = {
    'ComplexI': ('CI_subunits_ISO', 'CI_assembly_ISO'),
    'ComplexII': ('CII_subunits_ISO', 'CII_assembly_ISO'),
    'ComplexIII': ('CIII_subunits_ISO', 'CIII_assembly_ISO'),
    'ComplexIV': ('CIV_subunits_ISO', 'CIV_assembly_ISO'),
    'ComplexV': ('CV_subunits_ISO', 'CV_assembly_ISO'),
}


def genes(xl, sheet):
    return set(xl.parse(sheet)['Gene'].dropna().astype(str))


def main():
    xl = pd.ExcelFile(CURATED_PATH)

    flagged = xl.parse('Flagged_review')
    exclude_by_complex = {}
    for _, row in flagged.iterrows():
        exclude_by_complex.setdefault('Complex' + row['Complex'].replace('C', '', 1), set()).add(row['Gene'])
        print(f"Excluding {row['Gene']} from {row['Complex']}: {row['Reason']}")

    complex_genes = {}
    for col, (sub_sheet, asm_sheet) in COMPLEX_SHEETS.items():
        g = genes(xl, sub_sheet) | genes(xl, asm_sheet)
        g -= exclude_by_complex.get(col, set())
        complex_genes[col] = g
        print(f"  {col}: {len(g)} genes ({sub_sheet} + {asm_sheet}, minus flagged)")

    mito_ge_genes = genes(xl, 'Mito_gene_expression_ISO')
    fao_genes = genes(xl, 'Fatty_acid_oxidation_ISO')
    bcaa_genes = genes(xl, 'BCAA_metabolism_ISO')
    print(f"  Mito_gene_expression_ISO: {len(mito_ge_genes)} genes")
    print(f"  Fatty_acid_oxidation_ISO: {len(fao_genes)} genes")
    print(f"  BCAA_metabolism_ISO: {len(bcaa_genes)} genes")

    all_genes = set().union(*complex_genes.values(), mito_ge_genes, fao_genes, bcaa_genes)
    rows = []
    for g in sorted(all_genes):
        rows.append({
            'Gene names': g,
            'ComplexI': int(g in complex_genes['ComplexI']),
            'ComplexII': int(g in complex_genes['ComplexII']),
            'ComplexIII': int(g in complex_genes['ComplexIII']),
            'ComplexIV': int(g in complex_genes['ComplexIV']),
            'ComplexV': int(g in complex_genes['ComplexV']),
            'Translation': int(g in mito_ge_genes),
            'mtDNA_maintenance': 0,
            'mtRNA metabolism': 0,
            'Lipid_metabolism': int(g in fao_genes),
            'MitoCarta3.0_MitoPathways': 'Metabolism > Branched-chain' if g in bcaa_genes else '',
        })
    annot_df = pd.DataFrame(rows)
    annot_df.to_excel(ANNOT_OUT_PATH, index=False)
    print(f"Wrote {ANNOT_OUT_PATH} ({len(annot_df)} genes)")

    mitoproteome = genes(xl, 'Verified_mitoproteome_ISO')
    with open(MITOPROTEOME_OUT_PATH, 'w') as f:
        f.write('\n'.join(sorted(mitoproteome)))
    print(f"Wrote {MITOPROTEOME_OUT_PATH} ({len(mitoproteome)} genes)")


if __name__ == '__main__':
    main()
