import gsd.annotation
from gsd.gene_sets import GeneSetInfo, annotate_with_go

go_anno = gsd.annotation.read_go_anno_df("../annotation_data/entrezgene2go.tsv", "../annotation_data/go.tsv")

gene_sets_info_list = [
    GeneSetInfo(name="SetA",
                external_id="SetA",
                external_source="",
                summary="Lysyl hydroxylase (LH) (E.C. 1.14.11.4) is a dimeric enzyme that catalyzes the "
                        "formation of (2S,5R)-5-hydroxylysyl residues (5-Hyl) in proteins (reviewed in "
                        "Myllyharju & Kivirikko 2001) within a peptide linkage at the Y position of the "
                        "repeating X-Y-Gly sequence motif. The extent of 5-Hyl formation is much more "
                        "variable than that of hydroxyproline. It varies between collagen types, tissues "
                        "and by physiological state (Miller 1984). 5-Hyl content also differs between the "
                        "helical and telopeptide domains. This and the observation that purified lysyl "
                        "hydroxylase failed to hydroxylate Lys in the telopeptide domains has led to "
                        "speculation that there are separate enzymes responsible for Lys hydroxylation in "
                        "the helical and telopeptide domains (Royce & Barnes 1985, Gerriets et al. 1993). "
                        "The LH2b isoform may be the telopeptide-specific form (Pornprasertsuk et al. 2004)."
                        "<br>In human type I collagen, there are 38 residues of Lys in each alpha-1 chain "
                        "(36 in the helical domain, 1 each in the C- and N-telopeptide domains) and 31 in "
                        "each alpha-2 chain (30 in the helical domain,1 in the N-telopeptide and none in the "
                        "C-telopeptide domains) (Yamauchi & Shiiba 2002).  <br><br>LH requires ferrous iron, "
                        "oxygen, 2-oxoglutarate, and ascorbate. The hydroxylation reaction occurs during "
                        "collagen biosynthesis in the ER as a co- and post-translational event, before "
                        "triple helix formation. Three LH isoforms have been characterized in humans, "
                        "encoded by the genes PLOD1-3. The isoforms appear to have preferences for certain "
                        "collagen types, e.g. LH3 preferentially binds collagen types IV, VI, XI and XII "
                        "(Myllyla et al. 2007). LH3 has galactosyltransferase and glucosyltransferase "
                        "activities in addition to its lysyl hydroxylase activity (Heikkinen et al. 2000, "
                        "Wang et al. 2002), a multifunctionality also seen in the single C. elegans "
                        "orthologue (Wang et al. 2002a, b). Hydroxylysine residues can form stable "
                        "intermolecular cross-links between collagen molecules in fibrils and also "
                        "represent sites for glucosyl- and galactosyl- carbohydrate attachment. "
                        "<br><br>In this reaction all collagen subtypes are represented as having a single "
                        "hydroxylysine.",
                calculated=True,
                entrez_gene_ids={8908, 2998, 2997},
                gene_symbols={"GYG2", "GYS2", "GYS1"}),
    GeneSetInfo(name="SetB",
                external_id="SetB",
                external_source="",
                summary="Alignment of the C-terminal prodomains initiates triple helix formation, which"
                        "propagates in a zipper-like fashion in the C-to-N direction. This occurs in the "
                        "rough endoplasmic reticulum (Engel & Prockop 1991). Compared with the folding of "
                        "globular proteins and coiled-coil structures, the concentration-independent folding "
                        "steps of collagen are extremely slow (Bächinger et al. 1980). Triple helix "
                        "formation combines a fast process, interpreted as the folding of regions devoid "
                        "of cis residues, and a slow process, limited by the slow kinetics of cis to trans "
                        "prolyl-isomerization (Bächinger et al. 1978). Triple-helix formation in regions "
                        "devoid of cis-prolyl bonds is 3-4 times faster than formation limited by prolyl "
                        "isomerization reactions  (Bachmann et al. 2005). Conversion to trans is required "
                        "as only trans-peptide bonds can be incorporated into the collagen triple helix "
                        "(Zeng et al. 1998). Efficient helix folding requires the presence of the 3-prolyl "
                        "hydroxylation complex. This trimer of Prolyl 3-hydroxylase 1 (LEPRE1), "
                        "Cyclophilin B (CyPB), also called Peptidyl-prolyl cis-trans isomerase B (PPIB) and "
                        "CRTAP has 3-prolyl hydroxylase, PPIase and procollagen chaperone properties "
                        "(Ishikawa et al. 2009, van Dijk et al. 2009). Efficient folding involves additional "
                        "collagen-specific chaperones such as Serpin H1 (HSP47 - Smith et al. 1995). "
                        "CyPB belongs to the cyclophilins, a conserved class of intracellular and/or "
                        "secreted proteins originally identified as cellular binding proteins for the "
                        "immunosuppressive drug cyclosporin A. They are peptidyl-prolyl cis-trans isomerases "
                        "(PPIases), which catalyze the cis-trans isomerisation of peptide bonds. CyPB "
                        "localises to the rough ER but is also secreted extracellularly. It directly "
                        "interacts with procollagen and is believed to  be responsible for converting "
                        "procollagen cis- to trans-conformers (Zeng et al. 1998). CyPB and Serpin H1 are "
                        "also involved in procollagen export and secretion. Results obtained with collagen "
                        "peptides suggest that variations in the Gly-X-Y sequence are likely to result in a "
                        "non-uniform helical twist along the length of a collagen fibril. Sequences poor in "
                        "imino acids will have a symmetry close to 10 tripeptide units for every 3 turns of "
                        "the triple helix (10/3), while stretches of Gly-Pro-Hyp units may have 7/2 symmetry "
                        "(Brodsky & Persikov 2005).",
                calculated=True,
                entrez_gene_ids={5507, 8908, 2998},
                gene_symbols={"PPP1R3C", "GYG2", "GYS2"}),
    GeneSetInfo(name="SetC",
                external_id="SetC",
                external_source="",
                summary="Collagen was for many years considered the only source of 4-hydroxyproline (4-Hyp) "
                        "in animals. Though it is now known that other proteins such as C1q and elastin also "
                        "contain 4-Hyp, collagen is by far the major source (Adams & Frank 1980). 4-Hyp is "
                        "required for collagen stability at physiological temperatures. The abundance of Hyp "
                        "in animal proteins is ~4%, making it more abundant than the amino-acids Cys, Gln, "
                        "His, Met, Phe, Trp and Tyr (McCaldon & Argos 1988). In collagen Hyp abundance is "
                        "much higher at ~38% (Ramshaw et al. 1998). Full collagen proline hydroxylation "
                        "significantly raises the melting temperature (Tm) by stabilizing the collagen "
                        "triple helix (Berg & Prockop 1973a), a process that has been studied extensively "
                        "using synthetic collagen peptides (Sakakibara et al. 1973,  Holmgren et al. 1998) "
                        "and is well understood at the structural level (Shoulders & Raines 2009). Collagen "
                        "4-Hyp content is relatively stable, with small differences between collagen types. "
                        "Collagen type I has approximately 1 4-Hyp for every 10 residues, roughly 50% of "
                        "available prolines (Kivirikko et al. 1992). Conversion of Pro to "
                        "(2S,4R)-4-hydroxyproline (4-Hyp) is the most prevalent posttranslational "
                        "modification in humans, catalyzed by prolyl 4-hydroxylase (P4H).  Mammalian prolyl "
                        "4-hydroxylase is an alpha2 beta2 tetramer (Berg & Prockop 1973b). The 59-kDa alpha "
                        "subunit contains the substrate-binding domain and the enzymic active site "
                        "(Helaakoski et al. 1989). Humans and most other vertebrates have three isoforms of "
                        "the alpha subunit, isoform alpha-1 is the most prevalent. The pair of alpha "
                        "subunits can be any of the three isoforms (Gorres & Raines 2010). The 55-kDa beta "
                        "subunit is protein disulphide isomerase (PDI), which has additional functions in "
                        "collagen formation. As part of P4H it retains the tetramer in the ER lumen and "
                        "maintains the otherwise insoluble alpha subunit in an active form (Vuori et al. "
                        "1992, Nietfeld & Kemp 1981). P4H is a member of the non-heme iron(II), "
                        "alpha-ketoglutarate-dependent dioxygenase family. Molecular oxygen (O2), "
                        "2-oxoglutarate (alpha-ketoglutarate) and iron(II) are required for its activity "
                        "(Hutton & Udenfriend 1966). During the reaction, alpha-ketoglutarate is oxidatively "
                        "decarboxylated producing succinate and CO2 (Rhoads & Udenfriend 1968, Gorres & "
                        "Raines 2010). Ascorbate is required as a cofactor but not consumed "
                        "(Kivirikko et al 1989). The minimum substrate required for hydroxylation is an "
                        "Xaa-Pro-Gly tripeptide, with Pro preferred in the Xaa position, though "
                        "hydroxylation can occur at lower rates with a variety of residues at this "
                        "position (Kivirikko et al. 1972). A number of other peptides, notably elastin, "
                        "are substrates for P4H (Bhatnagar 1978).<br><br>For brevity, all forms of collagen "
                        "propeptide are shown as having 3X 4-Hyp residues following the action of P4H.",
                calculated=True,
                entrez_gene_ids={2998, 2997, 2992},
                gene_symbols={"GYS2", "GYS1", "GYG1"})
]

gene_sets = annotate_with_go(gene_sets_info_list, go_anno)
