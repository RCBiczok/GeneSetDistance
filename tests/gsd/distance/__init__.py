import gsd.gene_sets

gene_sets = [
    gsd.gene_sets.GeneSet(name="SetA",
                          external_id="SetA",
                          external_source="",
                          summary="",
                          calculated=True,
                          entrez_gene_ids={8908, 2998, 2997},
                          gene_symbols={"GYG2", "GYS2", "GYS1"}),
    gsd.gene_sets.GeneSet(name="SetB",
                          external_id="SetB",
                          external_source="",
                          summary="",
                          calculated=True,
                          entrez_gene_ids={5507, 8908, 2998},
                          gene_symbols={"PPP1R3C", "GYG2", "GYS2"}),
    gsd.gene_sets.GeneSet(name="SetC",
                          external_id="SetC",
                          external_source="",
                          summary="",
                          calculated=True,
                          entrez_gene_ids={2998, 2997, 2992},
                          gene_symbols={"GYS2", "GYS1", "GYG1"})
]