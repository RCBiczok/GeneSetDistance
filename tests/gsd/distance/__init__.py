import gsd.gene_sets

go_anno = gsd.gene_sets.read_go_anno_df("../annotation_data/entrezgene2go.tsv", "../annotation_data/go.tsv")


gene_sets = gsd.gene_sets.load_gene_sets(gene_sets_file="gsd/distance/fake_gene_sets.json",
                                         ncbi_gene_desc_file="gsd/distance/fake_ncbi_gene_desc.json",
                                         gwas_gene_traits_file="gsd/distance/fake_gwas_traits.json")
