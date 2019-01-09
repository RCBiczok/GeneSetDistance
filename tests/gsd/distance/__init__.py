import gsd.gene_sets
import gsd.annotation

go_anno = gsd.annotation.read_go_anno_df("../annotation_data/entrezgene2go.tsv", "../annotation_data/go.tsv")


gene_sets = gsd.gene_sets.load_gene_sets("gsd/distance/fake_gene_sets.json",
                                         "gsd/distance/fake_ncbi_gene_desc.json")
