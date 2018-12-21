import gsd.reactome
import gsd.annotation

human_tax_id = 9606
reactome_sub_tree_ids = ['R-HSA-8982491', 'R-HSA-1474290']

rule all:
    input:
        "evaluation_data/R-HSA-8982491",
        "evaluation_data/R-HSA-1474290",
        "annotation_data/entrezgene2gene_sym.tsv",
        "annotation_data/entrezgene2go.tsv",
        "annotation_data/go.tsv"

rule download_reactome_sub_tree:
    output:
        [directory("evaluation_data/%s" % reactome_id)
        for reactome_id in reactome_sub_tree_ids]
    run:
        for reactome_id in reactome_sub_tree_ids:
            gsd.reactome.download_reactome_sub_tree(human_tax_id, reactome_id,
                                                    "evaluation_data/%s" % reactome_id)

rule download_entrezgene2gene_sym_anno:
    output:
        anno_file = "annotation_data/entrezgene2gene_sym.tsv"
    run:
        gsd.annotation.download_biomart_anno(
            ["external_gene_name", "entrezgene"],
            output.anno_file)

rule download_entrezgene2go_anno:
    output:
        anno_file = "annotation_data/entrezgene2go.tsv"
    run:
        gsd.annotation.download_biomart_anno(
            ['entrezgene', 'go_id'],
            output.anno_file)

rule download_go_anno:
    output:
        anno_file = "annotation_data/go.tsv"
    run:
        gsd.annotation.download_biomart_anno(
            ['go_id', 'name_1006', 'definition_1006', 'go_linkage_type', 'namespace_1003'],
            output.anno_file)