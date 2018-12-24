import gsd
import gsd.reactome
import gsd.immune_cells
import gsd.annotation

from pandas import read_table

human_tax_id = 9606
reactome_sub_tree_ids = ['R-HSA-8982491', 'R-HSA-1474290']

rule all:
    input:
        "evaluation_data/R-HSA-8982491",
        "evaluation_data/R-HSA-1474290",
        "evaluation_data/immune_cells"

rule download_reactome_sub_tree:
    output:
        [directory("evaluation_data/%s" % reactome_id)
        for reactome_id in reactome_sub_tree_ids]
    run:
        for reactome_id in reactome_sub_tree_ids:
            node, gene_sets = gsd.reactome.download(human_tax_id, reactome_id)
            gsd.persist_reference_data(node, gene_sets, "evaluation_data/%s" % reactome_id)

rule etract_immuno_cell_data:
    input:
        raw_data = directory("raw_data/immune_cells"),
        entrezgene2gene_sym = "annotation_data/entrezgene2gene_sym.tsv"
    output:
        dir = directory("evaluation_data/immune_cells")
    run:
        gene_sym_hsapiens = read_table(input.entrezgene2gene_sym)
        node, gene_sets = gsd.immune_cells.extract_from_raw_data(input.raw_data, gene_sym_hsapiens)
        gsd.persist_reference_data(node, gene_sets, output.dir)

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
            [gsd.annotation.BIOMART_GO_ID,
             gsd.annotation.BIOMART_GO_NAME,
             gsd.annotation.BIOMART_GO_DEFINITION,
             gsd.annotation.BIOMART_GO_LINKAGE_TYPE,
             gsd.annotation.BIOMART_GO_NAMESPACE],
            output.anno_file)