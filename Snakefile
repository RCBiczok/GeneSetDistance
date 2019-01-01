import nltk
from pandas import read_table
from pathlib import Path

import gsd
import gsd.annotation
import gsd.distance
import gsd.reactome
import gsd.immune_cells
import gsd.annotation
import gsd.gene_sets

from gsd.distance.general import GENERAL_DISTS
from gsd.distance.nlp import NLP_DISTS
from gsd.distance.ppi import PPI_DISTS

## General Variables

HUMAN_TAX_ID = 9606
STOPWORD_FILE = "%s/nltk_data/corpora/stopwords" % str(Path.home())

## Variables for evaluation data

REACTOME_TARGETS = ['reactome/R-HSA-8982491', 'reactome/R-HSA-1474290']
#EVALUATION_TARGETS = REACTOME_TARGETS + ['immune_cells/all']
EVALUATION_TARGETS = REACTOME_TARGETS
EVALUATION_DATA_DIRS = ["evaluation_data/%s" % target_name for target_name in EVALUATION_TARGETS]

## Used distances

GENERAL_EVALUATION_OUTPUT = expand("experiment_data/general/{metric}/{evaluation_target}.json",
                                   metric=GENERAL_DISTS.keys(),
                                   evaluation_target=EVALUATION_TARGETS)

NLP_EVALUATION_OUTPUT = expand("experiment_data/nlp/{metric}/{evaluation_target}.json",
                               metric=NLP_DISTS.keys(),
                               evaluation_target=EVALUATION_TARGETS)

PPI_EVALUATION_OUTPUT = expand("experiment_data/ppi/{metric}/{evaluation_target}.json",
                               metric=PPI_DISTS.keys(),
                               evaluation_target=EVALUATION_TARGETS)

GO_DISTS = {
    'GO_SIM_BP_Wang_BMA': {'type': gsd.annotation.GOType.BIOLOGICAL_PROCESS, 'measure': "Wang", 'combine': "BMA"}
}

# {
#     'folder': "GO_SIM_BP_Resnik_BMA",
#     'distance': GOSimDistanceMetric(GOType.BIOLOGICAL_PROCESS, "Resnik", "BMA")
# }

GO_EVALUATION_OUTPUT = expand("experiment_data/go/{metric}/{evaluation_target}.json",
                               metric=GO_DISTS.keys(),
                               evaluation_target=EVALUATION_TARGETS)


###
# Default rule
###

rule all:
    input:
        GENERAL_EVALUATION_OUTPUT,
        NLP_EVALUATION_OUTPUT,
        PPI_EVALUATION_OUTPUT,
        GO_EVALUATION_OUTPUT

###
# Distance Measure Experiment
###

rule calc_general_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: file="experiment_data/general/{metric}/{target_category}/{evaluation_target}.json"
    run:
        dist = GENERAL_DISTS[wildcards.metric]
        gene_sets = gsd.gene_sets.load_gene_sets(input.file )
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)


rule calc_nlp_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json",
           stopwords_file=STOPWORD_FILE
    output: file="experiment_data/nlp/{metric}/{target_category}/{evaluation_target}.json"
    run:
        from gensim.models import KeyedVectors

        #TODO embeddings are loaded no
        print("loading w2v model")
        w2v_model = KeyedVectors.load_word2vec_format("__data/nlp/PubMed-Wilbur-2018/pubmed_s100w10_min.bin",
                                                      binary=True)

        dist = NLP_DISTS[wildcards.metric](w2v_model)
        print("Perform calculation for: %s" % dist.display_name)

        gene_sets = gsd.gene_sets.load_gene_sets(input.file)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)

rule calc_ppi_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: file="experiment_data/ppi/{metric}/{target_category}/{evaluation_target}.json"
    run:
        from gsd.distance.ppi import load_ppi_mitab

        #TODO embeddings are loaded no
        print("loading PPI data")
        ppi_data = load_ppi_mitab("__data/ppi/BioGrid/BIOGRID-ALL-3.5.166.mitab.txt", HUMAN_TAX_ID)

        dist = PPI_DISTS[wildcards.metric](ppi_data)
        print("Perform calculation for: %s" % dist.display_name)

        gene_sets = gsd.gene_sets.load_gene_sets(input.file)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)


rule calc_go_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: file="experiment_data/go/{metric}/{target_category}/{evaluation_target}.json"
    run:
        from gsd.distance.go import GOSimDistanceMetric

        dist_info = GO_DISTS[wildcards.metric]
        dist =  GOSimDistanceMetric(dist_info['type'], dist_info['measure'], dist_info['combine'])
        gene_sets = gsd.gene_sets.load_gene_sets(input.file)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)

###
# Data Download & initialization
###

rule download_stopwords:
    output:
        directory(STOPWORD_FILE)
    run:
        nltk.download("stopwords")

rule download_reactome_sub_tree:
    input:
        entrezgene2go = 'annotation_data/entrezgene2go.tsv',
        go = 'annotation_data/go.tsv'
    output:
        gene_set_file = "evaluation_data/reactome/{evaluation_target}/gene_sets.json",
        tree_file = "evaluation_data/reactome/{evaluation_target}/tree.json"
    run:
        go_anno = gsd.annotation.read_go_anno_df(input.entrezgene2go, input.go)
        node, gene_sets = gsd.reactome.download(HUMAN_TAX_ID, wildcards.evaluation_target, go_anno)
        gsd.persist_reference_data(node, gene_sets, output.tree_file, output.gene_set_file)

rule extract_all_immuno_cell_data:
    input:
        raw_data = directory("raw_data/immune_cells"),
        entrezgene2gene_sym = "annotation_data/entrezgene2gene_sym.tsv",
        entrezgene2go = 'annotation_data/entrezgene2go.tsv',
        go = 'annotation_data/go.tsv'
    output:
        gene_set_file = "evaluation_data/immune_cells/all/gene_sets.json",
        tree_file = "evaluation_data/immune_cells/all/tree.json"
    run:
        gene_sym_hsapiens = read_table(input.entrezgene2gene_sym)
        go_anno = gsd.annotation.read_go_anno_df(input.entrezgene2go, input.go)
        node, gene_sets = gsd.immune_cells.extract_from_raw_data(input.raw_data, gene_sym_hsapiens, go_anno)
        gsd.persist_reference_data(node, gene_sets, output.tree_file, output.gene_set_file)

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
