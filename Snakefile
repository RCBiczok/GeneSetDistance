import nltk
from pandas import read_table
from pathlib import Path
from anytree import PostOrderIter

import gsd
import gsd.annotation
import gsd.distance
import gsd.reactome
import gsd.immune_cells
import gsd.annotation
import gsd.gene_sets

from gsd.distance.general import GENERAL_DISTS
from gsd.distance.benchmark import BENCHMARK_DISTS
from gsd.distance.nlp import NLP_DISTS
from gsd.distance.ppi import PPI_DISTS

## General Variables

HUMAN_TAX_ID = 9606
STOPWORD_FILE = "%s/nltk_data/corpora/stopwords" % str(Path.home())

## Variables for evaluation data

REACTOME_TARGETS = ['reactome/R-HSA-8982491',
                    'reactome/R-HSA-1474290',
                    'reactome/R-HSA-373755',
                    'reactome/R-HSA-422475']
#EVALUATION_TARGETS = REACTOME_TARGETS + ['immune_cells/all']
EVALUATION_TARGETS = REACTOME_TARGETS + ['immune_cells/immune_only']
EVALUATION_DATA_DIRS = ["evaluation_data/%s" % target_name for target_name in EVALUATION_TARGETS]

## Used distances

GENERAL_EVALUATION_OUTPUT = expand("experiment_data/general/{metric}/{evaluation_target}.json",
                                   metric=GENERAL_DISTS.keys(),
                                   evaluation_target=EVALUATION_TARGETS)

BENCHMARK_EVALUATION_OUTPUT = expand("experiment_data/benchmark/{metric}/{evaluation_target}.json",
                                      metric=BENCHMARK_DISTS.keys(),
                                      evaluation_target=EVALUATION_TARGETS)

NLP_EVALUATION_OUTPUT = expand("experiment_data/nlp/{metric}/{evaluation_target}.json",
                               metric=NLP_DISTS.keys(),
                               evaluation_target=EVALUATION_TARGETS)

PPI_EVALUATION_OUTPUT = expand("experiment_data/ppi/{metric}/{evaluation_target}.json",
                               metric=PPI_DISTS.keys(),
                               evaluation_target=EVALUATION_TARGETS)

GO_DISTS = {
    'GO_SIM_BP_Wang_BMA': {'type': gsd.annotation.GOType.BIOLOGICAL_PROCESS, 'measure': "Wang", 'combine': "BMA"},
    'GO_SIM_CC_Wang_BMA': {'type': gsd.annotation.GOType.CELLULAR_COMPONENT, 'measure': "Wang", 'combine': "BMA"},
    'GO_SIM_MF_Wang_BMA': {'type': gsd.annotation.GOType.MOLECULAR_FUNCTION, 'measure': "Wang", 'combine': "BMA"}
}

GO_EVALUATION_OUTPUT = expand("experiment_data/go/{metric}/{evaluation_target}.json",
                               metric=GO_DISTS.keys(),
                               evaluation_target=EVALUATION_TARGETS)

TREE_PATH_OUTPUT = expand("experiment_data/tree_path/{evaluation_target}.json",
                          evaluation_target=EVALUATION_TARGETS)

###
# Default rule
###

rule all:
    input:
        GENERAL_EVALUATION_OUTPUT,
        BENCHMARK_EVALUATION_OUTPUT,
        NLP_EVALUATION_OUTPUT,
        PPI_EVALUATION_OUTPUT,
        GO_EVALUATION_OUTPUT,
        TREE_PATH_OUTPUT

###
# Distance Measure Experiment
###

rule calc_general_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json",
           gwas_gene_traits_file="evaluation_data/{target_category}/{evaluation_target}/gwas_gene_traits.json"
    output: file="experiment_data/general/{metric}/{target_category}/{evaluation_target}.json"
    run:
        dist = GENERAL_DISTS[wildcards.metric]
        gene_sets = gsd.gene_sets.load_gene_sets(input.file,
                                                 gwas_gene_traits_file=input.gwas_gene_traits_file)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)


rule calc_benchmark_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: file="experiment_data/benchmark/{metric}/{target_category}/{evaluation_target}.json"
    run:
        dist = BENCHMARK_DISTS[wildcards.metric]
        gene_sets = gsd.gene_sets.load_gene_sets(input.file)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)


rule calc_nlp_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json",
           ncbi_gene_desc_file="evaluation_data/{target_category}/{evaluation_target}/ncbi_gene_desc.json",
           stopwords_file=STOPWORD_FILE
    output: file="experiment_data/nlp/{metric}/{target_category}/{evaluation_target}.json"
    run:
        from gensim.models import KeyedVectors

        #TODO embeddings are not downloaded automatically
        print("Loading w2v model")
        w2v_model = KeyedVectors.load_word2vec_format("__data/nlp/PubMed-Wilbur-2018/pubmed_s100w10_min.bin",
                                                      binary=True)

        dist = NLP_DISTS[wildcards.metric](w2v_model)
        print("Perform calculation for: %s / %s" % (dist.display_name, wildcards.evaluation_target))

        gene_sets = gsd.gene_sets.load_gene_sets(input.file, input.ncbi_gene_desc_file)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)

rule calc_ppi_dists:
    input: file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: file="experiment_data/ppi/{metric}/{target_category}/{evaluation_target}.json"
    run:
        from gsd.distance.ppi import load_ppi_mitab

        #TODO ppi file should be downloaded automatically
        print("Loading PPI data")
        ppi_data = load_ppi_mitab("__data/ppi/BioGrid/BIOGRID-ALL-3.5.166.mitab.txt", HUMAN_TAX_ID)

        dist = PPI_DISTS[wildcards.metric](ppi_data)
        print("Perform calculation for: %s / %s" % (dist.display_name, wildcards.evaluation_target))

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


rule calc_tree_path_dists:
    input:
        gene_sets_file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json",
        tree_file="evaluation_data/{target_category}/{evaluation_target}/tree.json"
    output:
        file="experiment_data/tree_path/{target_category}/{evaluation_target}.json"
    run:
        root = gsd.gene_sets.load_tree(input.tree_file)
        gene_sets = gsd.gene_sets.load_gene_sets(input.gene_sets_file)
        dist =  gsd.distance.PairwiseTreePathDistanceMetric(root)
        gsd.distance.execute_and_persist_evaluation(dist, gene_sets, output.file)

###
# Data download & Data preparation
###

rule download_stopwords:
    output:
        directory(STOPWORD_FILE)
    run:
        nltk.download("stopwords")

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

rule extract_immune_cells_only:
    input:
        gene_set_file = "evaluation_data/immune_cells/all/gene_sets.json",
        tree_file = "evaluation_data/immune_cells/all/tree.json"
    output:
        gene_set_file = "evaluation_data/immune_cells/immune_only/gene_sets.json",
        tree_file = "evaluation_data/immune_cells/immune_only/tree.json"
    run:
        gene_sets = gsd.gene_sets.load_gene_sets(input.gene_set_file)
        root = gsd.gene_sets.load_tree(input.tree_file)

        filtered_root = root.children[0]
        filtered_names = [node.name for node in PostOrderIter(filtered_root)]
        filtered_gene_sets = [gene_set for gene_set in gene_sets if gene_set.general_info.name in filtered_names]

        gsd.persist_reference_data(filtered_root, filtered_gene_sets, output.tree_file, output.gene_set_file)

rule downlaod_ncbi_gene_desc:
    input: gene_set_file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: ncbi_gene_desc_file="evaluation_data/{target_category}/{evaluation_target}/ncbi_gene_desc.json"
    run:
        gene_sets = gsd.gene_sets.load_gene_sets(input.gene_set_file)
        gsd.annotation.downlaod_ncbi_gene_desc(gene_sets, output.ncbi_gene_desc_file)

rule extract_gwas_gene_traits:
    input: gene_set_file="evaluation_data/{target_category}/{evaluation_target}/gene_sets.json"
    output: gwas_gene_traits_file="evaluation_data/{target_category}/{evaluation_target}/gwas_gene_traits.json"
    run:
        #TODO mappings are not downloaded automatically
        gwas_gene_traigs = "__data/gwas/gwas_catalog_v1.0-associations_e93_r2018-12-21.tsv"

        gene_sets = gsd.gene_sets.load_gene_sets(input.gene_set_file)
        gsd.annotation.extract_gwas_traits(gwas_gene_traigs, gene_sets, output.gwas_gene_traits_file)