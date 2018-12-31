import nltk
from pandas import read_table
from pathlib import Path

import gsd
import gsd.distance
import gsd.reactome
import gsd.immune_cells
import gsd.annotation
import gsd.gene_sets

from gsd.distance.general import GENERAL_DISTS, GENERAL_DISTS_TITLES
from gsd.distance.nlp import NLP_DISTS, NLP_DISTS_TITLES
#from gsd.distance.go import GO_DISTS, GO_DISTS_TITLES
from gsd.distance.ppi import PPI_DISTS, PPI_DISTS_TITLES

## General Variables

HUMAN_TAX_ID = 9606
STOPWORD_FILE = "%s/nltk_data/corpora/stopwords" % str(Path.home())

## Variables for evaluation data

REACTOME_TARGETS = ['R-HSA-8982491', 'R-HSA-1474290']
#EVALUATION_TARGETS = REACTOME_TARGETS + ['immune_cells']
EVALUATION_TARGETS = REACTOME_TARGETS
EVALUATION_DATA_DIRS = ["evaluation_data/%s" % target_name for target_name in EVALUATION_TARGETS]

## Used distances

COSINE_DISTANCE_OVER_GENE_SYM = "CosineDistanceOverGeneSym"

#EVALUATION_METRICS = GENERAL_DISTS_TITLES + NLP_DISTS_TITLES + GO_DISTS_TITLES + PPI_DISTS_TITLES
EVALUATION_METRICS = GENERAL_DISTS_TITLES + NLP_DISTS_TITLES + PPI_DISTS_TITLES

TARGET_OUTPUT = expand("experiment_data/{metrics}/{evaluation_target}.json",
                       metrics=EVALUATION_METRICS,
                       evaluation_target=EVALUATION_TARGETS)

rule all:
    input:
        TARGET_OUTPUT,

###
# Distance Measure Experiment
###

rule calc_general_dists:
    input: EVALUATION_DATA_DIRS
    output:
        expand("experiment_data/{metrics}/{evaluation_target}.json",
               metrics=GENERAL_DISTS_TITLES,
               evaluation_target=EVALUATION_TARGETS)
    run:
        for dist_info in GENERAL_DISTS:
            for evaluation_target in EVALUATION_TARGETS:
                gene_sets = gsd.gene_sets.load("evaluation_data/%s/gene_sets.json" % evaluation_target)
                gsd.distance.execute_and_persist_evaluation(
                        dist_info['distance'],
                        gene_sets,
                        evaluation_target,
                        "experiment_data/%s" % dist_info['folder'])

#rule calc_go_dists:
#    input: EVALUATION_DATA_DIRS
#    output:
#        expand("experiment_data/{metrics}/{evaluation_target}.json",
#               metrics=GO_DISTS_TITLES,
#               evaluation_target=EVALUATION_TARGETS)
#    run:
#        for dist_info in GO_DISTS:
#            for evaluation_target in EVALUATION_TARGETS:
#                gene_sets = gsd.gene_sets.load("evaluation_data/%s/gene_sets.json" % evaluation_target)
#                gsd.distance.execute_and_persist_evaluation(
#                        dist_info['distance'],
#                        gene_sets,
#                        evaluation_target,
#                        "experiment_data/%s" % dist_info['folder'])

rule calc_nlp_dists:
    input: EVALUATION_DATA_DIRS, STOPWORD_FILE
    output:
        expand("experiment_data/{metrics}/{evaluation_target}.json",
               metrics=NLP_DISTS_TITLES,
               evaluation_target=EVALUATION_TARGETS)
    run:
        from gensim.models import KeyedVectors
        #TODO embeddings are loaded no
        w2v_model = KeyedVectors.load_word2vec_format("__data/nlp/PubMed-Wilbur-2018/pubmed_s100w10_min.bin",
                                                      binary=True)

        for dist_info in NLP_DISTS:
            for evaluation_target in EVALUATION_TARGETS:
                gene_sets = gsd.gene_sets.load("evaluation_data/%s/gene_sets.json" % evaluation_target)
                gsd.distance.execute_and_persist_evaluation(
                        dist_info['distance_factory'](w2v_model),
                        gene_sets,
                        evaluation_target,
                        "experiment_data/%s" % dist_info['folder'])

rule calc_ppi_dists:
    input: EVALUATION_DATA_DIRS
    output:
        expand("experiment_data/{metrics}/{evaluation_target}.json",
               metrics=PPI_DISTS_TITLES,
               evaluation_target=EVALUATION_TARGETS)
    run:
        from gsd.distance.ppi import load_ppi_mitab
        #TODO embeddings are loaded no
        ppi_data = load_ppi_mitab("__data/ppi/BioGrid/BIOGRID-ALL-3.5.166.mitab.txt", HUMAN_TAX_ID)


        for dist_info in PPI_DISTS:
            for evaluation_target in EVALUATION_TARGETS:
                gene_sets = gsd.gene_sets.load("evaluation_data/%s/gene_sets.json" % evaluation_target)
                gsd.distance.execute_and_persist_evaluation(
                        dist_info['distance_factory'](ppi_data),
                        gene_sets,
                        evaluation_target,
                        "experiment_data/%s" % dist_info['folder'])

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
        [directory("evaluation_data/%s" % reactome_id) for reactome_id in REACTOME_TARGETS]
    run:
        go_anno = gsd.annotation.read_go_anno_df(input.entrezgene2go, input.go)
        for reactome_id in REACTOME_TARGETS:
            node, gene_sets = gsd.reactome.download(HUMAN_TAX_ID, reactome_id, go_anno)
            gsd.persist_reference_data(node, gene_sets, "evaluation_data/%s" % reactome_id)

rule extract_immuno_cell_data:
    input:
        raw_data = directory("raw_data/immune_cells"),
        entrezgene2gene_sym = "annotation_data/entrezgene2gene_sym.tsv",
        entrezgene2go = 'annotation_data/entrezgene2go.tsv',
        go = 'annotation_data/go.tsv'
    output:
        dir = directory("evaluation_data/immune_cells")
    run:
        gene_sym_hsapiens = read_table(input.entrezgene2gene_sym)
        go_anno = gsd.annotation.read_go_anno_df(input.entrezgene2go, input.go)
        node, gene_sets = gsd.immune_cells.extract_from_raw_data(input.raw_data, gene_sym_hsapiens, go_anno)
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
