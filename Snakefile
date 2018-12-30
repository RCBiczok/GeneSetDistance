import nltk
from pandas import read_table
from pathlib import Path

import gsd
import gsd.distance
import gsd.reactome
import gsd.immune_cells
import gsd.annotation
import gsd.gene_sets
from gsd.distance.general import MinkowskiNormDistanceMetric, JaccardDistanceMetric, KappaDistanceMetric

## General Variables

HUMAN_TAX_ID = 9606
STOPWORD_FILE = "%s/nltk_data/corpora/stopwords" % str(Path.home())

## Variables for annotations

ENTREZGENE_2_GO_FILE = 'annotation_data/entrezgene2go.tsv'
GO_FILE = 'annotation_data/go.tsv'

ANNOTATION_FILES = [ENTREZGENE_2_GO_FILE,
                    GO_FILE,
                    'annotation_data/entrezgene2gene_sym.tsv']

## Variables for evaluation data

REACTOME_TARGETS = ['R-HSA-8982491', 'R-HSA-1474290']
EVALUATION_TARGETS = REACTOME_TARGETS + ['immune_cells']
EVALUATION_DATA_DIRS = ["evaluation_data/%s" % target_name for target_name in EVALUATION_TARGETS]

EVALUATION_INPUT_DATA = ANNOTATION_FILES + EVALUATION_DATA_DIRS

## Used distances

DIST_MINKOWSKI_P2 = "Minkowski_P2"

EVALUATION_METRICS = [DIST_MINKOWSKI_P2]

TARGET_OUTPUT = expand("experiment_data/{metrics}/{evaluation_target}.json",
                       metrics=EVALUATION_METRICS,
                       evaluation_target=EVALUATION_TARGETS)

rule all:
    input:
        TARGET_OUTPUT,
        STOPWORD_FILE


###
# Distance Measure Experiment
###

rule calc_minkowski_p2:
    input: EVALUATION_INPUT_DATA
    output:
        expand("experiment_data/{metrics}/{evaluation_target}.json",
               metrics=[DIST_MINKOWSKI_P2],
               evaluation_target=EVALUATION_TARGETS)
    run:
        go_anno = gsd.annotation.read_go_anno_df(ENTREZGENE_2_GO_FILE, GO_FILE)

        dist = MinkowskiNormDistanceMetric()

        for evaluation_target in EVALUATION_TARGETS:
            gene_sets = gsd.gene_sets.load("evaluation_data/%s/gene_sets.json" % evaluation_target, go_anno)
            gsd.distance.execute_and_persist_evaluation(
                        dist,
                        gene_sets,
                        evaluation_target,
                        "experiment_data/%s" % DIST_MINKOWSKI_P2)

###
# Data Download & initialization
###

rule download_stopwords:
    output:
        directory(STOPWORD_FILE)
    run:
        nltk.download("stopwords")

rule download_reactome_sub_tree:
    output:
        [directory("evaluation_data/%s" % reactome_id) for reactome_id in REACTOME_TARGETS]
    run:
        for reactome_id in REACTOME_TARGETS:
            node, gene_sets = gsd.reactome.download(HUMAN_TAX_ID, reactome_id)
            gsd.persist_reference_data(node, gene_sets, "evaluation_data/%s" % reactome_id)

rule extract_immuno_cell_data:
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
