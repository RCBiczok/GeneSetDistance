# ATTENTION: if tests on a new machine fail. run snakemake first to download nltk files first

from gensim.models import KeyedVectors

from gsd.distance.nlp import NLPCosineDistance, extract_gene_symbols, \
    extract_words_from_gene_set_summary
from tests import has_equal_elements
from tests.gsd.distance import gene_sets

w2v_model = KeyedVectors.load_word2vec_format("../__data/nlp/PubMed-Wilbur-2018/pubmed_s100w10_min.bin", binary=True)


def test_extract_gene_symbols():
    gene_set = gene_sets[0]
    assert has_equal_elements(set(extract_gene_symbols(gene_set, w2v_model)),
                              {'gys1', 'gys2', 'gyg2'})


def test_extract_words_from_gene_set_summary():
    summary_words = extract_words_from_gene_set_summary(gene_sets[0], w2v_model)
    assert has_equal_elements(summary_words[:3],
                              ['lysyl', 'hydroxylase', 'lh'])


def test_cosine_distance_over_gene_symbols():
    dist_metric = NLPCosineDistance(w2v_model)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.080, 0.015, 0.101], epsilon=0.001)
