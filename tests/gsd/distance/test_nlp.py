# ATTENTION: if tests on a new machine fail. run snakemake first to download nltk files first

from gensim.models import KeyedVectors

from gsd.annotation import GOType
from gsd.distance.nlp import NLPDistance, extract_gene_symbols, \
    extract_words_from_gene_set_summary, extract_words_from_go_descriptions, cosine_distance_of, \
    extract_words_from_go_names, extract_words_from_gene_symbols_and_summary_and_go_info, wm_distance_of
from tests import has_equal_elements
from tests.gsd.distance import gene_sets

w2v_model = KeyedVectors.load_word2vec_format("../__data/nlp/PubMed-Wilbur-2018/pubmed_s100w10_min.bin", binary=True)


def test_extract_gene_symbols():
    gene_set = gene_sets[0]
    assert set(extract_gene_symbols(gene_set, w2v_model)) == {'gys1', 'gys2', 'gyg2'}


def test_extract_words_from_gene_set_summary():
    summary_words = extract_words_from_gene_set_summary(gene_sets[0], w2v_model)
    assert has_equal_elements(summary_words[:3],
                              ['lysyl', 'hydroxylase', 'lh'])


def test_extract_words_from_go_name():
    words = extract_words_from_go_names(gene_sets[0],
                                        w2v_model,
                                        [GOType.BIOLOGICAL_PROCESS])
    assert {'glycogen', 'biosynthetic', 'process'}.issubset(words)


def test_extract_words_from_go_descriptions():
    words = extract_words_from_go_descriptions(gene_sets[0],
                                               w2v_model,
                                               [GOType.BIOLOGICAL_PROCESS])
    assert {'process', 'results', 'change'}.issubset(words)


def test_extract_words_from_gene_symbols_and_summary_and_go_info():
    words = extract_words_from_gene_symbols_and_summary_and_go_info(gene_sets[0],
                                                                    w2v_model)
    assert {'process', 'results', 'change', 'glycogen', 'biosynthetic',
            'process', 'lysyl', 'hydroxylase', 'lh'}.issubset(words)


def test_cosine_distance_over_gene_symbols():
    dist_metric = NLPDistance("Cosine Distance over gene symbols",
                              lambda x, y: cosine_distance_of(x, y, w2v_model),
                              lambda x: extract_gene_symbols(x, w2v_model))
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.080, 0.015, 0.101], epsilon=0.001)


def test_wm_distance_over_gene_symbols():
    dist_metric = NLPDistance("WM Distance over gene symbols",
                              lambda x, y: wm_distance_of(x, y, w2v_model),
                              lambda x: extract_gene_symbols(x, w2v_model))
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [1.013, 0.458, 1.472], epsilon=0.001)
