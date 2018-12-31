import re
from functools import reduce
from pprint import pprint
from typing import List, Callable

import numpy as np
from gensim.models.keyedvectors import Word2VecKeyedVectors
from scipy.spatial.distance import cosine

from gsd import flat_list
from gsd.annotation import GOType
from gsd.distance import DistanceMetric
from gsd.distance.general import calc_pairwise_distances
from gsd.gene_sets import GeneSet
import nltk.corpus


def _filter_stop_words(sentence: List[str]) -> List[str]:
    stop_words = nltk.corpus.stopwords.words('english')
    return [w for w in sentence if w not in stop_words]


def _extract_words_from_text(text: str) -> List[str]:
    return [word.lower() for word in re.split(r" |<br>|-|\.|\(|\)|,", text) if len(word) > 1]


def _filter_by_vocabulary(words: List[str], w2v_model: Word2VecKeyedVectors) -> List[str]:
    return [word for word in words if word in w2v_model.vocab]


def _extract_and_filter_words_from_text(text: str, w2v_model: Word2VecKeyedVectors) -> List[str]:
    raw_list = _extract_words_from_text(text)
    list_without_stopwords = _filter_stop_words(raw_list)
    return _filter_by_vocabulary(list_without_stopwords, w2v_model)


# Data Extraction

def extract_gene_symbols(gene_set: GeneSet, w2v_model: Word2VecKeyedVectors) -> List[str]:
    return _filter_by_vocabulary([gene_sym.lower() for gene_sym in gene_set.general_info.gene_symbols], w2v_model)


def extract_words_from_gene_set_summary(gene_set: GeneSet, w2v_model: Word2VecKeyedVectors) -> List[str]:
    return _extract_and_filter_words_from_text(gene_set.general_info.summary, w2v_model)


def extract_words_from_go_descriptions(
        gene_set: GeneSet,
        w2v_model: Word2VecKeyedVectors,
        go_filters: List[GOType]) -> List[str]:
    descriptions = flat_list([go_filter.select_category(gene_set.go_info).definitions for go_filter in go_filters])
    return flat_list([_extract_and_filter_words_from_text(description, w2v_model) for description in descriptions])


def extract_words_from_go_names(
        gene_set: GeneSet,
        w2v_model: Word2VecKeyedVectors,
        go_filters: List[GOType]) -> List[str]:
    names = flat_list([go_filter.select_category(gene_set.go_info).names for go_filter in go_filters])
    return flat_list([_extract_and_filter_words_from_text(name, w2v_model) for name in names])


def extract_words_from_gene_symbols_and_summary_and_go_info(
        gene_set: GeneSet,
        w2v_model: Word2VecKeyedVectors) -> List[str]:
    gene_symbols = list(extract_gene_symbols(gene_set, w2v_model))
    summary_words = extract_words_from_gene_set_summary(gene_set, w2v_model)
    go_name_words = extract_words_from_go_names(gene_set, w2v_model, [t for t in GOType])
    go_description_words = extract_words_from_go_descriptions(gene_set, w2v_model, [t for t in GOType])
    return gene_symbols + summary_words + go_name_words + go_description_words


# Distance similarities

def cosine_distance_of(words_a: List[str], words_b: List[str], w2v_model: Word2VecKeyedVectors) -> float:
    vectors_a = [w2v_model[word] for word in words_a]
    vectors_b = [w2v_model[word] for word in words_b]

    if len(vectors_a) == 0 or len(vectors_b) == 0:
        return np.inf

    return cosine(reduce(lambda x, y: x + y, vectors_a), reduce(lambda x, y: x + y, vectors_b))


def wm_distance_of(words_a: List[str], words_b: List[str], w2v_model: Word2VecKeyedVectors) -> float:
    return w2v_model.wmdistance(words_a, words_b)


# Distance implementations

class NLPDistance(DistanceMetric):
    def __init__(self,
                 name: str,
                 comparator: Callable[[List[str], List[str]], float],
                 extractor: Callable[[GeneSet], List[str]]):
        self.name = name
        self.comparator = comparator
        self.extractor = extractor

    @property
    def display_name(self) -> str:
        return self.name

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def calc_path_distance(gene_set_a: GeneSet, gene_set_b: GeneSet):
            gene_syms_a = self.extractor(gene_set_a)
            gene_syms_b = self.extractor(gene_set_b)

            return self.comparator(gene_syms_a, gene_syms_b)

        return calc_pairwise_distances(gene_sets, calc_path_distance)
