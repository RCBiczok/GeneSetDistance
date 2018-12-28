from functools import reduce
from typing import List

import numpy as np
from gensim.models.keyedvectors import Word2VecKeyedVectors
from scipy.spatial.distance import cosine

from gsd.distance import DistanceMetric
from gsd.distance.general import calc_pairwise_distances
from gsd.gene_sets import GeneSet
import nltk.corpus


def filter_stop_words(sentence: List[str]) -> List[str]:
    stop_words = nltk.corpus.stopwords.words('english')
    return [w for w in sentence if w not in stop_words]


def extract_gene_symbols(gene_set: GeneSet, w2v_model: Word2VecKeyedVectors) -> List[str]:
    return [gene_sym.lower() for gene_sym in gene_set.gene_symbols if gene_sym.lower() in w2v_model.vocab]


def cosine_distance_of(words_a: List[str], words_b: List[str], w2v_model: Word2VecKeyedVectors) -> float:
    gene_syms_vec_a = reduce(lambda x, y: x + y, w2v_model[words_a])
    gene_syms_vec_b = reduce(lambda x, y: x + y, w2v_model[words_b])

    return cosine(gene_syms_vec_a, gene_syms_vec_b)


class NLPCosineDistance(DistanceMetric):
    def __init__(self, w2v_model: Word2VecKeyedVectors):
        self.w2v_model = w2v_model

    @property
    def display_name(self) -> str:
        return "Cosine Distance over Gene Symbols"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def calc_path_distance(gene_set_a: GeneSet, gene_set_b: GeneSet):
            gene_syms_a = extract_gene_symbols(gene_set_a, self.w2v_model)
            gene_syms_b = extract_gene_symbols(gene_set_b, self.w2v_model)

            return cosine_distance_of(gene_syms_a, gene_syms_b, self.w2v_model)

        return calc_pairwise_distances(gene_sets, calc_path_distance)
