from typing import List
from scipy.spatial.distance import pdist
import numpy as np

from gsd import flat_list
from gsd.distance import DistanceMetric
from gsd.gene_sets import GeneSet
from sklearn.metrics import cohen_kappa_score


def to_binary_matrix(gene_sets: List[GeneSet]):
    all_genes = set(flat_list([gene_set.entregene_ids for gene_set in gene_sets]))

    return [[ref_gene in gene_set.entregene_ids for ref_gene in all_genes] for gene_set in gene_sets]


def calc_n_comparisons(gene_sets: List[GeneSet]) -> int:
    n = len(gene_sets) - 1
    return int(n * (n + 1) / 2)


class MinkowskiNormDistanceMetric(DistanceMetric):
    def __init__(self, p: float = 2):
        self.p = p

    def __repr__(self):
        return "%s(name=%s)" % (self.__class__.__name__, self.display_name)

    @property
    def display_name(self) -> str:
        return "Minkowski distance (p=%d)" % self.p

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        vector_list = to_binary_matrix(gene_sets)

        return pdist(np.array(vector_list), 'minkowski', self.p)


class JaccardDistanceMetric(DistanceMetric):

    def __repr__(self):
        return "%s(name=%s)" % (self.__class__.__name__, self.display_name)

    @property
    def display_name(self) -> str:
        return "Jaccard distance"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        vector_list = to_binary_matrix(gene_sets)

        return pdist(np.array(vector_list), 'jaccard')


class KappaDistanceMetric(DistanceMetric):
    def __repr__(self):
        return "%s(name=%s)" % (self.__class__.__name__, self.display_name)

    @property
    def display_name(self) -> str:
        return "Kappa distance"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        vector_list = to_binary_matrix(gene_sets)

        result = np.ndarray(shape=(calc_n_comparisons(gene_sets),), dtype=float)
        idx = 0

        for i in range(0, len(gene_sets) - 1):
            for j in range(i+1, len(gene_sets)):
                result[idx] = 1-cohen_kappa_score(vector_list[i], vector_list[j])
                idx += 1
        return result
