from typing import List
from scipy.spatial.distance import pdist
import numpy as np

from gsd.distance import DistanceMetric, to_binary_matrix, calc_n_comparisons
from gsd.gene_sets import GeneSet
from sklearn.metrics import cohen_kappa_score


class MinkowskiNormDistanceMetric(DistanceMetric):
    def __init__(self, p: float = 2):
        self.p = p

    @property
    def display_name(self) -> str:
        return "Minkowski distance (p=%d)" % self.p

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        vector_list = to_binary_matrix(gene_sets)

        return pdist(np.array(vector_list), 'minkowski', self.p)


class JaccardDistanceMetric(DistanceMetric):
    @property
    def display_name(self) -> str:
        return "Jaccard distance"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        vector_list = to_binary_matrix(gene_sets)

        return pdist(np.array(vector_list), 'jaccard')


class KappaDistanceMetric(DistanceMetric):
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
