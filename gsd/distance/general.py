from typing import List, TypeVar, Callable
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


T = TypeVar('T')


def calc_pairwise_distances(obj_list: List[T], dist_fun: Callable[[T, T], float]) -> np.ndarray:
    result = np.ndarray(shape=(calc_n_comparisons(obj_list),), dtype=float)
    idx = 0

    for i in range(0, len(obj_list) - 1):
        for j in range(i + 1, len(obj_list)):
            result[idx] = dist_fun(obj_list[i], obj_list[j])
            idx += 1
    return result


class KappaDistanceMetric(DistanceMetric):
    @property
    def display_name(self) -> str:
        return "Kappa distance"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        vector_list = to_binary_matrix(gene_sets)
        return calc_pairwise_distances(vector_list, lambda a, b: 1 - cohen_kappa_score(a, b))
