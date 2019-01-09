from typing import List
import numpy as np

from gsd.distance import DistanceMetric, calc_pairwise_distances
from gsd.gene_sets import GeneSet
import random


class RandomDistanceMetric(DistanceMetric):
    @property
    def display_name(self) -> str:
        return "Random (uniform, (0,1))"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        return calc_pairwise_distances(gene_sets, lambda a, b: random.uniform(0, 1))


def overlap_coefficient(list_a: List[bool], list_b: List[bool]) -> float:
    intersection = float(sum([a * b for a, b in zip(list_a, list_b)]))
    return intersection/min(sum(list_a), sum(list_b))


BENCHMARK_DISTS = {
    'Random_0_1': RandomDistanceMetric(),
}
