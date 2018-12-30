import os
import time

import jsonpickle
import numpy as np
from abc import abstractmethod
from typing import List, TypeVar

from gsd import flat_list
from gsd.gene_sets import GeneSet


class DistanceMetric:
    def __repr__(self):
        return "%s(name=%s)" % (self.__class__.__name__, self.display_name)

    @property
    @abstractmethod
    def display_name(self) -> str:
        """Returns a user-friendly name of this distance measure"""
        raise KeyError("Not implemented")

    @abstractmethod
    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        """Calculates a matrix of pairwise distances for the given gene sets"""
        raise KeyError("Not implemented")


class EvaluationResult:
    def __repr__(self):
        return "EvaluationResult(name=%s, exec_time=%f, results=%s)" \
               % (self.name, self.exec_time, self.results)

    def __init__(self,
                 name: str,
                 exec_time: float,
                 results: np.ndarray):
        self.name = name
        self.exec_time = exec_time
        self.results = results


def to_binary_matrix(gene_sets: List[GeneSet]):
    all_genes = set(flat_list([gene_set.general_info.entrez_gene_ids for gene_set in gene_sets]))

    return [[ref_gene in gene_set.general_info.entrez_gene_ids for ref_gene in all_genes] for gene_set in gene_sets]


T = TypeVar('T')


def calc_n_comparisons(gene_sets: List[T]) -> int:
    n = len(gene_sets) - 1
    return int(n * (n + 1) / 2)


def execute_and_persist_evaluation(
        metric: DistanceMetric,
        gene_sets: List[GeneSet],
        target_name: str,
        out_dir: str):
    os.makedirs(out_dir, exist_ok=True)

    time_begin = time.time()
    d = metric.calc(gene_sets)
    time_end = time.time()

    result = EvaluationResult(metric.display_name, time_end - time_begin, d)
    with open(os.path.join(out_dir, '%s.json' % target_name), "w") as gene_set_file:
        gene_set_file.write(jsonpickle.encode(result))
