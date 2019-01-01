import os
import time

import jsonpickle
import numpy as np
from abc import abstractmethod
from typing import List, TypeVar, Iterable, Callable

from anytree import Node, PostOrderIter
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from tqdm import tqdm

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
                 results: Iterable[float]):
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
        out_file: str):
    os.makedirs(os.path.dirname(out_file), exist_ok=True)

    time_begin = time.time()
    d = metric.calc(gene_sets)
    time_end = time.time()

    result = EvaluationResult(metric.display_name, time_end - time_begin, d.tolist())
    with open(out_file, "w") as gene_set_file:
        gene_set_file.write(jsonpickle.encode(result))


def calc_pairwise_distances(obj_list: List[T], dist_fun: Callable[[T, T], float]) -> np.ndarray:
    result = np.ndarray(shape=(calc_n_comparisons(obj_list),), dtype=float)
    idx = 0

    for i in tqdm(range(0, len(obj_list) - 1)):
        for j in range(i + 1, len(obj_list)):
            result[idx] = dist_fun(obj_list[i], obj_list[j])
            idx += 1
    return result


class PairwiseTreePathDistanceMetric(DistanceMetric):
    def __init__(self, root: Node):
        self.root = root

    @property
    def display_name(self) -> str:
        return "Pairwise path length in reference tree"

    @abstractmethod
    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        nodes = set([gene_set.general_info.name for gene_set in gene_sets])
        nodes_mapping = dict(zip(nodes, range(0, len(nodes))))
        dist_matrix = np.zeros((len(nodes_mapping), len(nodes_mapping)))

        for node in PostOrderIter(self.root):
            if node.parent is None:
                continue
            from_idx = nodes_mapping[node.name]
            to_idx = nodes_mapping[node.parent.name]
            dist_matrix[from_idx][to_idx] = 1

        graph = csr_matrix(dist_matrix)
        dist_matrix = shortest_path(csgraph=graph, directed=False, indices=range(0, len(gene_sets)))

        def calc_path_distance(node_name_a: str, node_name_b: str):
            i = nodes_mapping[node_name_a]
            j = nodes_mapping[node_name_b]
            return dist_matrix[i][j]

        return calc_pairwise_distances([gene_set.general_info.name for gene_set in gene_sets], calc_path_distance)
