import numpy as np
from abc import abstractmethod
from typing import List

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


def to_binary_matrix(gene_sets: List[GeneSet]):
    all_genes = set(flat_list([gene_set.entrez_gene_ids for gene_set in gene_sets]))

    return [[ref_gene in gene_set.entrez_gene_ids for ref_gene in all_genes] for gene_set in gene_sets]


def calc_n_comparisons(gene_sets: List[GeneSet]) -> int:
    n = len(gene_sets) - 1
    return int(n * (n + 1) / 2)
