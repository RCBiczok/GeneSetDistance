import numpy as np
from abc import abstractmethod
from typing import List

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
