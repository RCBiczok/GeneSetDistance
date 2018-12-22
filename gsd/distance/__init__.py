from abc import abstractmethod
from typing import List

from pandas import DataFrame

from gsd import GeneSet


class DistanceMeasure:
    def __repr__(self):
        return "%s(name=%s)" % (self.__class__.__name__, self.display_name)

    @property
    @abstractmethod
    def display_name(self) -> str:
        """Returns a user-friendly name of this distance measure"""
        raise KeyError("Not implemented")

    @property
    @abstractmethod
    def calc(self: List[GeneSet]) -> DataFrame:
        """Calculates a matrix of pairwise distances for the given gene sets"""
        raise KeyError("Not implemented")
