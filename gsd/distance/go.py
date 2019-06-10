from typing import List
import numpy as np

from gsd.distance import DistanceMetric, calc_pairwise_distances
from gsd.gene_sets import GeneSet, GOType

from rpy2.robjects.packages import importr

org_hs_en_db = importr("org.Hs.eg.db")
go_sem_sim = importr("GOSemSim")


class GOSimDistanceMetric(DistanceMetric):
    def __init__(self, go_type: GOType, measure="Wang", combine="BMA"):
        self.go_type = go_type
        self.measure = measure
        self.combine = combine
        self.hs_go_data = go_sem_sim.godata('org.Hs.eg.db', ont=go_type.value)

    @property
    def display_name(self) -> str:
        return "GO-distance (go_type=%s, measure=%s, combine=%s)" \
               % (self.go_type.value, self.measure, self.combine)

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def calc_dist(a: GeneSet, b: GeneSet):
            return 1 - go_sem_sim.mgoSim(list(self.go_type.select_category(a.go_info).ids),
                                         list(self.go_type.select_category(b.go_info).ids),
                                         self.hs_go_data,
                                         measure=self.measure,
                                         combine=self.combine)[0]

        return calc_pairwise_distances(gene_sets, calc_dist)
