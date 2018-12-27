from typing import List
import numpy as np
from pandas import DataFrame
from rpy2.robjects.packages import importr

from gsd.annotation import GOInfo, GOType
from gsd.distance import DistanceMetric
from gsd.distance.general import calc_pairwise_distances
from gsd.gene_sets import GeneSet

org_hs_en_db = importr("org.Hs.eg.db")
go_sem_sim = importr("GOSemSim")


class GOSimDistanceMetric(DistanceMetric):
    def __init__(self, go_anno: DataFrame, go_type: GOType, measure="Wang", combine="BMA"):
        self.go_anno = go_anno
        self.go_type = go_type
        self.measure = measure
        self.combine = combine
        self.hs_go_data = go_sem_sim.godata('org.Hs.eg.db', ont=go_type.value)

    @property
    def display_name(self) -> str:
        return "GO-based distance (measure=%s, combine=%s)" % (self.measure, self.combine)

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        go_sets = [GOInfo(genes=gene_set.entrez_gene_ids, go_anno=self.go_anno) for gene_set in gene_sets]

        return calc_pairwise_distances(go_sets,
                                       lambda a, b: 1 - go_sem_sim.mgoSim(list(self.go_type.select_category(a).ids),
                                                                          list(self.go_type.select_category(b).ids),
                                                                          self.hs_go_data,
                                                                          measure=self.measure,
                                                                          combine=self.combine)[0])
