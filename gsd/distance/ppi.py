from typing import List
from pandas import read_table, DataFrame
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from statistics import median
import copy

from gsd.distance import DistanceMetric
from gsd.distance.general import JaccardDistanceMetric, calc_pairwise_distances
from gsd.gene_sets import GeneSet


class DirectPPIDistanceMetric(JaccardDistanceMetric):
    def __init__(self, ppi_data: DataFrame):
        self.ppi_data = ppi_data

    @property
    def display_name(self) -> str:
        return "Jaccard distance over extended gene set"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def extend_gene_set(gene_set: GeneSet) -> GeneSet:
            directly_connected_genes = self.ppi_data.loc[
                self.ppi_data['FromId'].isin(gene_set.general_info.entrez_gene_ids),
                "ToId"].tolist()

            gene_set_copy = copy.deepcopy(gene_set)
            gene_set_copy.general_info.entrez_gene_ids |= set(directly_connected_genes)

            return gene_set_copy

        return super().calc([extend_gene_set(gene_set) for gene_set in gene_sets])


class ShortestPathPPI(DistanceMetric):
    def __init__(self, ppi_data: DataFrame):
        nodes = set(ppi_data['FromId'].tolist() + ppi_data['ToId'].tolist())
        self.nodes_mapping = dict(zip(nodes, range(0, len(nodes))))
        dist_matrix = np.zeros((len(self.nodes_mapping), len(self.nodes_mapping)))

        for index, row in ppi_data.iterrows():
            from_idx = self.nodes_mapping[row['FromId']]
            to_idx = self.nodes_mapping[row['ToId']]
            dist_matrix[from_idx][to_idx] = 1

            self.graph = csr_matrix(dist_matrix)

    @property
    def display_name(self) -> str:
        return "Shortest Path PPI"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def calc_path_distance(gene_set_a: GeneSet, gene_set_b: GeneSet):
            indices_a = [self.nodes_mapping[gene_id] for gene_id in gene_set_a.general_info.entrez_gene_ids]
            indices_b = [self.nodes_mapping[gene_id] for gene_id in gene_set_b.general_info.entrez_gene_ids]

            dist_matrix = shortest_path(csgraph=self.graph, directed=False, indices=indices_a)
            # return median([min([row[idx_b] for idx_b in indices_b]) for row in dist_matrix])
            if gene_set_a == gene_set_b:
                return 0
            return median([row[idx_b] for idx_b in indices_b for row in dist_matrix])

        return calc_pairwise_distances(gene_sets, calc_path_distance)


def load_ppi_mitab(ppi_file: str, tax_id) -> DataFrame:
    ppi_data = read_table(ppi_file)
    ppi_data = ppi_data.assign(
        FromId=[int(e.split(":")[1]) for e in ppi_data["#ID Interactor A"]],
        ToId=[int(e.split(":")[1]) for e in ppi_data["ID Interactor B"]],
        FromTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor A"]],
        ToTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor B"]])
    return ppi_data[(ppi_data['FromTaxID'] == tax_id) & (ppi_data['ToTaxID'] == tax_id)]
