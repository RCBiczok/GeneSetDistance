from typing import List
from pandas import read_table, DataFrame
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path
from statistics import mean

from scipy.spatial.distance import pdist
from tqdm import tqdm

from gsd.distance import DistanceMetric, calc_pairwise_distances
from gsd.distance.general import to_binary_matrix
from gsd.gene_sets import GeneSet


class DirectPPIDistanceMetric(DistanceMetric):
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

            return gene_set.general_info.entrez_gene_ids | set(directly_connected_genes)

        extended_id_map = {gene_set.general_info.name: extend_gene_set(gene_set) for gene_set in gene_sets}

        return pdist(np.array(to_binary_matrix(extended_id_map)), 'jaccard')


class ShortestPathPPI(DistanceMetric):
    def __init__(self, ppi_data: DataFrame):
        nodes = set(ppi_data['FromId'].tolist() + ppi_data['ToId'].tolist())
        self.nodes_mapping = dict(zip(nodes, range(0, len(nodes))))
        dist_matrix = np.zeros((len(self.nodes_mapping), len(self.nodes_mapping)))

        for index, row in tqdm(ppi_data.iterrows(), total=ppi_data.shape[0]):
            from_idx = self.nodes_mapping[row['FromId']]
            to_idx = self.nodes_mapping[row['ToId']]
            dist_matrix[from_idx][to_idx] = 1

        self.graph = csr_matrix(dist_matrix)

    @property
    def display_name(self) -> str:
        return "Dijkstra BMA PPI"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def calc_path_distance(gene_set_a: GeneSet, gene_set_b: GeneSet):
            indices_a = [self.nodes_mapping[gene_id] for gene_id in gene_set_a.general_info.entrez_gene_ids
                         if gene_id in self.nodes_mapping]
            indices_b = [self.nodes_mapping[gene_id] for gene_id in gene_set_b.general_info.entrez_gene_ids
                         if gene_id in self.nodes_mapping]

            if len(indices_a) == 0 or len(indices_b) == 0:
                return np.nan

            dist_matrix = shortest_path(csgraph=self.graph, directed=False, indices=indices_a)
            return mean([min([row[idx_b] for idx_b in indices_b]) for row in dist_matrix])

        return calc_pairwise_distances(gene_sets, calc_path_distance)


def load_ppi_mitab(ppi_file: str, tax_id) -> DataFrame:
    ppi_data = read_table(ppi_file)
    ppi_data = ppi_data.assign(
        FromId=[int(e.split(":")[1]) for e in ppi_data["#ID Interactor A"]],
        ToId=[int(e.split(":")[1]) for e in ppi_data["ID Interactor B"]],
        FromTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor A"]],
        ToTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor B"]])
    return ppi_data[(ppi_data['FromTaxID'] == tax_id) & (ppi_data['ToTaxID'] == tax_id)]


PPI_DISTS = {
    'Direct_PPI': lambda ppi_data: DirectPPIDistanceMetric(ppi_data),
    'Dijkstra_BMA_PPI': lambda ppi_data: ShortestPathPPI(ppi_data)
}
