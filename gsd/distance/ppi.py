from typing import List
from pandas import read_table, DataFrame
import numpy as np

from gsd.distance.general import JaccardDistanceMetric
from gsd.gene_sets import GeneSet


class DirectPPIDistanceMetric(JaccardDistanceMetric):
    def __init__(self, ppi_data: DataFrame):
        self.ppi_data = ppi_data

    @property
    def display_name(self) -> str:
        return "Jaccard distance over extended gene set"

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        def extend_gene_set(gene_set: GeneSet) -> GeneSet:
            directly_connected_genes = self.ppi_data.loc[self.ppi_data['FromId'].isin(gene_set.entrez_gene_ids),
                                                         "ToId"].tolist()
            return GeneSet(gene_set.name,
                           gene_set.external_id,
                           gene_set.external_source,
                           gene_set.summary,
                           gene_set.calculated,
                           set(list(gene_set.entrez_gene_ids) + directly_connected_genes),
                           gene_set.gene_symbols)

        return super().calc([extend_gene_set(gene_set) for gene_set in gene_sets])


def load_ppi_mitab(ppi_file: str) -> DataFrame:
    ppi_data = read_table(ppi_file)
    ppi_data = ppi_data.assign(
        FromId=[int(e.split(":")[1]) for e in ppi_data["#ID Interactor A"]],
        ToId=[int(e.split(":")[1]) for e in ppi_data["ID Interactor B"]],
        FromTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor A"]],
        ToTaxID=[int(e.split(":")[1]) for e in ppi_data["Taxid Interactor B"]])
    return ppi_data
