from typing import List
from pandas import DataFrame, read_table

BIOMART_GO_ID = "go_id"
BIOMART_GO_NAME = "name_1006"
BIOMART_GO_DEFINITION = "definition_1006"
BIOMART_GO_NAMESPACE = "namespace_1003"


class GOCategory:
    def __init__(self, df: DataFrame):
        self.ids = set(df[BIOMART_GO_ID].tolist())
        self.names = set(df[BIOMART_GO_NAME].tolist())
        self.definitions = set(df[BIOMART_GO_NAME].tolist())

    def __repr__(self):
        return "<GOCategory(n_ids=%d)>" % len(self.ids)


class GOInfo:
    def __init__(self, genes: List[int], go_anno: DataFrame):
        self.molecular_function = GOCategory(go_anno[(go_anno['entrezgene'].isin(genes))
                                                     & (go_anno[BIOMART_GO_NAMESPACE] == "molecular_function")])
        self.cellular_component = GOCategory(go_anno[(go_anno['entrezgene'].isin(genes))
                                                     & (go_anno[BIOMART_GO_NAMESPACE] == "cellular_component")])
        self.biological_process = GOCategory(go_anno[(go_anno['entrezgene'].isin(genes))
                                                     & (go_anno[BIOMART_GO_NAMESPACE] == "biological_process")])

    def __repr__(self):
        return "<GOInfo(molecular_function=%s, cellular_component=%s, biological_process=%s)>" % \
               (self.molecular_function, self.cellular_component, self.biological_process)


class GeneSet:
    def __init__(self,
                 name: str,
                 external_id: str,
                 external_source: str,
                 summary: str,
                 calculated: bool,
                 entregene_ids: List[float],
                 gene_symbols: List[str]):
        self.name = name
        self.external_id = external_id
        self.external_source = external_source
        self.summary = summary
        self.calculated = calculated
        self.entregene_ids = entregene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneSet(name='%s', n_external_id='%s')>" % (self.name, len(self.external_id))


def read_go_anno_df(entrezgene2go_file: str, go_file: str) -> DataFrame:
    entrezgene2go_df = read_table(entrezgene2go_file)
    go_df = read_table(go_file)
    return entrezgene2go_df.join(go_df.set_index("go_id"), on="go_id")
