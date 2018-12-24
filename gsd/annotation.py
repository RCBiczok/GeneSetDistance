from typing import List

from pandas import DataFrame, read_csv, read_table
from pandas.compat import cStringIO

from biomart import BiomartServer, BiomartDataset


def query_df(ds: BiomartDataset, params: dict) -> DataFrame:
    response = ds.search(params=params)
    return read_csv(cStringIO(response.text), sep='\t', names=params['attributes'], dtype=str)


def download_biomart_anno(attributes: List[str], out_file: str):
    server = BiomartServer("http://www.ensembl.org/biomart")
    ds = server.datasets["hsapiens_gene_ensembl"]
    df = query_df(ds, {'attributes': attributes}).dropna()
    df.to_csv(out_file, sep="\t", index=False)


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


def read_go_anno_df(entrezgene2go_file: str, go_file: str) -> DataFrame:
    entrezgene2go_df = read_table(entrezgene2go_file)
    go_df = read_table(go_file)
    return entrezgene2go_df.join(go_df.set_index("go_id"), on="go_id")