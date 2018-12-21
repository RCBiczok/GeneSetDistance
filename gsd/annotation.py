from typing import List

from pandas import DataFrame, read_csv
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
