import json
import urllib.request

import jsonpickle
from typing import List, Set, Dict, Any
from pandas import DataFrame, read_csv, read_table
from pandas.compat import cStringIO
from biomart import BiomartServer, BiomartDataset
from enum import Enum, unique

from tqdm import tqdm

from gsd import flat_list

BIOMART_GO_ID = "go_id"
BIOMART_GO_NAME = "name_1006"
BIOMART_GO_DEFINITION = "definition_1006"
BIOMART_GO_NAMESPACE = "namespace_1003"
BIOMART_GO_LINKAGE_TYPE = 'go_linkage_type'


class GOCategory:
    def __init__(self, df: DataFrame):
        self.ids = set(df[BIOMART_GO_ID].tolist())
        self.names = set(df[BIOMART_GO_NAME].tolist())
        self.definitions = set(df[BIOMART_GO_DEFINITION].tolist())

    def __repr__(self):
        return "<GOCategory(n_ids=%d)>" % len(self.ids)


class GOInfo:
    def __init__(self, genes: Set[int], go_anno: DataFrame):
        self.molecular_function = GOCategory(go_anno[(go_anno['entrezgene'].isin(genes))
                                                     & (go_anno[BIOMART_GO_NAMESPACE] == "molecular_function")])
        self.cellular_component = GOCategory(go_anno[(go_anno['entrezgene'].isin(genes))
                                                     & (go_anno[BIOMART_GO_NAMESPACE] == "cellular_component")])
        self.biological_process = GOCategory(go_anno[(go_anno['entrezgene'].isin(genes))
                                                     & (go_anno[BIOMART_GO_NAMESPACE] == "biological_process")])

    def __repr__(self):
        return "<GOInfo(molecular_function=%s, cellular_component=%s, biological_process=%s)>" % \
               (self.molecular_function, self.cellular_component, self.biological_process)


@unique
class GOType(Enum):
    MOLECULAR_FUNCTION = "MF"
    CELLULAR_COMPONENT = "CC"
    BIOLOGICAL_PROCESS = "BP"

    def select_category(self, go_info: GOInfo) -> GOCategory:
        if self.value == "MF":
            return go_info.molecular_function
        elif self.value == "CC":
            return go_info.cellular_component
        return go_info.biological_process


class NCBIGeneInfo:
    def __init__(self,
                 gene_set_name: str,
                 gene_infos: Dict[str, Any]):
        self.gene_set_name = gene_set_name
        self.gene_infos = gene_infos

    def __repr__(self):
        return "<NCBIGeneInfo(gene_set_name='%s', gene_infos='%s')>" % (self.gene_set_name, self.gene_infos)


def read_go_anno_df(entrezgene2go_file: str, go_file: str) -> DataFrame:
    entrezgene2go_df = read_table(entrezgene2go_file)
    go_df = read_table(go_file)
    return entrezgene2go_df.join(go_df.set_index("go_id"), on="go_id")


def query_df(ds: BiomartDataset, params: dict) -> DataFrame:
    response = ds.search(params=params)
    return read_csv(cStringIO(response.text), sep='\t', names=params['attributes'], dtype=str)


def download_biomart_anno(attributes: List[str], out_file: str):
    server = BiomartServer("http://www.ensembl.org/biomart")
    ds = server.datasets["hsapiens_gene_ensembl"]
    df = query_df(ds, {'attributes': attributes}).dropna()
    df.to_csv(out_file, sep="\t", index=False)


def get_json_from(url):
    with urllib.request.urlopen(url) as con:
        data = json.loads(con.read().decode())
    return data


def chunks(l: List, n: int):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def _get_ncbi_gene_dscr(entrezgene_list: List[int]):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s&retmode=json;' % \
          ','.join([str(x) for x in entrezgene_list])
    return get_json_from(url)


def _get_all_ncbi_gene_dscr(entrezgene_list: List[int]):
    chunk_list = [chunk for chunk in chunks(entrezgene_list, 300)]
    result = {}
    for chunk in tqdm(chunk_list):
        part_result = _get_ncbi_gene_dscr(chunk)['result']
        result = {**result, **part_result}
    return result


def create_gene_info(gene_set, data):
    genes = {int(gene_id): entry for gene_id, entry in data.items()
             if int(gene_id) in gene_set.general_info.entrez_gene_ids}

    return NCBIGeneInfo(gene_set.general_info.name, genes)


def downlaod_ncbi_gene_desc(gene_sets, ncbi_gene_desc_file: str):
    target_genes = set(flat_list([gene_set.general_info.entrez_gene_ids for gene_set in gene_sets]))
    data = _get_all_ncbi_gene_dscr(list(target_genes))
    del data['uids']

    gene_info_list = [create_gene_info(gene_set, data) for gene_set in gene_sets]
    with open(ncbi_gene_desc_file, "w") as out_file:
        out_file.write(jsonpickle.encode(gene_info_list))
