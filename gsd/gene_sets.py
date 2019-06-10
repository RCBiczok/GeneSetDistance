import json
import urllib.request
from enum import unique, Enum
from typing import Set, List, Any, Dict
import jsonpickle
from anytree import Node
from biomart import BiomartDataset, BiomartServer
from pandas import DataFrame, read_table, read_csv

from anytree.importer import JsonImporter
from pandas.compat import cStringIO
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


class GWASGeneTraitInfo:
    def __init__(self,
                 gene_set_name: str,
                 gene_traits: Dict[str, List[str]]):
        self.gene_set_name = gene_set_name
        self.gene_traits = gene_traits

    def __repr__(self):
        return "<GWASGeneTraitInfo(gene_set_name='%s', gene_traits='%s')>" % (self.gene_set_name, self.gene_traits)


class GeneSetInfo:
    def __init__(self,
                 name: str,
                 external_id: str,
                 external_source: str,
                 summary: str,
                 calculated: bool,
                 entrez_gene_ids: Set[int],
                 gene_symbols: Set[str]):
        self.name = name
        self.external_id = external_id
        self.external_source = external_source
        self.summary = summary
        self.calculated = calculated
        self.entrez_gene_ids = entrez_gene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneralInfo(name='%s', n_entrez_gene_ids='%s')>" % (self.name, len(self.entrez_gene_ids))


class GeneSet:
    def __init__(self,
                 general_info: GeneSetInfo,
                 go_info: GOInfo,
                 ncbi_gene_desc: Dict[str, Any] = None,
                 gwas_gene_traigs: Dict[str, List[str]] = None):
        self.general_info = general_info
        self.go_info = go_info
        self.ncbi_gene_desc = ncbi_gene_desc
        self.gwas_gene_traigs = gwas_gene_traigs

    def __repr__(self):
        return "<GeneSet(general_info='%s')>" % self.general_info


def annotate_with_go(gene_set_info_list: List[GeneSetInfo], go_anno: DataFrame) -> [GeneSet]:
    return [GeneSet(gene_set_info,
                    GOInfo(genes=gene_set_info.entrez_gene_ids, go_anno=go_anno))
            for gene_set_info in gene_set_info_list]


# TODO Compose GeneSet instead of modifying GeneSet instance
def load_gene_sets(gene_sets_file: str,
                   ncbi_gene_desc_file: str = None,
                   gwas_gene_traits_file: str = None) -> List[GeneSet]:
    with open(gene_sets_file) as f:
        gene_sets = jsonpickle.decode(f.read())

    if ncbi_gene_desc_file is not None:
        with open(ncbi_gene_desc_file) as f:
            gene_set_anno = jsonpickle.decode(f.read())

        gene_set_anno = {elem.gene_set_name: elem for elem in gene_set_anno}

        for gene_set in gene_sets:
            if gene_set.general_info.name not in gene_set_anno:
                raise KeyError("Gene set not found in gene annotation file: %s" % gene_set.general_info.name)
            gene_set.ncbi_gene_desc = gene_set_anno[gene_set.general_info.name]

    if gwas_gene_traits_file is not None:
        with open(gwas_gene_traits_file) as f:
            gene_set_anno = jsonpickle.decode(f.read())

        gene_set_anno = {elem.gene_set_name: elem for elem in gene_set_anno}

        for gene_set in gene_sets:
            if gene_set.general_info.name not in gene_set_anno:
                raise KeyError("Gene set not found in gene annotation file: %s" % gene_set.general_info.name)
            gene_set.gwas_gene_traigs = gene_set_anno[gene_set.general_info.name]

    return gene_sets


def load_tree(tree_file: str) -> Node:
    importer = JsonImporter()
    with open(tree_file) as f:
        json_str = f.read()
        return importer.import_(json_str)


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


def _get_ncbi_gene_desc(entrezgene_list: List[int]):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s&retmode=json;' % \
          ','.join([str(x) for x in entrezgene_list])
    return get_json_from(url)


def _get_all_ncbi_gene_dscr(entrezgene_list: List[int]):
    chunk_list = [chunk for chunk in chunks(entrezgene_list, 300)]
    result = {}
    for chunk in tqdm(chunk_list):
        part_result = _get_ncbi_gene_desc(chunk)['result']
        result = {**result, **part_result}
    return result


def create_gene_info(gene_set: GeneSet, data):
    genes = {int(gene_id): entry for gene_id, entry in data.items()
             if int(gene_id) in gene_set.general_info.entrez_gene_ids}

    return NCBIGeneInfo(gene_set.general_info.name, genes)


def downlaod_ncbi_gene_desc(gene_sets: List[GeneSet], ncbi_gene_desc_file: str):
    target_genes = set(flat_list([gene_set.general_info.entrez_gene_ids for gene_set in gene_sets]))
    data = _get_all_ncbi_gene_dscr(list(target_genes))
    del data['uids']

    gene_info_list = [create_gene_info(gene_set, data) for gene_set in gene_sets]
    with open(ncbi_gene_desc_file, "w") as out_file:
        out_file.write(jsonpickle.encode(gene_info_list))


def create_gwas_traits_info(gene_set: GeneSet, gwas_file: DataFrame):
    tbl_slice = gwas_file.loc[gwas_file['MAPPED_GENE'].isin(gene_set.general_info.gene_symbols),
                              ["DISEASE/TRAIT", 'MAPPED_GENE']]

    mapping = {idx: row for idx, row in
               tbl_slice.groupby('MAPPED_GENE')['DISEASE/TRAIT'].apply(list).iteritems()}

    return GWASGeneTraitInfo(gene_set.general_info.name, mapping)


def extract_gwas_traits(gwas_file: str, gene_sets: List[GeneSet], gwas_mapping_out: str):
    gwas_file = read_table(gwas_file)
    gene_traits_list = [create_gwas_traits_info(gene_set, gwas_file) for gene_set in gene_sets]

    with open(gwas_mapping_out, "w") as out_file:
        out_file.write(jsonpickle.encode(gene_traits_list))
