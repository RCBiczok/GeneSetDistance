from typing import Set, List, Any, Dict
import jsonpickle
from anytree import Node
from pandas import DataFrame

from anytree.importer import JsonImporter
from gsd.annotation import GOInfo


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
        return "<AnnotatedGeneSet(general_info='%s')>" % self.general_info


def annotate_with_go(gene_set_info_list: List[GeneSetInfo], go_anno: DataFrame) -> [GeneSet]:
    return [GeneSet(gene_set_info,
                    GOInfo(genes=gene_set_info.entrez_gene_ids, go_anno=go_anno))
            for gene_set_info in gene_set_info_list]


# TODO Compose GeneSet instead of modifying GeneSet instance
def load_gene_sets(gene_sets_file: str,
                   ncbi_gene_desc_file: str = None,
                   gwas_gene_traits_file: str = None) -> [GeneSet]:
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
