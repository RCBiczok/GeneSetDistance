import json
import sys
import urllib.request
from functools import reduce
from typing import Tuple, List, Set
from anytree import Node
from pandas import DataFrame

from gsd import flat_list
from gsd.gene_sets import GeneSetInfo, annotate_with_go, GeneSet


def _get_json_from(url):
    with urllib.request.urlopen(url) as con:
        data = json.loads(con.read().decode())
    return data


def _get_event_hierarchy(tax_id):
    url = 'https://reactome.org/ContentService/data/eventsHierarchy/%s' % tax_id
    return _get_json_from(url)


def _get_reactome_information(reactome_id):
    url = 'https://reactome.org/ContentService/data/query/%s' % reactome_id
    return _get_json_from(url)


def _get_reactome_reference_entities(reactome_id):
    url = 'https://www.reactome.org/ContentService/data/participants/%s/referenceEntities' % reactome_id
    return _get_json_from(url)


def _get_node_by_reactome_id(reactome_node, reactome_id):
    if reactome_node["stId"] == reactome_id:
        return reactome_node

    if 'children' in reactome_node:
        for child in reactome_node['children']:
            match = _get_node_by_reactome_id(child, reactome_id)
            if match is not None:
                return match

    return None


def _extract_entrezgene_ident(gene_product):
    entrezgene_prefix = "EntrezGene:"
    if 'otherIdentifier' not in gene_product:
        return []
    return [int(ident[len(entrezgene_prefix):]) for ident in gene_product['otherIdentifier'] if
            ident.startswith(entrezgene_prefix)]


def _extract_reactome_gene_set(reactome_id) -> GeneSetInfo:
    reference_entities = _get_reactome_reference_entities(reactome_id)

    gene_products = [elem for elem in reference_entities if elem["className"] == "ReferenceGeneProduct"]

    entrezgene_id_list = [_extract_entrezgene_ident(gene_product) for gene_product in gene_products]

    symbol_list = {gene_product["name"][0] for gene_product in gene_products}

    if any([len(ids) > 1 for ids in entrezgene_id_list]):
        print("%s has multiple entrezgene IDs for the same gene" % reactome_id, file=sys.stderr)

    if any([ids == [] for ids in entrezgene_id_list]):
        print("%s has gene products without an entrezgene ID" % reactome_id, file=sys.stderr)

    reactome_genes = flat_list(entrezgene_id_list)

    if len(reactome_genes) != len(set(reactome_genes)):
        print("%s contains redundant entrezegene IDs" % reactome_id, file=sys.stderr)

    reactome_genes = set(reactome_genes)

    reactome_info = _get_reactome_information(reactome_id)

    reactome_name = reactome_info['displayName']
    reactome_summation = reactome_info['summation']

    if len(reactome_summation) > 1:
        print("%s has more than one summary" % reactome_id, file=sys.stderr)

    if len(reactome_summation) == 0:
        print("%s has no summary" % reactome_id, file=sys.stderr)

    reactome_summary = reduce(lambda a, b: a + " <br><br><br> " + b, [x['text'] for x in reactome_summation])

    return GeneSetInfo("%s (%s)" % (reactome_name, reactome_info['stId']),
                       reactome_info['stId'],
                       "Reactome",
                       reactome_summary,
                       False,
                       reactome_genes,
                       symbol_list)


def _dump_gene_sets(reactome_node) -> List[GeneSetInfo]:
    gene_set = [_extract_reactome_gene_set(reactome_node['stId'])]
    if 'children' not in reactome_node:
        return gene_set
    return gene_set + reduce(lambda a, b: a + b, [_dump_gene_sets(node) for node in reactome_node['children']])


def _reactome_to_anytree(node, filters: Set[str], parent: Node = None) -> Node:
    n = Node(name="%s (%s)" % (node['name'], node['stId']), parent=parent)
    if 'children' in node:
        for child in node['children']:
            if child['stId'] not in filters:
                _reactome_to_anytree(child, filters, n)
    return n


def download(tax_id: int, reactome_id: str, go_anno: DataFrame) -> Tuple[Node, List[GeneSet]]:
    reactome_tree = _get_event_hierarchy(tax_id)

    reactome_pseudo_tree = {'stId': "FAKE", 'name': "PseudoRoot", 'children': reactome_tree}

    sub_tree = _get_node_by_reactome_id(reactome_pseudo_tree, reactome_id)

    gene_sets = _dump_gene_sets(sub_tree)
    reactome_ids = [gene_set.external_id for gene_set in gene_sets]
    duplicated_ids = set([x for x in reactome_ids if reactome_ids.count(x) > 1])

    if len(duplicated_ids) > 0:
        print("%s has %d redundant child pathways. Gene sets %s are removed from analysis" %
              (reactome_id, len(duplicated_ids), ','.join(duplicated_ids)),
              file=sys.stderr)

    unique_gene_sets = [gene_set for gene_set in gene_sets if gene_set.external_id not in duplicated_ids]

    root = _reactome_to_anytree(sub_tree, duplicated_ids)

    return root, annotate_with_go(unique_gene_sets, go_anno)
