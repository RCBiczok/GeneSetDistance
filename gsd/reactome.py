import json
import sys
import os.path

import urllib.request
from functools import reduce


def get_json_from(url):
    with urllib.request.urlopen(url) as con:
        data = json.loads(con.read().decode())
    return data


def get_event_hierarchy(tax_id):
    url = 'https://reactome.org/ContentService/data/eventsHierarchy/%s' % tax_id
    return get_json_from(url)


def get_reactome_information(reactome_id):
    url = 'https://reactome.org/ContentService/data/query/%s' % reactome_id
    return get_json_from(url)


def get_reactome_reference_entities(reactome_id):
    url = 'https://www.reactome.org/ContentService/data/participants/%s/referenceEntities' % reactome_id
    return get_json_from(url)


def get_node_by_reactome_id(reactome_node, reactome_id):
    if reactome_node["stId"] == reactome_id:
        return reactome_node

    if 'children' in reactome_node:
        for child in reactome_node['children']:
            match = get_node_by_reactome_id(child, reactome_id)
            if match is not None:
                return match

    return None


def flat_list(l):
    return [item for sublist in l for item in sublist]


def extract_reactome_data(reactome_id):
    reference_entities = get_reactome_reference_entities(reactome_id)

    gene_products = [elem for elem in reference_entities if elem["className"] == "ReferenceGeneProduct"]

    entrezgene_prefix = "EntrezGene:"
    entregene_id_list = [[int(ident[len(entrezgene_prefix):]) for ident in gene_product['otherIdentifier'] if
                          ident.startswith(entrezgene_prefix)]
                         for gene_product in gene_products]

    symbol_list = [gene_product["name"][0] for gene_product in gene_products]

    if any([len(ids) > 1 for ids in entregene_id_list]):
        print("%s has multiple entregene IDs for the same gene" % reactome_id, file=sys.stderr)

    if any([ids == [] for ids in entregene_id_list]):
        print("%s has gene products without an entregene ID" % reactome_id, file=sys.stderr)

    reactome_genes = flat_list(entregene_id_list)

    if len(reactome_genes) != len(set(reactome_genes)):
        print("%s contains redundant entrezegene IDs" % reactome_id, file=sys.stderr)

    reactome_genes = list(set(reactome_genes))

    reactome_info = get_reactome_information(reactome_id)

    reactome_name = reactome_info['displayName']
    reactome_summation = reactome_info['summation']

    if len(reactome_summation) > 1:
        print("%s has more than one summary" % reactome_id, file=sys.stderr)

    if len(reactome_summation) == 0:
        print("%s has no summary" % reactome_id, file=sys.stderr)

    reactome_summary = reduce(lambda a, b: a + " --- " + b, [x['text'] for x in reactome_summation])

    geneset_info = {'name': reactome_name,
                    'externalId': reactome_info['stId'],
                    'externalSource': "Reactome",
                    'calculated': False,
                    'summary': reactome_summary,
                    'genes': reactome_genes,
                    'gene_symbols': symbol_list}
    return geneset_info


def dump_gene_sets(reactome_node):
    gene_set = [extract_reactome_data(reactome_node['stId'])]
    if 'children' not in reactome_node:
        return gene_set
    return gene_set + reduce(lambda a, b: a + b, [dump_gene_sets(node) for node in reactome_node['children']])


def download_reactome_sub_tree(tax_id : int, reactome_id: str, out_dir: str):
    reactome_tree = get_event_hierarchy(tax_id)

    reactome_pseudo_tree = {'stId': "FAKE", 'name': "PseudoRoot", 'children': reactome_tree}

    sub_tree = get_node_by_reactome_id(reactome_pseudo_tree, reactome_id)

    gene_sets = dump_gene_sets(sub_tree)
    unique_ids = set([gene_set['externalId'] for gene_set in gene_sets])
    unique_gene_sets = [gene_set for gene_set in gene_sets if gene_set['externalId'] in unique_ids]

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'tree.json'), "w") as tree_file:
        json.dump(sub_tree, indent=2, fp=tree_file)
    with open(os.path.join(out_dir, 'gene_sets.json'), "w") as gene_set_file:
        json.dump(unique_gene_sets, indent=2, fp=gene_set_file)