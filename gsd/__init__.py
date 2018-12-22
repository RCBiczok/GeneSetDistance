import json
import os.path
from typing import List

from anytree import Node
from anytree.exporter import JsonExporter


class GeneSet:
    def __init__(self, name, external_id, external_source, summary, entregene_ids, gene_symbols):
        self.name = name
        self.external_id = external_id
        self.external_source = external_source
        self.summary = summary
        self.entregene_ids = entregene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneSet(name='%s', external_id='%s')>" % (self.name, self.external_id)


def flat_list(l):
    return [item for sublist in l for item in sublist]


def persist_reference_data(tree: Node, gene_sets: List, out_dir: str):
    exporter = JsonExporter(indent=2, sort_keys=True)

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'tree.json'), "w") as tree_file:
        exporter.write(tree, filehandle=tree_file)
    with open(os.path.join(out_dir, 'gene_sets.json'), "w") as gene_set_file:
        json.dump(gene_sets, indent=2, fp=gene_set_file)
