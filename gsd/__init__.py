import jsonpickle
import os.path
from typing import List, TypeVar

from anytree import Node
from anytree.exporter import JsonExporter

T = TypeVar('T')


def flat_list(l: List[List[T]]) -> List[T]:
    return [item for sublist in l for item in sublist]


def persist_reference_data(tree: Node, gene_sets: List, tree_out: str, gene_set_out: str):
    exporter = JsonExporter(indent=2, sort_keys=True)

    os.makedirs(os.path.basename(tree_out), exist_ok=True)
    os.makedirs(os.path.basename(gene_set_out), exist_ok=True)

    with open(tree_out, "w") as tree_file:
        exporter.write(tree, filehandle=tree_file)
    with open(gene_set_out, "w") as gene_set_file:
        gene_set_file.write(jsonpickle.encode(gene_sets))
