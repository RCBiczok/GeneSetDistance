import os.path
from functools import reduce
from typing import TypeVar, Iterable, Set, Tuple, List, Dict

from anytree import Node
from pandas import read_excel, read_table

from gsd import flat_list
from gsd.gene_sets import GeneSet

T = TypeVar('T')


def intersect_lists(lists: Iterable[T]):
    return list(reduce(lambda a, b: set(a) & set(b), lists))


def extract_immune_cell_tree(immune_cell_data_dir: str) -> Node:
    cell_type_tree = read_excel(os.path.join(immune_cell_data_dir, "cell_type_mapping.xlsx"),
                                sheet="controlled_vocabulary")[["parent", "cell_type"]]

    node_list = {'cell': Node('cell')}

    for index, row in cell_type_tree.iterrows():
        if row['cell_type'] not in node_list:
            node_list[row['cell_type']] = Node(row['cell_type'], parent=node_list[row['parent']])

    return node_list['cell']


def filter_missing_sub_trees(node: Node, gene_set_names: Set[str]):
    filtered_children = [filter_missing_sub_trees(child, gene_set_names) for child in node.children]
    filtered_children = [children for children in filtered_children if children != []]

    if node.name not in gene_set_names and filtered_children == []:
        return list()

    node_cpy = Node(node.name)
    if filtered_children:
        node_cpy.children = filtered_children

    return node_cpy


def to_gene_set(generated_immune_marker_genes, cell_type) -> GeneSet:
    tbl_slice = generated_immune_marker_genes[generated_immune_marker_genes.cell_type == cell_type]

    return GeneSet(cell_type,
                   cell_type,
                   'Literature Review',
                   "",
                   False,
                   set(tbl_slice['entrezgene']),
                   set(tbl_slice['gene_symbol']))


def extract_genes_from(node: Node, gene_sets: Dict[str, GeneSet], cell_types_with_genes) -> List[GeneSet]:
    children_gene_sets = flat_list([extract_genes_from(children, gene_sets, cell_types_with_genes)
                                    for children in node.children])
    if node.name in cell_types_with_genes:
        return [gene_sets[node.name]] + children_gene_sets

    genes = flat_list([children_gene_set.entregene_ids for children_gene_set in children_gene_sets])
    gene_symbol = flat_list([children_gene_set.gene_symbols for children_gene_set in children_gene_sets])

    parent_gene_set = GeneSet(node.name,
                              node.name,
                              'Literature Review',
                              "",
                              True,
                              set(genes),
                              set(gene_symbol))

    return [parent_gene_set] + children_gene_sets


def extract_from_raw_data(immune_cell_data_dir: str, gene_sym_hsapiens) -> Tuple[Node, List]:
    generated_immune_marker_genes = read_table(os.path.join(immune_cell_data_dir, "signatures_all.txt"))
    generated_immune_marker_genes = generated_immune_marker_genes.join(
        gene_sym_hsapiens.set_index("external_gene_name"), on="gene_symbol").dropna()
    cell_types_with_genes = set(generated_immune_marker_genes["cell_type"].tolist())

    immune_cell_tree = extract_immune_cell_tree(immune_cell_data_dir)
    immune_cell_tree = filter_missing_sub_trees(immune_cell_tree, cell_types_with_genes)

    gene_sets = [to_gene_set(generated_immune_marker_genes, cell_type) for cell_type in cell_types_with_genes]
    gene_sets = [gene_set for gene_set in gene_sets if len(gene_set.entregene_ids) > 0]

    all_gene_sets = extract_genes_from(immune_cell_tree,
                                       {gene_set.name: gene_set for gene_set in gene_sets},
                                       cell_types_with_genes)
    return immune_cell_tree, all_gene_sets
