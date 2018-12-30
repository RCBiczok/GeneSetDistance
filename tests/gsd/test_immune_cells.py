from pandas import read_table

from gsd.immune_cells import extract_from_raw_data
from tests.gsd.distance import go_anno

the_immune_cell_data_dir = "../raw_data/immune_cells"
gene_sym_hsapiens = read_table("../annotation_data/entrezgene2gene_sym.tsv")


def test_load_cell_tree():
    node, gene_sets = extract_from_raw_data(the_immune_cell_data_dir, gene_sym_hsapiens, go_anno)
    assert node is not None
    assert list is not None
