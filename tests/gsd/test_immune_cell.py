from pandas import read_table

from gsd.immune_cell import etract_immuno_cell_data

the_immune_cell_data_dir = "../../raw_data/immune_cells"
gene_sym_hsapiens = read_table("../../annotation_data/entrezgene2gene_sym.tsv")


def test_load_cell_tree():
    node, gene_sets = etract_immuno_cell_data(the_immune_cell_data_dir, gene_sym_hsapiens)
    assert node is not None
    assert list is not None
