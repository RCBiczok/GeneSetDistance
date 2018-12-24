import json

from gsd.gene_sets import GOInfo, read_go_anno_df

with open("../evaluation_data/immune_cells/gene_sets.json") as fp:
    gene_sets = json.load(fp)

go_anno = read_go_anno_df("../annotation_data/entrezgene2go.tsv", "../annotation_data/go.tsv")


def test_load_cell_tree():
    go_info = GOInfo(gene_sets[0]['genes'], go_anno)
    assert len(go_info.molecular_function.ids) > 100
    assert len(go_info.cellular_component.ids) > 100
    assert len(go_info.biological_process.ids) > 100
