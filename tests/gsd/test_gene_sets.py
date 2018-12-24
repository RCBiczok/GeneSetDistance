import gsd.gene_sets


go_anno = gsd.gene_sets.read_go_anno_df("../annotation_data/entrezgene2go.tsv", "../annotation_data/go.tsv")


def test_load_cell_tree():
    entrezgene_ids = [8192, 8200, 22]
    go_info = gsd.gene_sets.GOInfo(entrezgene_ids, go_anno)
    assert len(go_info.molecular_function.ids) >= 10
    assert len(go_info.cellular_component.ids) >= 10
    assert len(go_info.biological_process.ids) >= 10
