import gsd.reactome

reactome_id = "R-HSA-8982491"
out_dir = "test/R-HSA-8982491"

human_tax_id = 9606


def test_load():
    gsd.reactome.download_reactome_sub_tree(human_tax_id, reactome_id, out_dir)
