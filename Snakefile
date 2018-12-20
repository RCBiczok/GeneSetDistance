import gsd.reactome

human_tax_id = 9606
reactome_sub_tree_ids = ['R-HSA-8982491', 'R-HSA-1474290']

rule all:
    input:
        "evaluation_data/R-HSA-8982491",
        "evaluation_data/R-HSA-1474290"

rule download_reactome_sub_tree:
    output:
        [directory("evaluation_data/%s" % reactome_id)
        for reactome_id in reactome_sub_tree_ids]
    run:
        for reactome_id in reactome_sub_tree_ids:
            gsd.reactome.download_reactome_sub_tree(human_tax_id, reactome_id,
                                                    "evaluation_data/%s" % reactome_id)
