human_tax_id = 9606

import gsd.reactome

rule targets:
    input:
        "evaluation_data/R-HSA-8982491",
        "evaluation_data/R-HSA-1474290"

rule download_reactome:
    output:
        directory("evaluation_data/{reactome_id}")
    run:
        gsd.reactome.download_reactome_sub_tree(human_tax_id,
                                                wildcards.reactome_id,
                                                output[0])
