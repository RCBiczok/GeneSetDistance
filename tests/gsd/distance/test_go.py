from gsd.distance.go import GOSimMetric
from tests import has_equal_elements

import gsd.annotation
from tests.gsd.distance import gene_sets

go_anno = gsd.annotation.read_go_anno_df("../annotation_data/entrezgene2go.tsv", "../annotation_data/go.tsv")


def test_go_sim_anno():
    dist_metric = GOSimMetric(go_anno, gsd.annotation.GOType.CELLULAR_COMPONENT)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.073, 0.118, 0.2], epsilon=0.001)
