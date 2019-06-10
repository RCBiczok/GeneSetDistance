import gsd.gene_sets
from gsd.distance.go import GOSimDistanceMetric
from tests import has_equal_elements

from tests.gsd.distance import gene_sets


def test_go_sim_anno():
    dist_metric = GOSimDistanceMetric(gsd.gene_sets.GOType.CELLULAR_COMPONENT)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.073, 0.118, 0.2], epsilon=0.001)
