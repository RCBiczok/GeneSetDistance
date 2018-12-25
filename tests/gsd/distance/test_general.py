from tests import has_equal_elements

from gsd.distance.general import MinkowskiNormDistanceMetric, JaccardDistanceMetric, KappaDistanceMetric
from tests.gsd.distance import gene_sets


def test_euclidean():
    dist_metric = MinkowskiNormDistanceMetric(2)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [1.414, 1.414, 2], epsilon=0.001)


def test_jaccard():
    dist_metric = JaccardDistanceMetric()
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.5, 0.5, 0.8], epsilon=0.001)


def test_kappa():
    dist_metric = KappaDistanceMetric()
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.833, 0.833, 1.666], epsilon=0.001)
