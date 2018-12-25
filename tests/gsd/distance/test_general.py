from tests import has_equal_elements

import gsd.gene_sets
from gsd.distance.general import MinkowskiNormDistanceMetric, JaccardDistanceMetric, KappaDistanceMetric

gene_sets = [
    gsd.gene_sets.GeneSet(name="SetA",
                          external_id="SetA",
                          external_source="",
                          summary="",
                          calculated=True,
                          entregene_ids={8908, 2998, 2997},
                          gene_symbols={"GYG2", "GYS2", "GYS1"}),
    gsd.gene_sets.GeneSet(name="SetB",
                          external_id="SetB",
                          external_source="",
                          summary="",
                          calculated=True,
                          entregene_ids={5507, 8908, 2998},
                          gene_symbols={"PPP1R3C", "GYG2", "GYS2"}),
    gsd.gene_sets.GeneSet(name="SetC",
                          external_id="SetC",
                          external_source="",
                          summary="",
                          calculated=True,
                          entregene_ids={2998, 2997, 2992},
                          gene_symbols={"GYS2", "GYS1", "GYG1"})
]


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
