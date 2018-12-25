from gsd.distance.ppi import DirectPPIDistanceMetric, load_ppi_mitab
from tests import has_equal_elements
from tests.gsd.distance import gene_sets

ppi_data = load_ppi_mitab("../__data/ppi/BioGrid/BIOGRID-ALL-3.5.166.mitab.txt")


def test_direct_ppi():
    dist_metric = DirectPPIDistanceMetric(ppi_data)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.631, 0.550, 0.833], epsilon=0.001)
