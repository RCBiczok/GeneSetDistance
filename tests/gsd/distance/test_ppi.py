from pandas import read_table

from gsd.distance.ppi import DirectPPIDistanceMetric, load_ppi_mitab, ShortestPathPPI
from tests import has_equal_elements
from tests.gsd.distance import gene_sets
from tests.gsd.test_reactome import human_tax_id

ppi_data = load_ppi_mitab("../__data/ppi/BioGrid/BIOGRID-ALL-3.5.166.mitab.txt", human_tax_id)

fake_ppi_data = read_table("gsd/distance/fake_ppi.tsv")


def test_direct_ppi():
    dist_metric = DirectPPIDistanceMetric(ppi_data)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.631, 0.550, 0.833], epsilon=0.001)


def test_shortest_path_ppi():
    dist_metric = ShortestPathPPI(fake_ppi_data)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.333, 0.666, 1], epsilon=0.001)
