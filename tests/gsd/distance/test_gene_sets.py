import gsd.gene_sets
import gsd.annotation

from tests.gsd.distance import go_anno


def test_load():
    gene_sets = gsd.gene_sets.load("gsd/distance/fake_gene_sets.json", go_anno)
    assert len(gene_sets) == 3
    assert gene_sets[0].general_info.name == 'Secretion of collagens'
    assert len(gene_sets[0].go_info.biological_process.ids) > 10
