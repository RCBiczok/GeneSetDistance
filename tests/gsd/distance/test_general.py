from anytree import Node
from scipy.spatial.distance import pdist

from gsd.distance import PairwiseTreePathDistanceMetric
from tests import has_equal_elements

from gsd.distance.general import overlap_coefficient, to_binary_matrix, to_gene_id_map, to_gene_trait_map, \
    MatrixBasedDistanceMetric, kappa_distance, overlap_distance, to_gene_trait_freq, to_freq_matrix
from tests.gsd.distance import gene_sets


def test_euclidean():
    dist_metric = MatrixBasedDistanceMetric("Minkowski distance (p=2)",
                                            lambda x: to_binary_matrix(to_gene_id_map(x)),
                                            lambda y: pdist(y, 'minkowski', 2))
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [1.414, 1.414, 2], epsilon=0.001)


def test_jaccard():
    dist_metric = MatrixBasedDistanceMetric("Jaccard Distance",
                                            lambda x: to_binary_matrix(to_gene_id_map(x)),
                                            lambda y: pdist(y, 'jaccard'))
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.5, 0.5, 0.8], epsilon=0.001)


def test_kappa():
    dist_metric = MatrixBasedDistanceMetric("Kappa distance over genes",
                                            lambda x: to_binary_matrix(to_gene_id_map(x)),
                                            kappa_distance)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.833, 0.833, 1.666], epsilon=0.001)


def test_overlap_coefficient():
    assert overlap_coefficient([True, False, False], [True, True, False]) == 1
    assert overlap_coefficient([True, True, False], [True, True, False]) == 1
    assert overlap_coefficient([True, True, True], [True, True, False]) == 1
    assert overlap_coefficient([False, False, True], [True, True, False]) == 0
    assert overlap_coefficient([False, True, True], [True, True, False]) == 0.5


def test_overlap_distance():
    dist_metric = MatrixBasedDistanceMetric("Overlap distance over genes",
                                            lambda x: to_binary_matrix(to_gene_id_map(x)),
                                            overlap_distance)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [0.333, 0.333, 0.666], epsilon=0.001)


def test_pairwise_tree_path_distanc():
    root = Node(gene_sets[0].general_info.name)
    Node(gene_sets[1].general_info.name, parent=root)
    Node(gene_sets[2].general_info.name, parent=root)

    dist_metric = PairwiseTreePathDistanceMetric(root)
    d = dist_metric.calc(gene_sets)
    assert has_equal_elements(d, [1, 1, 2], epsilon=0.001)


def test_to_gene_id_map():
    id_mapping = to_gene_id_map(gene_sets)
    assert id_mapping == {'SetA': {8908, 2997, 2998},
                          'SetB': {5507, 8908, 2998},
                          'SetC': {2992, 2997, 2998}}


def test_to_gene_trait_map():
    id_mapping = to_gene_trait_map(gene_sets)
    assert id_mapping == {'SetA': {'Cough in response to angiotensin-converting enzyme inhibitor drugs',
                                   'Type 2 diabetes'},
                          'SetB': {'Cough in response to angiotensin-converting enzyme inhibitor drugs',
                                   'Type 2 diabetes'},
                          'SetC': {'Cough in response to angiotensin-converting enzyme inhibitor drugs',
                                   'Type 2 diabetes'}}


def test_to_gene_trait_dist_freq():
    id_mapping = to_gene_trait_freq(gene_sets)
    assert id_mapping == {'SetA': {'Cough in response to angiotensin-converting enzyme inhibitor drugs': 1,
                                   'Type 2 diabetes': 1},
                          'SetB': {'Cough in response to angiotensin-converting enzyme inhibitor drugs': 1,
                                   'Type 2 diabetes': 1},
                          'SetC': {'Cough in response to angiotensin-converting enzyme inhibitor drugs': 1,
                                   'Type 2 diabetes': 1}}


def test_to_binary_matrix():
    id_mapping = to_gene_id_map(gene_sets)
    assert to_binary_matrix(id_mapping) == [[False, True, False, True, True],
                                            [True, True, False, False, True],
                                            [False, False, True, True, True]]


def test_to_freq_matrix():
    freqs = {
        'SetA': {
            'TraitA': 2,
            'TraitB': 1,
            'TraitC': 3
        },
        'SetB': {
            'TraitB': 2,
            'TraitD': 3
        }
    }
    assert to_freq_matrix(freqs) == [[2, 1, 3, 0],
                                     [0, 2, 0, 3]]
