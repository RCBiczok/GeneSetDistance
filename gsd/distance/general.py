from typing import List, Dict, Set, Callable

from Cython.Utils import OrderedSet
from scipy.spatial.distance import pdist
import numpy as np
from sklearn.metrics import cohen_kappa_score

from gsd import flat_list, quote
from gsd.distance import DistanceMetric, calc_pairwise_distances
from gsd.gene_sets import GeneSet


def to_gene_trait_map(gene_sets: List[GeneSet]) -> Dict[str, Set[str]]:
    return {gene_set.general_info.name: set(
        flat_list([traits for gene_name, traits in gene_set.gwas_gene_traigs.gene_traits.items()]))
        for gene_set in gene_sets}


def to_gene_trait_freq(gene_sets: List[GeneSet]) -> Dict[str, Dict[str, int]]:
    def extrect_freq(gene_set: GeneSet):
        traits = flat_list([traits for gene_name, traits in gene_set.gwas_gene_traigs.gene_traits.items()])
        trait_freq = {trait: traits.count(trait) for trait in set(traits)}
        return trait_freq

    return {gene_set.general_info.name: extrect_freq(gene_set) for gene_set in gene_sets}


def to_gene_id_map(gene_sets: List[GeneSet]) -> Dict[str, Set[int]]:
    return {gene_set.general_info.name: gene_set.general_info.entrez_gene_ids for gene_set in gene_sets}


def to_binary_matrix(id_map: Dict[str, Set]) -> List[List]:
    all_ids = set(flat_list([id_set for gene_set_name, id_set in id_map.items()]))

    return [[ref_gene in id_set for ref_gene in all_ids] for gene_set_name, id_set in id_map.items()]


def to_freq_matrix(id_map: Dict[str, Dict[str, int]]) -> List[List]:
    all_keys = OrderedSet(flat_list([freqs.keys() for gene_set_name, freqs in id_map.items()]))

    ret = [[freqs[ref_key] if ref_key in freqs else 0 for ref_key in all_keys]
           for gene_set_name, freqs in id_map.items()]

    return ret


def kappa_distance(data: np.array) -> np.array:
    return calc_pairwise_distances(np.array(data), lambda a, b: 1 - cohen_kappa_score(a, b))


def overlap_coefficient(list_a: List[bool], list_b: List[bool]) -> float:
    intersection = float(sum([a * b for a, b in zip(list_a, list_b)]))
    return intersection / min(sum(list_a), sum(list_b))


def overlap_distance(data: np.array) -> np.array:
    return calc_pairwise_distances(np.array(data), lambda a, b: 1 - overlap_coefficient(a, b))


class MatrixBasedDistanceMetric(DistanceMetric):
    def __init__(self,
                 name: str,
                 extractor: Callable[[List[GeneSet]], List[List]],
                 dist_fun: Callable[[np.ndarray], np.ndarray]):
        self.name = name
        self.extractor = extractor
        self.dist_fun = dist_fun

    @property
    def display_name(self) -> str:
        return self.name

    def calc(self, gene_sets: List[GeneSet]) -> np.ndarray:
        return self.dist_fun(np.array(self.extractor(gene_sets)))


_GENERAL_DISTS = [
    MatrixBasedDistanceMetric("Minkowski distance (p=1) over genes",
                              lambda x: to_binary_matrix(to_gene_id_map(x)),
                              lambda y: pdist(y, 'minkowski', 1)),
    MatrixBasedDistanceMetric("Minkowski distance (p=2) over genes",
                              lambda x: to_binary_matrix(to_gene_id_map(x)),
                              lambda y: pdist(y, 'minkowski', 2)),
    MatrixBasedDistanceMetric("Jaccard distance over genes",
                              lambda x: to_binary_matrix(to_gene_id_map(x)),
                              lambda y: pdist(y, 'jaccard')),
    MatrixBasedDistanceMetric("Kappa distance over genes",
                              lambda x: to_binary_matrix(to_gene_id_map(x)),
                              kappa_distance),
    MatrixBasedDistanceMetric("Overlap distance over genes",
                              lambda x: to_binary_matrix(to_gene_id_map(x)),
                              overlap_distance),

    MatrixBasedDistanceMetric("Minkowski distance (p=1) over gene traits",
                              lambda x: to_binary_matrix(to_gene_trait_map(x)),
                              lambda y: pdist(y, 'minkowski', 1)),
    MatrixBasedDistanceMetric("Minkowski distance (p=2) over gene traits",
                              lambda x: to_binary_matrix(to_gene_trait_map(x)),
                              lambda y: pdist(y, 'minkowski', 2)),
    MatrixBasedDistanceMetric("Jaccard distance over gene traits",
                              lambda x: to_binary_matrix(to_gene_trait_map(x)),
                              lambda y: pdist(y, 'jaccard')),
    MatrixBasedDistanceMetric("Kappa distance over gene traits",
                              lambda x: to_binary_matrix(to_gene_trait_map(x)),
                              kappa_distance),
    MatrixBasedDistanceMetric("Overlap distance over gene traits",
                              lambda x: to_binary_matrix(to_gene_trait_map(x)),
                              overlap_distance),

    MatrixBasedDistanceMetric("Minkowski distance (p=1) over gene trait frequency",
                              lambda x: to_freq_matrix(to_gene_trait_freq(x)),
                              lambda y: pdist(y, 'minkowski', 1)),
    MatrixBasedDistanceMetric("Minkowski distance (p=2) over gene trait frequency",
                              lambda x: to_freq_matrix(to_gene_trait_freq(x)),
                              lambda y: pdist(y, 'minkowski', 2)),
    MatrixBasedDistanceMetric("Cosine distance over gene trait frequency",
                              lambda x: to_freq_matrix(to_gene_trait_freq(x)),
                              lambda y: pdist(y, 'cosine')),
]

GENERAL_DISTS = {quote(dist.display_name): dist for dist in _GENERAL_DISTS}
