from typing import Set


class GeneSet:
    def __init__(self,
                 name: str,
                 external_id: str,
                 external_source: str,
                 summary: str,
                 calculated: bool,
                 entrez_gene_ids: Set[float],
                 gene_symbols: Set[str]):
        self.name = name
        self.external_id = external_id
        self.external_source = external_source
        self.summary = summary
        self.calculated = calculated
        self.entrez_gene_ids = entrez_gene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneSet(name='%s', n_entrez_gene_ids='%s')>" % (self.name, len(self.entrez_gene_ids))
