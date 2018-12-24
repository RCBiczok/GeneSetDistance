from typing import Set


class GeneSet:
    def __init__(self,
                 name: str,
                 external_id: str,
                 external_source: str,
                 summary: str,
                 calculated: bool,
                 entregene_ids: Set[float],
                 gene_symbols: Set[str]):
        self.name = name
        self.external_id = external_id
        self.external_source = external_source
        self.summary = summary
        self.calculated = calculated
        self.entregene_ids = entregene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneSet(name='%s', n_external_id='%s')>" % (self.name, len(self.external_id))
