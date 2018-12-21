class GeneSet:
    def __init__(self, name, external_id, external_source, summary, entregene_ids, gene_symbols):
        self.name = name
        self.external_id = external_id
        self.external_source = external_source
        self.summary = summary
        self.entregene_ids = entregene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneSet(name='%s', external_id='%s')>" % (self.name, self.external_id)


def flat_list(l):
    return [item for sublist in l for item in sublist]
