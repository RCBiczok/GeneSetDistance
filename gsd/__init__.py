class GeneSet:
    def __init__(self, name, externalId, externalSource, summary, entregene_ids, gene_symbols):
        self.name = name
        self.externalId = externalId
        self.externalSource = externalSource
        self.summary = summary
        self.entregene_ids = entregene_ids
        self.gene_symbols = gene_symbols

    def __repr__(self):
        return "<GeneSet(name='%s')>" % (self.name)
