import pytest
from anytree import RenderTree

import gsd.reactome

human_tax_id = 9606
reactome_id = "R-HSA-416550"


@pytest.mark.skip(reason="just takes long")
def test_download_reactome():
    node, gene_sets = gsd.reactome.download(human_tax_id, reactome_id)
    print(RenderTree(node))
    assert node is not None
    assert list is not None
