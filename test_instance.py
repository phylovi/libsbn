import numpy as np
import sbn

def test_instance():
    inst = sbn.instance('charlie')
    inst.parse_file('data/four_taxon.tre')
    inst.print_status()
    assert inst.tree_count() == 5
    sbn.f(np.array([3,4]))
    # m = sbn.make_tag_to_bitset_map(inst.trees[0])
    # print(m)

