from typing import Set
from ..hin import HIN
from .hash import HashMotif
from .count_dict import CountDict


def count_triangle_based_4_node_motifs(hin: HIN,
                                       edge_id: int,
                                       Si: Set[int],
                                       Sj: Set[int],
                                       Tij: Set[int],
                                       counts: CountDict,
                                       hf: HashMotif,
                                       comb: bool):

    """
    Derive triangle-based 4-node motifs.

    :param hin: HIN
        the underlying graph
    :param edge_id: int
        edge ID of the current edge between nodes i and j
    :param Si: Set[int]
        set of node IDs that are connected to i (and not j)
    :param Sj: Set[int]
        set of node IDs that are connected to j (and not i)
    :param Tij: Set[int]
        set of node IDs that are connected to i and j
    :param counts : CountDict
        maintain local and global motif counts, as well as local orbit counts
    :param hf: HashMotif
        class that can en- and decode motifs to hash strings
    :param comb : bool
        Flag to signal the use of combinatorial relationships for efficiency
    """

    i, j = hin.edges[edge_id]
    t_i, t_j = hin.nodes[i].type, hin.nodes[j].type

    for k in Tij:
        t_k = hin.nodes[k].type
        for r in hin.neighbors[k]:
            if r != i and r != j:
                t_r = hin.nodes[r].type
                if r in Tij and r < k:  # 4-clique
                    mh, oh = hf.hash_motif(12, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
                elif r in Si or r in Sj:    # chordal-cycle (edge orbit)
                    mh, oh = hf.hash_motif(10, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
                elif r not in Si and r not in Sj and r not in Tij:  # tailed-triangle (center orbit)
                    mh, oh = hf.hash_motif(8, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)

        if not comb:    # these will be derived with combinatorial relationships instead

            for r in Tij:
                if r < k and not hin.connected(r, k):  # chordal cycle (center orbit)
                    t_r = hin.nodes[r].type
                    mh, oh = hf.hash_motif(11, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)

            for r in Si:
                if r != j and not hin.connected(r, k):  # tailed triangle (tri-edge orbit)
                    t_r = hin.nodes[r].type
                    mh, oh = hf.hash_motif(9, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)

            for r in Sj:
                if r != i and not hin.connected(r, k):  # tailed triangle (tri-edge orbit)
                    t_r = hin.nodes[r].type
                    mh, oh = hf.hash_motif(9, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
