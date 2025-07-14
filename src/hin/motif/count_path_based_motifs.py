from typing import Set
from ..hin import HIN
from .hash import HashMotif
from .count_dict import CountDict


def count_path_based_4_node_motifs(hin: HIN,
                                   edge_id: int,
                                   Si: Set[int],
                                   Sj: Set[int],
                                   counts: CountDict,
                                   hf: HashMotif,
                                   comb: bool):

    """
    Derive path-based 4-node motifs.

    :param hin: HIN
        the underlying graph
    :param edge_id: int
        edge ID of the current edge between nodes i and j
    :param Si: Set[int]
        set of node IDs that are connected to i (and not j)
    :param Sj: Set[int]
        set of node IDs that are connected to j (and not i)
    :param counts : CountDict
        maintain local and global motif counts, as well as local orbit counts
    :param hf: HashMotif
        class that can en- and decode motifs to hash strings
    :param comb : bool
        Flag to signal the use of combinatorial relationships for efficiency
    """

    i, j = hin.edges[edge_id]
    t_i, t_j = hin.nodes[i].type, hin.nodes[j].type

    for k in Si:

        t_k = hin.nodes[k].type
        for r in hin.neighbors[k]:
            if r != i and r != j:
                t_r = hin.nodes[r].type
                if not hin.connected(r, i) and not hin.connected(r, j):  # 4-path (edge orbit)
                    mh, oh = hf.hash_motif(3, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
                elif r in Si and r < k:  # tailed-triangle (tail orbit)
                    mh, oh = hf.hash_motif(7, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)

        if not comb:    # these will be derived with combinatorial relationships instead

            for r in Si:
                if r != j and r < k and not hin.connected(r, k):  # 4-star
                    t_r = hin.nodes[r].type
                    mh, oh = hf.hash_motif(5, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)

            for r in Sj:    # 4-path (center orbit)
                if r != i and not hin.connected(r, k):
                    t_r = hin.nodes[r].type
                    mh, oh = hf.hash_motif(4, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)

    for k in Sj:

        t_k = hin.nodes[k].type
        for r in hin.neighbors[k]:  # all neighbours of k
            if r != i and r != j:
                t_r = hin.nodes[r].type
                if not hin.connected(r, i) and not hin.connected(r, j):  # 4-path (edge orbit)
                    mh, oh = hf.hash_motif(3, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
                elif r in Sj and r < k:   # tailed-triangle (tail orbit)
                    mh, oh = hf.hash_motif(7, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
                elif r in Si:  # 4-cycle
                    mh, oh = hf.hash_motif(6, t_i, t_j, t_k, t_r)

                    counts.update(edge_id, mh, oh)

        if not comb:    # these will be derived with combinatorial relationships instead

            for r in Sj:
                if r != i and r < k and not hin.connected(r, k):  # 4-star
                    t_r = hin.nodes[r].type
                    mh, oh = hf.hash_motif(5, t_i, t_j, t_k, t_r)
                    counts.update(edge_id, mh, oh)
