from typing import Set
from ..hin import HIN
from .hash import HashMotif
from .count_dict import CountDict
from .count_path_based_motifs import count_path_based_4_node_motifs
from .count_triangle_based_motifs import count_triangle_based_4_node_motifs
from .comb_relationships import derive_comb_counts


def count_per_edge(hin: HIN,
                   edge_id: int,
                   counts: CountDict,
                   hf: HashMotif,
                   comb: bool = True):
    """
    Count all 3- and 4-node motifs that a given edge in its respective HIN participates in
    (without combinatorial relationships).

    The motif counts are managed in two dictionaries M and X. M assigns each edge ID a set of all motifs the edge
    participates in (motifs are encoded by a string hash) and X assigns each edge ID a dictionary of counts for each
    motif (encoded by string hash).

    :param hin : HIN
        the underlying graph
    :param edge_id : int
        edge id of the current edge between nodes i and j
    :param counts : CountDict
        maintain local and global motif counts, as well as local orbit counts
    :param hf : HashMotif
        class that can en- and decode motifs to hash strings
    :param comb : bool
        Flag to signal the use of combinatorial relationships for efficiency (Default: True)
    """

    # add an empty count entry for the new edge ID
    counts.orbit_count[edge_id] = {}
    counts.local_count[edge_id] = {}

    # ID of node i, ID of node j
    i, j = hin.edges[edge_id]
    # types of nodes i and j
    t_i, t_j = hin.nodes[i].type, hin.nodes[j].type

    # set of node IDs that...
    Si: Set[int] = set()      # ... form 3-paths centered at node i
    Sj: Set[int] = set()      # ... form 3-paths centered at node j
    Tij: Set[int] = set()     # ... form triangles with nodes i and j

    for k in hin.neighbors[i]:
        if k != j:
            Si.add(k)   # k may later be moved to Tij

    for k in hin.neighbors[j]:
        t_k = hin.nodes[k].type
        if k != i:
            if k in Si:  # if k is neighbour of i
                Si.remove(k)
                Tij.add(k)
                mh, oh = hf.hash_motif(2, t_i, t_j, t_k, '--')
                counts.update(edge_id, mh, oh)
            else:  # 3-star
                Sj.add(k)
                mh, oh = hf.hash_motif(1, t_i, t_j, t_k, '--')
                counts.update(edge_id, mh, oh)

    for k in Si:
        t_k = hin.nodes[k].type
        mh, oh = hf.hash_motif(1, t_i, t_j, t_k, '--')
        counts.update(edge_id, mh, oh)

    count_path_based_4_node_motifs(hin, edge_id, Si, Sj, counts, hf, comb)

    count_triangle_based_4_node_motifs(hin, edge_id, Si, Sj, Tij, counts, hf, comb)

    if comb:    # derive remaining motif counts from combinatorial relationships (for g4, g5, g9, g11)
        derive_comb_counts(hin, edge_id, Si, Sj, Tij, counts, hf)


def count_motifs(hin, comb: bool = True) -> CountDict:
    """ Count all 3- and 4-node motifs in an HIN.

    :param hin: HIN
        the graph for which all 3- and 4-node motifs are to be counted
    :param comb: bool
        flag to signal whether to utilize combinatorial relationships (default: True)
    :return: CountDict
        a CountDict object that contains orbit counts, as well as local and global motif counts
    """

    counts = CountDict()
    hf = HashMotif(hin.node_types)
    for e_ij, _ in enumerate(hin.edges):
        count_per_edge(hin, e_ij, counts, hf, comb=comb)
    counts.correct_global_counts()
    return counts
