from typing import Set
from ..hin import HIN
from .hash import HashMotif
from .count_dict import CountDict
from math import comb


def derive_comb_counts(hin: HIN,
                       edge_id: int,
                       Si: Set[int],
                       Sj: Set[int],
                       Tij: Set[int],
                       counts: CountDict,
                       hf: HashMotif):
    """
    Derive remaining motif counts from combinatorial relationships (for g4, g5, g9, g11).

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
    """

    i, j = hin.edges[edge_id]
    t_i, t_j = hin.nodes[i].type, hin.nodes[j].type

    L: Set[str] = hin.node_types

    for t1 in L:
        for t2 in L:

            if t1 > t2:
                continue    # to avoid duplicate enumerations

            Si1: Set[int] = set()
            Si2: Set[int] = set()
            for v in Si:
                if v == j:
                    continue
                if hin.nodes[v].type == t1:
                    Si1.add(v)
                elif hin.nodes[v].type == t2:
                    Si2.add(v)

            Sj1: Set[int] = set()
            Sj2: Set[int] = set()
            for v in Sj:
                if v == i:
                    continue
                if hin.nodes[v].type == t1:
                    Sj1.add(v)
                elif hin.nodes[v].type == t2:
                    Sj2.add(v)

            Tij1: Set[int] = set()
            Tij2: Set[int] = set()
            for v in Tij:
                if hin.nodes[v].type == t1:
                    Tij1.add(v)
                elif hin.nodes[v].type == t2:
                    Tij2.add(v)

            # derive g_4 (4-path center orbit)
            _, h6 = hf.hash_motif(6, t_i, t_j, t1, t2)
            if h6 in counts.orbit_count[edge_id]:
                g_6 = counts.orbit_count[edge_id][h6]
            else:
                g_6 = 0
            if t1 == t2:  # cf. eq. 19 of "Heterogeneous Graphlets", Rossi et al. (TKDD'2020)
                count_4 = len(Si1) * len(Sj1) - g_6
            else:
                count_4 = len(Si1) * len(Sj2) + len(Si2) * len(Sj1) - g_6
            if count_4 > 0:
                mh, oh = hf.hash_motif(4, t_i, t_j, t1, t2)
                counts.update(edge_id, mh, oh, count=count_4)
            
            # derive g_5 (4-star)
            _, h7 = hf.hash_motif(7, t_i, t_j, t1, t2)
            if h7 in counts.orbit_count[edge_id]:
                g_7 = counts.orbit_count[edge_id][h7]
            else:
                g_7 = 0
            if t1 == t2:    # cf. eq. 23 of "Heterogeneous Graphlets", Rossi et al. (TKDD'2020)
                count_5 = int(comb(len(Si1), 2)) + int(comb(len(Sj1), 2)) - g_7
            else:
                count_5 = len(Si1) * len(Si2) + len(Sj1) * len(Sj2) - g_7
            if count_5 > 0:
                mh, oh = hf.hash_motif(5, t_i, t_j, t1, t2)
                counts.update(edge_id, mh, oh, count=count_5)

            # derive g_9 (tailed triangle tri-edge orbit)
            _, h10 = hf.hash_motif(10, t_i, t_j, t1, t2)
            if h10 in counts.orbit_count[edge_id]:
                g_10 = counts.orbit_count[edge_id][h10]
            else:
                g_10 = 0
            if t1 == t2:    # cf. eq. 26 of "Heterogeneous Graphlets", Rossi et al. (TKDD'2020)
                count_9 = len(Tij1) * (len(Si1) + len(Sj1)) - g_10
            else:
                count_9 = len(Tij1) * (len(Si2) + len(Sj2)) + len(Tij2) * (len(Si1) + len(Sj1)) - g_10
            if count_9 > 0:
                mh, oh = hf.hash_motif(9, t_i, t_j, t1, t2)
                counts.update(edge_id, mh, oh, count=count_9)

            # derive g_11 (chordal cycle center orbit)
            _, h12 = hf.hash_motif(12, t_i, t_j, t1, t2)
            if h12 in counts.orbit_count[edge_id]:
                g_12 = counts.orbit_count[edge_id][h12]
            else:
                g_12 = 0
            if t1 == t2:    # # cf. eq. 30 of "Heterogeneous Graphlets", Rossi et al. (TKDD'2020)
                count_11 = int(comb(len(Tij1), 2)) - g_12
            else:
                count_11 = len(Tij1) * len(Tij2) - g_12
            if count_11 > 0:
                mh, oh = hf.hash_motif(11, t_i, t_j, t1, t2)
                counts.update(edge_id, mh, oh, count=count_11)
