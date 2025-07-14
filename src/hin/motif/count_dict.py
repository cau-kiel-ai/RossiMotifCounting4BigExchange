from __future__ import annotations
from typing import Dict
import os
import json


class CountDict:

    def __init__(self):
        """
        Initialize the count dictionary to maintain orbit counts (always local) and local and global motif counts.
        Note that the global motif counts need to be corrected once (!) after finishing the counting.

        """
        self.orbit_count: Dict[int, Dict[str, int]] = {}
        self.local_count: Dict[int, Dict[str, int]] = {}
        self.global_count: Dict[str, int] = {}

    def load_from_json(self, directory: str):
        """ Initialize self with counts stored at specified directory.

        Directory must include the 3 files 'orbit_counts.json', 'local_counts.json', and 'global_counts.json'.

        :param directory: str
            path to the directory where files are stored
        """
        tmp = json.load(open(os.path.join(directory, "orbit_counts.json"), 'r'))
        for e in tmp:
            self.orbit_count[int(e)] = {}
            for orbit in tmp[e]:
                self.orbit_count[int(e)][orbit] = int(tmp[e][orbit])

        tmp = json.load(open(os.path.join(directory, "local_counts.json"), 'r'))
        for e in tmp:
            self.local_count[int(e)] = {}
            for motif in tmp[e]:
                self.local_count[int(e)][motif] = int(tmp[e][motif])

        tmp = json.load(open(os.path.join(directory, "global_counts.json"), 'r'))
        for motif in tmp:
            self.global_count[motif] = int(tmp[motif])

        del tmp

    def dump_to_json(self, directory: str):
        """ Dump the resp. counts into json files at the specified directory.

        Directory will then contain the 3 files 'orbit_counts.json', 'local_counts.json', and 'global_counts.json'.

        :param directory: str
            path to the directory where files will be stored
        """
        json.dump(self.orbit_count, open(os.path.join(directory, 'orbit_counts.json'), 'w'))
        json.dump(self.local_count, open(os.path.join(directory, 'local_counts.json'), 'w'))
        json.dump(self.global_count, open(os.path.join(directory, 'global_counts.json'), 'w'))

    def update(self, edge_id: int, motif_hash: str = None, orbit_hash: str = None, count: int = 1):
        """ Update orbit count, and local/global motif counts.

        :param edge_id: int
            id of the edge for which orbit count and local motif count are updated
        :param motif_hash: str (optional)
            hash value that encodes the motif
        :param orbit_hash: str (optional)
            hash value that encodes the orbit
        :param count: int (Default: 1)
            value by which the respective counts are increased
        """

        if orbit_hash is not None:
            if orbit_hash not in self.orbit_count[edge_id]:
                self.orbit_count[edge_id][orbit_hash] = count
            else:
                self.orbit_count[edge_id][orbit_hash] += count

        if motif_hash is not None:
            if motif_hash not in self.local_count[edge_id]:
                self.local_count[edge_id][motif_hash] = count
            else:
                self.local_count[edge_id][motif_hash] += count

            if motif_hash not in self.global_count:
                self.global_count[motif_hash] = count
            else:
                self.global_count[motif_hash] += count

    def correct_global_counts(self):
        """ Correct the global motif count since each motif is counted once for each edge in the motif. """

        for motif_hash in self.global_count:
            m_id = int(motif_hash[0:2])
            if m_id == 1:   # 3-star has 2 edges
                self.global_count[motif_hash] //= 2
            elif m_id in [2, 3, 4]:   # 3-clique, 4-path, 4-star have 3 edges
                self.global_count[motif_hash] //= 3
            elif m_id in [5, 6]:    # 4-cycle, tailed triangle have 4 edges
                self.global_count[motif_hash] //= 4
            elif m_id == 7:     # chordal cycle has 5 edges
                self.global_count[motif_hash] //= 5
            elif m_id == 8:     # 4-clique has 6 edges
                self.global_count[motif_hash] //= 6

    def get_total_count(self, edge_id: int = None) -> int:
        """
        Return the total number of motifs (local or global).

        :param edge_id: int (optional)
            edge id for which to return the total count. If None, return global total
        :return: int
            total count of motifs (local or global)
        """
        total = 0
        if edge_id is None:
            for h in self.global_count:
                total += self.global_count[h]
        else:
            for h in self.local_count[edge_id]:
                total += self.local_count[edge_id][h]
        return total

    def derive_untyped_dict(self) -> CountDict:
        """ Return an untyped version of the CountDict. """

        newDict: CountDict = CountDict()

        for edge in self.orbit_count:
            newDict.orbit_count[edge] = {}
            for key in self.orbit_count[edge]:
                new_key = key[0:2]  # we have 12 potential undirected orbits
                if new_key not in newDict.orbit_count[edge]:
                    newDict.orbit_count[edge][new_key] = 0
                newDict.orbit_count[edge][new_key] += self.orbit_count[edge][key]

        for edge in self.local_count:
            newDict.local_count[edge] = {}
            for key in self.local_count[edge]:
                new_key = key[0:2]
                if new_key not in newDict.local_count[edge]:
                    newDict.local_count[edge][new_key] = 0
                newDict.local_count[edge][new_key] += self.local_count[edge][key]

        for key in self.global_count:
            new_key = key[0:2]
            if new_key not in newDict.global_count:
                newDict.global_count[new_key] = 0
            newDict.global_count[new_key] += self.global_count[key]

        return newDict
