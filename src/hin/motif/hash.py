from typing import Dict, Set


class HashMotif:

    def __init__(self, n_types: Set[str]):
        """ Initialize the hash function based on the number of node types in the HIN schema.

        :param n_types: List[str]
            list of the node type names
        """
        self.n_types: Dict[str, int] = {'--': -1}
        self.n_types_num: Dict[int, str] = {-1: '--'}
        for i, n_type in enumerate(n_types):
            if n_type.isdigit():    # so that the hash value matches the original type
                i = int(n_type)
            self.n_types[n_type] = i
            self.n_types_num[i] = n_type

    def decode(self, hash_str: str) -> (int, str, str, str, str):
        """ Decode hash string into a more readable format. """
        g = int(hash_str[0:2])
        i = self.n_types_num[int(hash_str[2:4])]
        j = self.n_types_num[int(hash_str[4:6])]
        k = self.n_types_num[int(hash_str[6:8])]
        if hash_str[8:10] == '--':
            r = '--'
        else:
            r = self.n_types_num[int(hash_str[8:10])]
        return g, i, j, k, r

    def hash_motif(self, g: int, i: str, j: str, k: str, r: str) -> (str, str):
        """
        Compute a hash that encodes the edge types in the motif (for max. 100 types in the schema).
        (cf. section 4.4 of 'Heterogeneous Graphlets' by Rossi et al.)

        :param g: int
            orbit ID
        :param i: str
            type of node i
        :param j: str
            type of node j
        :param k: str
            type of node k
        :param r: str
            type of node r
        :return: (str, str)
            motif hash and orbit hash that encode the node types
        """

        # derive motif ID m from orbit ID g
        if g == 1 or g == 2:
            m = g  # 3-star/3-path (1) or triangle/3-clique (2)
        elif g == 3 or g == 4:
            m = 3  # 4-path
        elif g == 5:
            m = 4  # 4-star
        elif g == 6:
            m = 5  # 4-cycle
        elif g == 7 or g == 8 or g == 9:
            m = 6  # tailed triangle
        elif g == 10 or g == 11:
            m = 7  # chordal cycle
        elif g == 12:
            m = 8  # 4-clique
        else:
            raise Exception(f"Invalid orbit ID ({g}). ID must be between 1 and 12.")

        if r == '--':
            types = sorted([self.n_types[i], self.n_types[j], self.n_types[k]])
            type_sum = types[0] * pow(10, 4) + types[1] * pow(10, 2) + types[2]
            g = g * pow(10, 6)
            m = m * pow(10, 6)
            orbit_hash = f"{g + type_sum:08}--"
            motif_hash = f"{m + type_sum:08}--"
        else:
            types = sorted([self.n_types[i], self.n_types[j],
                            self.n_types[k], self.n_types[r]])
            type_sum = types[0] * pow(10, 6) + types[1] * pow(10, 4) + types[2] * pow(10, 2) + types[3]
            g = g * pow(10, 8)
            m = m * pow(10, 8)
            orbit_hash = f"{g + type_sum:010}"
            motif_hash = f"{m + type_sum:010}"

        return motif_hash, orbit_hash
