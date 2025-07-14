from __future__ import annotations
from typing import List, Set


HINEdge = (int, int)


class HINNode:
    """
    A node of an HIN. Consisting of the ID, the value (e.g. the name of the author) and the node type.
    """

    def __init__(self, n_id: int, n_type: str):
        """ Initializes a node in an HIN.

        :param n_id: int
            node id in the HIN
        :param n_type: str
            node type name
        """
        self.n_id: int = n_id
        self.type: str = n_type


class HIN:

    def __init__(self, nodes: List[HINNode], edges: List[HINEdge]):
        """
        Initializes an HIN from a List of HINNodes and a list of HINEdges.

        :param nodes : List[HINNode]
            list of nodes
        :param edges : List[HINEdge]
            list of links

        Attributes
        ----------
        nodes : List[HINNode]
            store HINNode object, where list index corresponds to node ID
        edges : List[(int, int)]
            stores the node IDs of the connected nodes for an edge ID, which corresponds to list index
        neighbors: List[Set[int]]
            stores the neighboring node IDs, where the list index corresponds to node ID
        """

        self.nodes: List[HINNode] = []
        self.neighbors: List[Set[int]] = []
        self.edges: List[(int, int)] = []
        self.node_types: Set[str] = set()

        for v in nodes:
            self.nodes.append(v)
            self.neighbors.append(set())
            self.node_types.add(v.type)

        for i, j in edges:
            self.edges.append((i, j))
            self.neighbors[i].add(j)
            self.neighbors[j].add(i)

    def connected(self, i: int, j: int) -> bool:
        """ Return true if node i and node j are connected by an edge in the network. """
        if j in self.neighbors[i]:
            return True
        else:
            return False
