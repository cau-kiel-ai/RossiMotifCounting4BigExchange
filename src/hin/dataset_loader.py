import os
from .hin import HIN, HINNode


def load_dataset(path: str) -> HIN:
    curr_dir = os.getcwd()

    with open(os.path.join(curr_dir, path, 'nodes.csv'), 'r') as node_file:
        nodes = []
        lines = node_file.readlines()
        for node_id, line in enumerate(lines):
            node_type = line.strip()
            node = HINNode(node_id, node_type)
            nodes.append(node)
        del lines

    with open(os.path.join(curr_dir, path, 'edges.csv'), 'r') as edge_file:
        edges = []
        for line in edge_file.readlines():
            chunks = line.strip().split(',')
            edge = int(chunks[0]), int(chunks[2])
            edges.append(edge)

    return HIN(nodes, edges)
