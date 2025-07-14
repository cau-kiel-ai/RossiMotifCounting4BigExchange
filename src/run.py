import argparse
from hin.hin import HIN
from hin.dataset_loader import load_dataset
from hin.motif.count_3_4_node_motifs import count_motifs
import cProfile
import pstats
from hin.motif.count_dict import CountDict
import os


descr = '''
    Count 3- and 4-node motifs in an HIN.\n
    Usage example on Windows:\n
    python run.py --dataset ..\\data\\DBLP4areas --output ..\\results --no_comb
    Usage example on Linux:\n
    python run.py --dataset ../data/DBLP4areas --output ../results
    '''


parser = argparse.ArgumentParser(description=descr)
parser.add_argument("-d", "--dataset",
                    help="Path to dataset folder",
                    required=True)
parser.add_argument("-o", "--output",
                    help="Path to a folder where results will be placed",
                    required=True)
parser.add_argument("--no_comb",
                    help="Turns off the use combinatorial relationships",
                    action="store_true")
args = parser.parse_args()

path_to_dataset: str = args.dataset
path_to_output: str = args.output


if not os.path.exists(path_to_dataset):
    raise FileNotFoundError(f"Dataset path does not exist: {path_to_dataset}")
if not os.path.exists(path_to_output):
    raise FileNotFoundError(f"Output path does not exist: {path_to_output}")


profiler = cProfile.Profile()

profiler.enable()
hin: HIN = load_dataset(path_to_dataset)
counts: CountDict = count_motifs(hin, comb=not args.no_comb)
profiler.disable()

stats = pstats.Stats(profiler).sort_stats('tottime')
stats.dump_stats(os.path.join(path_to_output, 'timing.pstats'))
counts.dump_to_json(path_to_output)
