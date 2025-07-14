# Counting Heterogeneous Graphlets

This implementation counts 3- and 4-node motifs (also referred to as graphlets)
in heterogeneous information networks (HINs).
It was created by Mattis thor Straten and Christian Beth and is based on the
paper _"Heterogeneous Graphlets"_ by Rossi et al., which can be found
[here](https://dl.acm.org/doi/abs/10.1145/3418773?casa_token=W4FnTuzKs5IAAAAA:VC3vCOfR6vS-3LB7XkfbEzgZBAaaB_y9eek_jdcAmsixqrI6OP2y0ts8T8gWTqkBD6lNLaK0Q3Vihg).

When using this implementation, credit should be given to the paper's
original authors:
````
@article{rossi2020heterogeneous,
  title={Heterogeneous graphlets},
  author={Rossi, Ryan A and Ahmed, Nesreen K and Carranza, Aldo and Arbour, David and Rao, Anup and Kim, Sungchul and Koh, Eunyee},
  journal={ACM Transactions on Knowledge Discovery from Data (TKDD)},
  volume={15},
  number={1},
  pages={1--43},
  year={2020},
  publisher={ACM New York, NY, USA}
}
````


### Usage Guide

In order to run experiments, run the `run.py` file as follows:
````
python run.py --dataset <path to dataset> --output <path to output directory> [--no_comb]
````
or
````
python run.py --d <path to dataset> --o <path to output directory> [--no_comb]
````
where `--no_comb` signals that no combinatorial relationships should be used.
This would result in the code running significantly slower though.

A concrete running example (if run out of the box) could look as follows:
````
python run.py --dataset ..\data\BigExchange --output ..\results\BigExchange
````
which runs the counting on the Big Exchange dataset on a Windows system, or
````
python run.py --dataset ../data/BigExchange --output ../results/BigExchange
````
which runs the counting on the DBLP dataset on a Linux system.

### Interpretation of the Results

After a successful run the algorithm creates the following files in the
specified output directory:
- `orbit_counts.json`:
This contains the per-edge _orbit_ counts in a nested dictionary. The outer
level is indexed by the edge ID, the inner layer by the orbit hash value.
The innermost level contains the respective orbit counts.
- `local_counts.json`:
This contains the per-edge _motif_ counts in a nested dictionary. The outer
level is indexed by the edge ID, the inner layer by the motif hash value.
The innermost level contains the respective motif counts.
- `global_counts.json`:
This contains the global motif counts in a nested dictionary. The outer
level is indexed by motif hash value. The innermost level contains the
respective motif counts
- `timing.pstats`:
This contains a breakdown of how much time has been consumed by individual
function calls and how frequently functions have been called. It can conveniently 
be visualized with [snakeviz](https://jiffyclub.github.io/snakeviz/) in the
browser by running the following in the command line: `snakeviz timing.pstats`

Each orbit and Motif are encoded by an invertible hash function,
which is modelled in the `HashMotif` class.
Each node type is assigned an integer value, which can be accessed in the dictionary
`HashMotif.n_types`.
The first to digits of the hash string encode the motif (resp. orbit) ID as it was
used in the work by Rossi et al.
This ID is followed up to 4 node types encoded by their respective integer value
that take up two digits each.
For example, the motif hash `'02001005--'` encodes the motif `02`  (a triangle) with
node types `00`, `10`, and `05`.
`'--''`is used in 3-node motifs (wedge and triangle) to replace the missing 4th node
type.
Note that the node type integers are sorted in ascending order, which provides a canonical
encoding, but also discards the information of the topology of the types within the motif.

The motif IDs represent the following motifs:
- `01`: Wedge (sometimes also 2-star or 2-path)
- `02`: Triangle (sometimes also 3-clique)
- `03`: 3-path (in Rossi et al.'s paper referred to as 4-path!)
- `04`: 3-star (in Rossi et al.'s paper referred to as 4-star!)
- `05`: 4-cycle
- `06`: Tailed triangle
- `07`: Diamond (sometimes also chordal cycle)
- `08`: 4-clique

The which orbit IDs correspond to which motif IDs is described in `hin.motif.hash.py`.

For questions please contact [Christian](mailto:cbe@informatik.uni-kiel.de).