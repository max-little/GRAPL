# GRAPL
GRAPL: A computational library for nonparametric structural causal modelling, analysis and inference

(CC BY-SA 4.0) 2022. If you use this code, please cite: M.A. Little, R. Badawy, 2019, "Causal bootstrapping", arXiv:1910.09648

## Installation and getting started
- Download and unzip the package to some local directory.
- The files tutorial1.py, tutorial2.py, tutorial3.py and tutorial4.py contain demonstrations, graduated by complexity, for how to use the package.

## Release notes, v1.3
- Tutorial introduction to the library
- Comprehensive unit tests
- Factorized, marginalized joint distributions for DAGs and ADMGs
- Marginalized truncated factorization for DAGs
- List all local Markov conditional independence relations for DAGs
- Topological sorting and DAG/ADMG acyclicity tests
- All directed relationships e.g. child, parent, ancestor, descendent now operate on sets of nodes, giving e.g. inclusive parents for a given set of nodes
- All factorization/identification algorithms now take a set of effect nodes if specified, otherwise the effect is all non-cause nodes
- All factorization/identification algorithms now output a (potentially) simplified  Expr object for further processing, as well as a corresponding conditional distribution string
