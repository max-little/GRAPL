![GRAPL logo](https://raw.githubusercontent.com/max-little/GRAPL/main/grapl.png)
### GRAPL: A computational library for nonparametric structural causal modelling, analysis and inference

*Structural causal models* (SCMs) provide a probabilistic language for describing directed relationships between random variables. SCMs are widely used in science and engineering to capture **causal relationships** between quantitative, measured phenomena in the real world. Two SCM formalisms, *directed acyclic graphs* (DAGs) and *acyclic directed mixed graphs* (ADMGs) have been extensively studied. In these formalisms, the conditions under which causal dependence between variables occurs is well understood. Furthermore, analytical techniques have been developed which allow manipulation of the model so as to perform causal adjustment, that is, the isolation of desired *causal relationships* from the SCM. The GRAPL library brings together the most important and useful of such algorithms in one convenient software package. Using this library it is possible to represent, analyze and manipulate DAGs and ADMGs of arbitrary complexity.

## Features
- A simple text-based domain-specific language for representing both directed acyclic (DAG) and mixed, directed acyclic models (ADMGs) of interacting variables
- Derivation of factorized, marginalized nonparametric distributional models for arbitrary DAGs
- Computation of nonparametric, causal interventional distributions in arbitrarily complex models with (ADMG) and without (DAG) hidden variables, where such an interventional distribution can be computed in principle
- Various algorithms for the analysis of causal influence in DAGs/ADMGs (e.g. c-components/districts, node interventions, local Markov conditional independence relations, topological sorting)
- Latex format output distributions which can be easily dropped into documents for publication

## Installation and getting started

We currently offer seamless installation with  `pip`. Download the current distribution of the package, and run:
```
pip install .
```
in the root directory of the decompressed package.

Tutorials are included as part of the package in increasing order of complexity. These can be found in `grapl/tutorials`.
To run these directly, use:
```python
>>> import grapl.tutorials.tutorial1
```
and similarly for tutorials 2-4.

## Usage example
*Computing the front-door adjusted distribution of the causal effect of `X` on `Y` with mediator `M` (e.g. `X` &rarr; `M` &rarr; `Y`) with hidden/latent confounding back-door path `X` &mdash; `Y`.*

First, create a GRAPL DSL object, then the GRAPL description of the graph in a string, and parse the string using `grapl.dsl.readgrapl` to create the graph object `G`:

```python
>>> import grapl.algorithms as algs
>>> import grapl.dsl as dsl

>>> grapl_obj = dsl.GraplDSL()
>>> dag_grapl = ' "Front door adjustment"; \
>>>   X; Y; M; \
>>>   X -> M; \
>>>   M -> Y; \
>>>   X <-> Y; '
>>> G = grapl_obj.readgrapl(dag_grapl)
```

Next, invoke the `grapl.algorithms.idfixing` algorithm to find the interventional distribution (if it can be identified):

```python
>>> id_str, expr, isident = algs.idfixing(G, {'X'}, {'Y'})
>>> if isident:
>>>     print(id_str) # p_{X}(Y)=\sum_{M,X'}[p(Y|M,X')p(M|X)p(X')]
>>> else:
>>>     print('Interventional distribution not identifiable')
```

## Testing

```python
>>> import grapl.test.unit_tests as tests
>>> tests.run_all()
```

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
