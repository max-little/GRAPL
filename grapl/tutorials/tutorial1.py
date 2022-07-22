# -*- coding: utf-8 -*-
"""
Tutorial illustrating the use of the GRAPL library, for representing, analyzing and processing
acyclic directed mixed graphs for structural causal modelling.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.algorithms as algs
import grapl.dsl as dsl

# Create a GRAPL DSL parser
grapl_obj = dsl.GraplDSL()

# Create a DAG from a GRAPL string
dag_grapl = ' "Simple back-door graph"; \
    C; X; Y; \
    C -> X; \
    C -> Y; \
    X -> Y; '
G = grapl_obj.readgrapl(dag_grapl)
G.display()

print('\nCheck if the DAG is acyclic:')
print(G.isdag()) # True

print('\nFactorized joint distribution:')
fac_str, expr, isdag = algs.dagfactor(G, simplify=False)
print(fac_str) # p(X,C,Y)=[p(Y|X,C)p(X|C)p(C)]

print('\nInterventional (cause-effect) distribution of X on Y:')
id_str, expr, isdag = algs.truncfactor(G, {'X'}, {'Y'})
print(id_str) # p_{X}(Y)=\sum_{C}[p(Y|X,C)p(C)]
