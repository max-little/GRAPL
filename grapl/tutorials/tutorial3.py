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

# Load an ADMG description from a GRAPL file and display it
G = grapl_obj.readgrapl(open('grapl/graphs/front_door.grapl', 'r').read())
G.display()

print('\nCheck if the ADMG is acyclic:')
print(G.isayclic()) # True

print('\nCompute all ADMG districts:')
print(G.districts()) # [{'M'}, {'X', 'Y'}]

print('\nInterventional (cause-effect) distribution of X on Y:')
id_str, expr, isident = algs.idfixing(G, {'X'}, {'Y'})
if isident:
    print(id_str) # p_{X}(Y)=\sum_{M,X'}[p(Y|M,X')p(M|X)p(X')]
else:
    print('Interventional distribution not identifiable')
