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

# Load an ADMG description from a GRAPL file
G = grapl_obj.readgrapl(open('grapl/graphs/richardson_2017.grapl', 'r').read())
G.display()

print('\nCheck if the ADMG is acyclic:')
print(G.isayclic()) # True

print('\nTopological sort of observed nodes:')
print(G.topsort()) # ['X1', 'X2', 'X3', 'X4']

print('\nCompute all ADMG districts:')
print(G.districts()) # [{'X3'}, {'X2', 'X4'}, {'X1'}]

print('\nInterventional (cause-effect) distribution of X2 on X4:')
id_str, expr, isident = algs.idfixing(G, {'X2'}, {'X4'})
if isident:
    print(id_str) # p_{X2}(X4)=\sum_{X3,X2',X1',X1}[p(X4|X3,X2',X1')p(X3|X2,X1)p(X2',X1')p(X1)]
else:
    print('Interventional distribution not identifiable')
