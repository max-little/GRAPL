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
dag_grapl = ' "Symptoms"; \
    Sinus; Headache; Nose; Flu; Allergy; \
    Sinus -> Nose; \
    Flu -> Sinus; \
    Sinus -> Headache; \
    Allergy -> Sinus;'
G = grapl_obj.readgrapl(dag_grapl)
G.display()

print('\nCheck if the DAG acyclic:')
print(G.isdag()) # True

print('\nTopological sort of DAG nodes:')
print(G.topsort()) # ['Flu', 'Allergy', 'Sinus', 'Nose', 'Headache']

print('\nAll ancestors of Headache:')
print(G.an({'Headache'})) # {'Allergy', 'Headache', 'Sinus', 'Flu'}

print('\nAll DAG local Markov conditional independence relations:')
ci, isdag = algs.localmarkov(G) # {'Flu⊥Allergy', 'Nose⊥Headache,Allergy,Flu|Sinus', 'Allergy⊥Flu', 'Headache⊥Nose,Allergy,Flu|Sinus'}
print(ci)

print('\nFactorized joint distribution over Headache,Allergy,Flu:')
factor_str, expr, isdag = algs.dagfactor(G, {'Headache','Allergy','Flu'})
print(factor_str) # p(Flu,Allergy,Headache)=\sum_{Sinus}[p(Sinus|Flu,Allergy)p(Headache|Sinus)p(Allergy)p(Flu)]

print('\nInterventional (cause-effect) distribution of Sinus on Headache,Allergy,Flu:')
id_str, expr, isdag = algs.truncfactor(G, {'Sinus'}, {'Headache','Allergy','Flu'})
print(id_str) # p_{Sinus}(Flu,Allergy,Headache)=[p(Headache|Sinus)p(Allergy)p(Flu)]
