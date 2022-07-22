# -*- coding: utf-8 -*-
"""
Unit tests for the GRAPL library, for representing, analyzing and processing
acyclic directed mixed graphs for structural causal modelling.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.algorithms as algs
import grapl.dsl as dsl

def run_all():
    test_cyclic()
    test_topsort()
    test_tianfactor()
    test_fixing()
    test_identify1()
    test_identify2()
    test_non_identify1()
    test_non_identify2()
    test_complex_dag()
    test_markov()

    print('All unit tests passed')


def test_cyclic():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/cyclic_graph.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert not G.isdag()
    assert not G.isayclic()

def test_topsort():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/topsort_graph.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    node_list = G.topsort()
    for node in node_list:
        for parent in G.pa({node}):
            assert node_list.index(node) > node_list.index(parent)

def test_tianfactor():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/bareinboim_2020.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert G.isayclic()
    assert all([(D in G.districts()) for D in [{'A'},{'F'},{'B','D'},{'E','C'}]])
    assert all([(v in G.an({'D'})) for v in {'A','B','C','D'}])
    assert all([(v in G.de({'B','F'})) for v in {'B','D','E','F'}])

def test_fixing():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/richardson_2017.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixing(G, {'X2'}, {'X4'})
    assert not G.isdag()
    assert G.isayclic()
    assert isident

def test_identify1():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/shpitser_thesis1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixing(G, {'X'}, {'Y1','Y2'})
    assert not G.isdag()
    assert G.isayclic()
    assert isident
    assert all([(term in expr.num) for term in [{'X', "W1'", 'Y1'},{"W1'"},{'Y2'}]])
    assert {'X', "W1'"} in expr.den
    assert "W1'" in expr.mrg

def test_identify2():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/shpitser_thesis2.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixing(G, {'X'}, {'Y'})
    assert not G.isdag()
    assert G.isayclic()
    assert isident

def test_non_identify1():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/non_identify1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixing(G, {'X'}, {'Y'})
    assert not G.isdag()
    assert G.isayclic()
    assert not isident

def test_non_identify2():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/non_identify2.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixing(G, {'X'}, {'Y'})
    assert not G.isdag()
    assert G.isayclic()
    assert not isident

def test_complex_dag():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/complex_dag.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_str, expr, isdag = algs.truncfactor(G, {'X1'}, {'X3','X4','Y1'})
    assert all([(term in expr.num) for term in [{'Y1','C','X1'},{'X3'},{'X4'},{'C'}]])
    assert {'C','X1'} in expr.den
    assert 'C' in expr.mrg

def test_markov():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = 'grapl/graphs/markov1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert G.nd({'D'}) == {'S','I'}
    assert G.nd({'G'}) == {'D','S','I'}
    assert G.nd({'L'}) == {'D','S','I','G'}
    assert G.nd({'S'}) == {'D','I','G','L'}
    assert G.nd({'I'}) == {'D'}
