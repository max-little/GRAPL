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
    test_mconn()
    test_msep()
    test_fixall_greedy()
    test_fixall_full()

    print('All unit tests passed')


def test_cyclic():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/cyclic_graph.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert not G.isdag()
    assert not G.isayclic()

def test_topsort():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/topsort_graph.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    node_list = G.topsort()
    for node in node_list:
        for parent in G.pa({node}):
            assert node_list.index(node) > node_list.index(parent)

def test_tianfactor():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/bareinboim_2020.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert G.isayclic()
    assert all([(D in G.districts()) for D in [{'A'},{'F'},{'B','D'},{'E','C'}]])
    assert all([(v in G.an({'D'})) for v in {'A','B','C','D'}])
    assert all([(v in G.de({'B','F'})) for v in {'B','D','E','F'}])

def test_fixing():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/richardson_2017.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixing(G, {'X2'}, {'X4'})
    assert not G.isdag()
    assert G.isayclic()
    assert isident

def test_identify1():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/shpitser_thesis1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, id_eqn, isident = algs.idfixing(G, {'X'}, {'Y1','Y2'})
    assert not G.isdag()
    assert G.isayclic()
    assert isident
    assert all([(term in id_eqn.rhs.num) for term in [{'X', "W1'", 'Y1'},{"W1'"},{'Y2'}]])
    assert {'X', "W1'"} in id_eqn.rhs.den
    assert "W1'" in id_eqn.rhs.mrg

def test_identify2():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/shpitser_thesis2.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, id_eqn, isident = algs.idfixing(G, {'X'}, {'Y'})
    assert not G.isdag()
    assert G.isayclic()
    assert isident

def test_non_identify1():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/non_identify1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, id_eqn, isident = algs.idfixing(G, {'X'}, {'Y'})
    assert not G.isdag()
    assert G.isayclic()
    assert not isident

def test_non_identify2():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/non_identify2.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, id_eqn, isident = algs.idfixing(G, {'X'}, {'Y'})
    assert not G.isdag()
    assert G.isayclic()
    assert not isident

def test_complex_dag():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/complex_dag.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_str, id_eqn, isdag = algs.truncfactor(G, {'X1'}, {'X3','X4','Y1'})
    assert all([(term in id_eqn.rhs.num) for term in [{'Y1','C','X1'},{'X3'},{'X4'},{'C'}]])
    assert {'C','X1'} in id_eqn.rhs.den
    assert 'C' in id_eqn.rhs.mrg

def test_markov():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/markov1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    cis, cisstr, isdag = algs.localmarkov(G)
    for ci in cis.ciset:
        if 'I' in ci.X:
            assert ci.Y == {'D'}
            assert ci.Z == set()
        if 'S' in ci.X:
            assert ci.Y == {'G','D','L'}
            assert ci.Z == {'I'}
        if 'L' in ci.X:
            assert ci.Y == {'I','S','D'}
            assert ci.Z == {'G'}
        if 'G' in ci.X:
            assert ci.Y == {'S'}
            assert ci.Z == {'I', 'D'}
        if 'D' in ci.X:
            assert ci.Y == {'I','S'}
            assert ci.Z == set()

def test_dsep():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/markov1.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    isdsep, ci, cistr, isdag = algs.dseparate(G,{'D'},{'S'},{'G','I'})
    assert isdsep
    assert ci.X == {'D'}
    assert ci.Y == {'S'}
    assert ci.Z == {'G','I'}
    isdsep, ci, cistr, isdag = algs.dseparate(G,{'D'},{'S'},{'G'})
    assert not isdsep

def test_mconn():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/m_sep_graph.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert G.mconn('Y') == {'N', 'W', 'T', 'Z', 'V', 'Y', 'X'}
    assert G.mconn('W', {'X'}) == {'W', 'T', 'Z', 'N'}
    assert G.mconn('Y', {'V'}) == {'W', 'T', 'Z', 'Y', 'X', 'N'}
    
def test_msep():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/m_sep_graph.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    assert G.ismsep("W", "Y", {"X"})
    assert G.ismsep("Y", "V", {"X"})
    assert G.ismsep("W", "Z")
    assert not G.ismsep("W", "N", {"V"})
    assert not G.ismsep("W", "T", {"X"})
    
    
def test_fixall_greedy():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/complex_dag.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixall(G, {'X2'}, {'Y1','Y2'}, mode = "mostmrg", greedy = True)
    assert G.isdag()
    assert G.isayclic()
    assert isident
    
    id_expr, expr, isident = algs.idfixall(G, {'X2'}, {'Y1','Y2'}, mode = "random")
    assert G.isdag()
    assert G.isayclic()
    assert isident
    
    id_expr, expr, isident = algs.idfixall(G, {'X2'}, {'Y1','Y2'}, mode = "shortest", greedy = True)
    assert G.isdag()
    assert G.isayclic()
    assert isident
    
def test_fixall_full():
    grapl_obj = dsl.GraplDSL()
    grapl_filename = '/grapl/graphs/general_scenario.grapl'
    G = grapl_obj.readgrapl(open(grapl_filename, 'r').read())
    id_expr, expr, isident = algs.idfixall(G, {'Y'}, {'X'}, mode = "all", greedy = False)
    assert not G.isdag()
    assert G.isayclic()
    assert isident
    
    id_expr, expr, isident = algs.idfixall(G, {'Y'}, {'X'}, mode = "shortest", greedy = False)
    assert not G.isdag()
    assert G.isayclic()
    assert isident
    
    id_expr, expr, isident = algs.idfixall(G, {'Y'}, {'X'}, mode = "mostmrg", greedy = False)
    assert not G.isdag()
    assert G.isayclic()
    assert isident
    
