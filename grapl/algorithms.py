# -*- coding: utf-8 -*-
"""
Implements algorithms for processing of structural causal models represented
by acyclic directed mixed graphs (ADMGs).

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.expr as expr
import grapl.eqn as eqn
import grapl.condind as condind
from grapl.fixtree import FixSeqTree
import copy as cp
import numpy as np
from itertools import product
import warnings


def preintvar(G, X = {}, Y = {}):
    """
    Find Y*, the set of potential causes of effects Y, after taking into account interventions X in an ADMG G.
    
    Parameters:
    G (ADMG)           -- A mixed causal graph object
    X (Set of Strings) -- Set of interventional variables, where each string is a random variable name (must not be empty)
    Y (Set of Strings) -- Set of variables on which interventional expression is represented, where each string is a random
                          variable name.
    Returns:
    Y_star (Set of Strings) -- Set of potential causes of those effects after removing effect vars and the desired intervention vars
    """

    # Collate variable names
    V = G.nodes()

    # Handle effects
    if not Y:
        Y = V.difference(X)

    # Find Y* = an_G(V)_{V\A}(Y)    
    A = X
    VdA = V.difference(A)
    G_sub_VdA = G.sub(VdA)
    Y_star = G_sub_VdA.an(Y)
    
    return Y_star


def fixseqsdistr(G, D, X = {}, Y = {}, degrad = False):
    """
    Given a graph G, district D and desired effect and intervention variables, find all possible fixing sequences and
    corresponding fixed graphs.
    
    Parameters:
    G (ADMG)           -- A mixed causal graph object
    D (Set of Strings) -- A district, where each string is a random variable name
    X (Set of Strings) -- Set of interventional variables, where each string is a random variable name (must not be empty).
    Y (Set of Strings) -- Set of variables on which interventional expression is represented, where each string is a random
                          variable name.
    degrad (Boolean)   -- If True, degrade the search tree to linear structure of searching. Default False.
    
    Returns:
    fix_node_seqs (List of List of Strings) -- A list of all possible fixing sequences, where each sequence is a list of strings
                                               representing the name of the node to be fixed. 
    fix_graph_seqs (List of List of ADMG)   -- A list of all corresponding fixed graphs, where each graph is an ADMG object.
    """

    # Initialise fixing sequence searcher
    seq_searcher = FixSeqTree(G = G, X = X, Y = Y, degrad = degrad)

    # Find all possible fixing sequences and their fixed graphs
    fix_node_seqs, fix_graph_seqs, identifiable = seq_searcher.searchdfixseq(D)
    
    if identifiable:
        return fix_node_seqs, fix_graph_seqs, identifiable
    else:
        return None, None, identifiable


def exprfixseq(G, fix_node_seq, fix_graph_seq):
    """
    Given the graph G, a fixing node sequence and its corresponding the fixed graphs, compute a corresponding
    distributional expression.
    
    Parameters:
    G (ADMG)                       -- A mixed causal graph object
    fix_node_seq (List of Strings) -- A list of strings representing the name of the node to be fixed.
    fix_graph_seq (List of ADMG)   -- A list of corresponding fixed graphs, where each graph is an ADMG object.
    
    Returns:
    fixexpr_str (String) -- The simplied expression as a string
    fixexpr (Expr)       -- The expression as an Expr object
    """

    # Start with the joint distribution for all nodes
    fixexpr = expr.Expr(num=[G.nodes()])
    graph_seq = [G] + fix_graph_seq[:-1]
    for i in range(len(fix_node_seq)):
        fix = fix_node_seq[i]
        G_fixed = graph_seq[i]

        # Fix or marginalize (Corollary 33) the accumulated expression
        # explain it more based on the paper
        dispa_fix = G_fixed.dispa(fix)
        fix_ch = G_fixed.ch({fix})

        # if no children, then marginalize...
        if len(fix_ch) == 0:
            fixexpr.fixmarginal(fix)

        # if there are children, then fix... 
        else:
            fixexpr.fix(fix, dispa_fix)

    fixexpr.simplify()
    fixexpr_str = fixexpr.tostr()
    return fixexpr_str, fixexpr 


def idfixall(G, X = {}, Y ={}, mode = "shortest", greedy = True):
    """
    Search for interventional probability expressions for the given causal graph
    and retrieve results based on the selected mode.
    A general implementation of Richardson et al.'s Theorem 60 for causal effect 
    identification from given ADMG. See: T.S. Richardson, J.M. Robins, I. Shpitser, 2012:
    "Nested Markov properties for acyclic directed mixed graphs", UAI'12.
    
    Parameters:
    G (ADMG)           -- A mixed causal graph object
    X (Set of Strings) -- Set of interventional variables, where each string is a random variable name (must not be empty)
    Y (Set of Strings) -- Set of variables on which interventional expression is represented, where each string is a random
                          variable name (if empty, then all variables in G except for those in X)
    mode (String)      -- The mode for selecting interventional equations {'shortest', 'mostmrg', 'random', or 'all'}, default 'shortest'.
    greedy (Boolean)   -- If True, use greedy strategy to select the shortest or most marginalized expression. Otherwise, global search applied. Default True.
    
    Returns:
    Depending on the selected 'mode', it returns different results:
    - 'shortest': a Latex string of identified equation with the least number of distributions, its corresponding Eqn object, tracking information of fixing, and a boolean indicating if the causal effect is identifiable.
    - 'mostmrg': a Latex string of identified equation with the most number marginalized variables, its corresponding Eqn object, tracking information of fixing, and a boolean indicating if the causal effect is identifiable.
    - 'random': A randomly selected Latex string of identified equation, its corresponding Eqn object, tracking information of fixing, and a boolean indicating if the causal effect is identifiable.
    - 'all': A list of all Latex string of identified equations, a list of their corresponding Eqn objects, tracking information of fixing, and a boolean indicating if the causal effect is identifiable.
    If the causal effect is not identifiable: returns '', None, None, False
    """

    if mode not in ["shortest", "mostmrg", "random", "all"]:
        raise ValueError("'{}' is not a valid value for mode. Valid values are: 'shortest', 'mostmrg', 'random', 'all'. Default: 'shortest'. ".format(mode))
    
    # Find the set of potential causes of those effects after removing effect var and the desired intervention var
    Y_star = preintvar(G, X, Y)

    # Find the districts of G_Y*
    G_sub_Ystar = G.sub(Y_star)
    D_y = G_sub_Ystar.districts()
    D_exprs = []
    
    for D in D_y:
        # For each district, find fixing sequences and corresponding fixed graphs
        fix_node_seqs, fix_graph_seqs, identifiable = fixseqsdistr(G = G, D = D, X = X, Y = Y, degrad = mode == "random")
        if not identifiable:
            return None, None, False
        D_exprs.append([])

        # According to the selected mode and greedy strategy, formulate the expression sequences
        if mode == "shortest" and greedy:
            shortest_dist_cnt = np.inf
            shortest_D_expr = None
            for i in range(len(fix_node_seqs)):
                _, D_expr  = exprfixseq(G, fix_node_seqs[i], fix_graph_seqs[i])
                D_num = D_expr.num
                D_den = D_expr.den
                dist_cnt = len(D_num) + len(D_den)
                if dist_cnt < shortest_dist_cnt:
                    shortest_dist_cnt = dist_cnt
                    shortest_D_expr = D_expr
            D_exprs[-1].append(shortest_D_expr)
        elif mode == "mostmrg" and greedy:
            most_mrg_cnt = -1
            most_mrg_D_expr = None
            for i in range(len(fix_node_seqs)):
                _, D_expr  = exprfixseq(G, fix_node_seqs[i], fix_graph_seqs[i])
                mrg_cnt = len(D_expr.mrg)
                if mrg_cnt > most_mrg_cnt:
                    most_mrg_cnt = mrg_cnt
                    most_mrg_D_expr = D_expr
            D_exprs[-1].append(most_mrg_D_expr)
        else:
            D_expr_num = []
            D_expr_den = []
            for i in range(len(fix_node_seqs)):
                _, D_expr  = exprfixseq(G, fix_node_seqs[i], fix_graph_seqs[i])
                D_num = D_expr.num
                D_den = D_expr.den
                if (D_num not in D_expr_num) and (D_den not in D_expr_den):
                    D_exprs[-1].append(D_expr)
                    D_expr_num.append(D_num)
                    D_expr_den.append(D_den)
    D_exprs_cartprod = list(product(*D_exprs))

    id_str_all =[]
    id_eqn_all = []

    # Combining all expressions of fixing districts
    for exprs in D_exprs_cartprod:
        id_expr = expr.Expr()
        id_expr.combine(exprs)
        Y_mrg = Y_star.difference(Y)
        id_expr.mrg = id_expr.mrg.union(Y_mrg)
        id_expr.simplify()
        lhs = expr.Expr()
        lhs.addvars(num=Y, dov=X)
        id_eqn = eqn.Eqn(lhs, id_expr)
        id_str = id_eqn.tocondstr()
        id_str_all.append(id_str)
        id_eqn_all.append(id_eqn)
    
    # Modes for selecting identified equations. Options: Shortest, Most Marginalized, Random, All
    if mode == "shortest":
        shortest_dist_cnt = np.inf
        shortest_id_eqn = None
        shortest_id_str = None
        for i,id_eqn in enumerate(id_eqn_all):
            w_denom = id_eqn.rhs.den
            w_nom = id_eqn.rhs.num
            dist_cnt = len(w_denom) + len(w_nom)
            if dist_cnt < shortest_dist_cnt:
                shortest_dist_cnt = dist_cnt
                shortest_id_eqn = id_eqn
                shortest_id_str = id_str_all[i]
        return shortest_id_str, shortest_id_eqn, identifiable
    
    elif mode == "mostmrg":
        most_mrg_cnt = -1
        most_mrg_id_eqn = None
        most_mrg_id_str = None
        for i, id_eqn in enumerate(id_eqn_all):
            mrg_cnt = len(id_eqn.rhs.mrg)
            if mrg_cnt > most_mrg_cnt:
                most_mrg_cnt = mrg_cnt
                most_mrg_id_eqn = id_eqn
                most_mrg_id_str = id_str_all[i]
        return most_mrg_id_str, most_mrg_id_eqn, identifiable
    
    elif mode == "random":
        rd_id_eqn = id_eqn_all[0]
        rd_id_str = id_str_all[0]
        return rd_id_str, rd_id_eqn, identifiable
    
    elif mode == "all":
        return id_str_all, id_eqn_all, identifiable


def idfixing(G, X={}, Y={}):
    """A special case of the idfixall function which returns a randomly selected Latex string of 
       identified equation, its corresponding Eqn object, and a boolean indicating if the causal effect is identifiable.
        A general implementation of Richardson et al.'s Theorem 60 for causal effect identification from given ADMG.
        See: T.S. Richardson, J.M. Robins, I. Shpitser, 2012: "Nested Markov properties for acyclic directed mixed graphs", UAI'12.

       Parameters:
       G (ADMG)           -- A mixed causal graph object
       X (Set of Strings) -- Set of interventional variables, where each string is a random variable name (must not be empty)
       Y (Set of Strings) -- Set of variables on which interventional expression is represented, where each string is a random
                             variable name (if empty, then all variables in G except for those in X)

       Returns: (String, Eqn, Boolean)
       If interventional distribution is identifiable, returns a randomly selected identified 
       equation as a Latex String, with the corresponding Eqn object, and True. Otherwise, 
       returns '', None, False
    """

    id_str, id_eqn, identifiable = idfixall(G, X, Y, mode = "random")
    if identifiable:
        return id_str, id_eqn, identifiable
    else:
        return "", None, identifiable


def dagfactor(G, Y={}, simplify=True):
    """Factorized distribution for DAGs (ADMGs with no bidirects).

       Parameters:
       G (ADMG)           -- DAG representing the causal graph (must not have bidirects)
       Y (Set of Strings) -- Set of effect variables (if empty, all variables in G)
       simplify (Boolean) -- If True, final expression is explicitly simplified

       Returns: (String, Expr, Boolean)
       If G is a DAG, factored, Latex-format equation String, corresponding Eqn object, and True.
       If G is not a DAG then '', None, False.
    """

    # Check to make sure there are no latent variables
    if not G.isdag():
        return "", None, False

    # Collate variable names
    V = G.nodes()

    # Handle marginals
    if not Y:
        Y = V
    Y_mrg = V.difference(Y)
    fac_expr = expr.Expr(mrg=Y_mrg)

    # Construct factorized expression over all nodes
    for node in V:
        fac_expr.num = fac_expr.num + [G.pa({node}).union({node})]
        fac_expr.den = fac_expr.den + [G.pa({node})]

    # Final expression simplification
    if simplify:
        fac_expr.simplify()

    # Construct Latex expression
    lhs = expr.Expr()
    lhs.addvars(num=Y)
    fac_eqn = eqn.Eqn(lhs, fac_expr)
    dist_str = fac_eqn.tocondstr()

    return dist_str, fac_eqn, True


def truncfactor(G, X={}, Y={}, prefactor=True):
    """Truncated factorization ("g-formula") for DAGs (ADMGs with no bidirects).

       Parameters:
       G (ADMG) -- DAG object representing the causal graph (must not have bidirects)
       X (Set of Strings)  -- Interventional variables where each string is a random variable name (must not be empty)
       Y (Set of Strings)  -- Effect variables, each sring is a random variable name (if empty, all variables in G other than the set X)
       prefactor (Boolean) -- If True, joint distribution is chain factored before fixing

       Returns: (String, Expr, Boolean)
       If G is a DAG, factored, Latex-formated interventonal distribution string, corresponding Eqn object,
       and True. Otherwise, returns '', None, False.
    """

    # Must supply the set of interventional variables
    if not X:
        return '', None, False

    # Check there are no latent variables
    if not G.isdag():
        return '', None, False

    # Collate variable names
    V = G.nodes()

    # Initialize accumulated DAG and interventional expression, p(V)
    G_fixed = cp.deepcopy(G)
    if prefactor:
        dummy, D_eqn, dummy_bool = dagfactor(G)        # Start with factorized joint distribution
        D_expr = cp.copy(D_eqn.rhs)
        # dummy, D_expr, dummy_bool = dagfactor(G)        # Start with factorized joint distribution
    else:
        D_expr = expr.Expr(num=[V])                     # Start with the unfactored joint distribution

    # Apply fixing sequence
    nodes_fix = cp.deepcopy(X)
    while (len(nodes_fix) > 0):
        
        # Select next fixing node in sequence
        fix = nodes_fix.pop()
        
        # Fix or marginalize (Corollary 33) the accumulated expression
        dispa_fix = G_fixed.dispa(fix)
        fix_ch = G_fixed.ch({fix})
        if len(fix_ch) == 0:
            D_expr.fixmarginal(fix)
        else:
            D_expr.fix(fix, dispa_fix)

        # Fix node (remove all its incoming edges) in accumulated DAG,
        # remove from remaining fixing sequence
        G_fixed.fix(fix)
        nodes_fix = nodes_fix.difference({fix})

    # Simplify final expression
    D_expr.simplify()

    # Handle marginals
    if not Y:
        Y = V.difference(X)
    Y_mrg = (V.difference(Y)).difference(X)
    id_expr = expr.Expr(mrg=Y_mrg)

    # Avoid naming clashes with dummy marginal variables
    clash_vars = id_expr.mrg.difference(id_expr.mrg.difference(D_expr.mrg))
    for old_var in clash_vars:
        new_var = old_var + chr(39)
        D_expr.subsvar(old_var, new_var)
    id_expr.num = id_expr.num + D_expr.num
    id_expr.den = id_expr.den + D_expr.den
    id_expr.mrg = id_expr.mrg.union(D_expr.mrg)

    # Final simplification
    id_expr.simplify()

    # Construct Latex expression
    lhs = expr.Expr()
    lhs.addvars(num=Y, dov=X)
    id_eqn = eqn.Eqn(lhs, id_expr)
    id_str = id_eqn.tocondstr()

    return id_str, id_eqn, True


def admgfactor(G, Y={}):
    """Factorized distribution for ADMGs.

       Parameters:
       G (ADMG)           -- A mixed causal graph
       Y (Set of Strings) -- Set of variables on which expression is represented, where each string is a random variable name
                             (if empty, then all variables in G)

       Returns: (String, Eqn)
       Factored Latex-format distribution equation string, corresponding Eqn object.
    """

    # Collate variable names
    V = G.nodes()

    # Find a topological ordering of the observed nodes
    sort_nodes = G.topsort()

    # Handle marginals
    if not Y:
        Y = V
    Y_mrg = V.difference(Y)
    fac_expr = expr.Expr(mrg=Y_mrg)

    # Factor distribution for each observed node in turn
    for node in V:

        # Set of all nodes preceding and including this one in topological order
        nodes_prec = set(sort_nodes[:sort_nodes.index(node) + 1])

        # Subgraph of all such ordered nodes
        G_prec = G.sub(nodes_prec)

        # District for node in G_prec
        D_node = G_prec.district(node)

        # Restrict ordered nodes to those within the district
        D_node_prec = nodes_prec.intersection(D_node)

        # Find (inclusive) parents of all such nodes
        D_parents = G.pa(D_node_prec).union(D_node_prec)

        # Remove node to find 'shielding' nodes
        shield = D_parents.difference({node})

        fac_expr.num = fac_expr.num + [shield.union({node})]
        fac_expr.den = fac_expr.den + [shield]

    # Final expression simplification
    fac_expr.simplify()

    # Construct Latex expression
    lhs = expr.Expr()
    lhs.addvars(num=Y)
    fac_eqn = eqn.Eqn(lhs, fac_expr)
    dist_str = fac_eqn.tocondstr()

    return dist_str, fac_eqn


def localmarkov(G):
    """Compute all local Markov independences for DAGs (ADMGs with no bidirects).

       Parameters:
       G (ADMG) -- DAG object representing the causal graph (must not have bidirects)

       Returns: (Set of Strings, Boolean)
       If G is a DAG, set of Unicode-format strings representing conditional independences,
       a CondIndSet object containing these, and True. Otherwise returns empty string,
       empty set and False.
    """

    # Check to make sure there are no latent variables
    if not G.isdag():
        return "", set(), False

    # Collate variable names
    nodes = set(G.vars)

    cis = condind.CondIndSet()
    for node in nodes:
        parents = G.pa({node})
        non_descend = G.nd({node})
        conds = non_descend.difference({node}).difference(parents)
        if conds:
            ci = condind.CondInd({node},conds,parents)
            cis.ciset.add(ci)

    return cis, cis.tostr(),  True


def dseparate(G,X,Y,Z={}):
    """For DAGs, test whether a set of source nodes are conditionally independent from a set of targets,
       conditioned on another set of nodes.

       Parameters:
       G (ADMG)            -- DAG object representing the causal graph (must not have bidirects)
       X (Set of Strings)  -- Source variables where each string is a random variable name
       Y (Set of Strings)  -- Target variables, each string is a random variable name
       Z (Set of Strings)  -- Conditioned variables, each string a random variable name

       Returns: (Boolean, String, Boolean)
       If G is a DAG, last parameter is True, otherwise False. If X is d-separated from Y given Z,
       first parameter is True, second parameter is a CondInd object representing the conditional
       independence, third parameter contains the corresponding Unicode-format string.
    """

    # Check to make sure there are no latent variables
    if not G.isdag():
        return None, None, "", False

    dsep = True
    for nodeX in X:
        for nodeY in Y:
            dsep = dsep and G.isdsep(nodeX,nodeY,Z)

    if dsep:
        conds = Y.difference(X)
        ci = condind.CondInd(X,conds,Z)
        return True, ci, ci.tostr(), True
    else:
        return False, None, "", True
