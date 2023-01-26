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
import copy as cp


def idfixing(G, X={}, Y={}):
    """Implementation of Richardson et al.'s Theorem 60 for causal effect identification
       from given ADMG. See: T.S. Richardson, J.M. Robins, I. Shpitser, 2012:
       "Nested Markov properties for acyclic directed mixed graphs", UAI'12.

       Parameters:
       G (ADMG)           -- A mixed causal graph object
       X (Set of Strings) -- Set of interventional variables, where each string is a random variable name (must not be empty)
       Y (Set of Strings) -- Set of variables on which interventional expression is represented, where each string is a random
                             variable name (if empty, then all variables in G except for those in X)

       Returns: (String, Eqn, Boolean)
       If interventional distribution is identifiable, returns the identified equation as a Latex
       String, with the corresponding Eqn object, and True. Otherwise, returns '', None, False
    """

    # Must supply the set of interventional (cause) variables
    if not X:
        return '', None, False

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

    # Find subgraph G(V)_Y* and its districts
    G_sub_Ys = G.sub(Y_star)
    D_y = G_sub_Ys.districts()

    # For each district, determine post-fixed ADMG and corresponding expression
    D_exprs = []
    for D in D_y:

        # Nodes to fix for this district, fixing transform: phi_{V\D}(p(V;G_V))
        VdD = V.difference(D)

        # Initialize accumulated ADMG, G(V), and district expression, p(V)
        G_fixed = cp.deepcopy(G)
        D_expr = expr.Expr(num=[V])

        # Try to apply fixing sequence
        nodes_fix = VdD.copy()
        while (len(nodes_fix) > 0):
            
            # Find a valid fixing node in sequence
            fixable_nodes = {node for node in nodes_fix if G_fixed.isfixable(node)}
            node_valid = (len(fixable_nodes) > 0)
            if node_valid:
                fix = fixable_nodes.pop()
            else:
                break
            
            # Fix or marginalize (Corollary 33) the accumulated expression
            dispa_fix = G_fixed.dispa(fix)
            fix_ch = G_fixed.ch({fix})
            if len(fix_ch) == 0:
                D_expr.fixmarginal(fix)
            else:
                D_expr.fix(fix, dispa_fix)

            # Fix node (remove all its incoming edges) in accumulated ADMG,
            # remove from remaining fixing sequence
            G_fixed.fix(fix)
            nodes_fix = nodes_fix.difference({fix})

        if not node_valid:
            break
        
        D_expr.simplify()
        D_exprs.append(D_expr)

    # At this point, if no valid fixing node, intervention is not identifiable
    if not node_valid:
        return '', None, False
    else:
        Y_mrg = Y_star.difference(Y)
        id_expr = expr.Expr(mrg=Y_mrg)
        for ex in D_exprs:
            clash_vars = id_expr.mrg.difference(id_expr.mrg.difference(ex.mrg))
            for old_var in clash_vars:
                new_var = old_var + chr(39)
                ex.subsvar(old_var, new_var)
            id_expr.num = id_expr.num + ex.num
            id_expr.den = id_expr.den + ex.den
            id_expr.mrg = id_expr.mrg.union(ex.mrg)

        id_expr.simplify()

        # Construct Latex equation for the identified interventional distribution
        lhs = expr.Expr()
        lhs.addvars(num=Y, dov=A)
        id_eqn = eqn.Eqn(lhs, id_expr)
        id_str = id_eqn.tocondstr()

        return id_str, id_eqn, True


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
