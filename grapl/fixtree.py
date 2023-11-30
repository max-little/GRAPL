# -*- coding: utf-8 -*-
"""
Implement a Breadth-First Search (BFS) tree structure to search possible fixing
sequences to fix a given causal graph (ADMG) if identifiable.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""
import grapl.expr as expr
import grapl.eqn as eqn
import copy as cp
import grapl.util as util
from random import randint

class FixNode:
    """
    A FixNode object to store the intermediate variables in the recursive tree growing and search process.
    The member variables are: G_fixed (fixed graph), an ADMG object, the updated graph after a certain node is fixed;
    D_expr (district expression), an Expr object, tracking the accumlated district expression; D_exprs (district
    expressions), a List of Expr object, each Expr object is the district expression for a district; nodes_fix (nodes 
    to be fixed), a List of Strings, each string being the name of the node that needs to be fixed; fixing_seq (fixing
    sequence information), A String, recording the fixing sequence and process details; childern (child nodes), a List
    of Node object, being the child nodes of a given Node.
    """
    def __init__(self, G_fixed, nodes_fix, fix_node_seq, fix_graph_seq):

        self.G_fixed = G_fixed
        self.nodes_fix = nodes_fix
        self.fix_node_seq = fix_node_seq
        self.fix_graph_seq = fix_graph_seq
        self.if_dfixable = False
        self.childern = []


class FixSeqTree:
    """
    A Class object that generates a tree to search all possible fixing sequences for a district.
    
    Parameters:
    G (ADMG)                    -- A mixed causal graph object
    X (Set of Strings)          -- Set of interventional variables, where each string is a random variable name (must not be empty)
    Y (Set of Strings)          -- Set of variables on which interventional expression is represented, where each string is a random
                                   variable name (if empty, then all variables in G except for those in X)
    """
            
    def __init__(self, G, X = {}, Y = {}, degrad = False):
        self.G = G
        self.X = X
        self.Y = Y
        self.degrad = degrad
        
        self.V = G.nodes()
        self.A = X
        
        self.id_eqn_all = []
        self.id_str_all = []
        self.num_all = []
        self.denum_all = []
        
        # Must supply the set of interventional (cause) variables
        if not self.X:
            print("Not idenfiable without interventional variable assigned!")
            return
        
    def expand_tree_in_D(self, cur_node):
        """
        A recursive function for tree growing to search all possible 
        fixing orders within a given district D.
        
        Parameters:
        cur_node (FixNode): The current search node.
        """

        G_fixed = cur_node.G_fixed
        nodes_fix = cur_node.nodes_fix
        fix_node_seq = cur_node.fix_node_seq
        fix_graph_seq = cur_node.fix_graph_seq
        
        if len(nodes_fix) > 0:

            # Find a valid fixing node in sequence
            fixable_nodes = [node for node in nodes_fix if G_fixed.isfixable(node)]
            self.node_valid = (len(fixable_nodes) > 0)
            if self.node_valid:

                # Fix the graph by picking every possible fixable node
                if not self.degrad:
                    fixable_node_pool = [i for i in range(len(fixable_nodes))]
                else:
                    fixable_node_pool = [randint(0, len(fixable_nodes)-1)]
                for fix_node_index in fixable_node_pool:             
                    G_fixed_cp = cp.deepcopy(G_fixed)
                    nodes_fix_cp = cp.deepcopy(nodes_fix)
                    fix_node_seq_cp = cp.deepcopy(fix_node_seq)
                    fix_graph_seq_cp = cp.deepcopy(fix_graph_seq)
                    fix = fixable_nodes[fix_node_index]

                    # Fix node (remove all its incoming edges) in accumulated ADMG
                    G_fixed_cp.fix(fix)

                    # Remove from remaining fixing sequence
                    nodes_fix_cp = nodes_fix_cp.difference({fix})
                    fix_node_seq_cp.append(fix)
                    fix_graph_seq_cp.append(G_fixed_cp)
                    child_node = FixNode(G_fixed = G_fixed_cp, nodes_fix = nodes_fix_cp, 
                                         fix_node_seq = fix_node_seq_cp, fix_graph_seq = fix_graph_seq_cp)             
                    cur_node.childern.append(child_node)
            else:
                return
        else:
            if self.node_valid:
                cur_node.if_dfixable = True
            else:
                cur_node.if_dfixable = False
            return
        
        # Expand the tree for every fixable node
        for child in cur_node.childern:
            self.expand_tree_in_D(child)
        
    def get_leaf(self, tree):
        """
        A recursive function to retrieve leaf nodes in the fixing tree structure.

        Parameter:
        tree (FixNode): The tree to search the leaf modes.
        
        Returns:
        A List of leaf Nodes in the tree.
        """

        if len(tree.childern) == 0:
            if tree.if_dfixable:
                return [tree]
            else:
                return []
        
        leaf = []
        for child in tree.childern:
            leaf.extend(self.get_leaf(child))
        
        return leaf
    
    def searchdfixseq(self, D):
        fix_node_seqs = []
        fix_graph_seqs = []
        
        # Nodes to fix for this district, fixing transform: phi_{V\D}(p(V;G_V))
        VdD = self.V.difference(D)
        
        # Initialize accumulated ADMG, G(V), and district expression, p(V)
        G_fixed = cp.deepcopy(self.G)
        nodes_fix = VdD.copy()
        root_node = FixNode(G_fixed = G_fixed, nodes_fix = nodes_fix, fix_node_seq = [], fix_graph_seq = [])
        valid_leaf = [root_node]
        for i, leaf in enumerate(valid_leaf):
            leaf.G_fixed = G_fixed
            leaf.nodes_fix = nodes_fix
            self.expand_tree_in_D(cur_node = leaf)
            # retrieve leaf nodes for this district
            if i == 0:
                valid_leaf = self.get_leaf(leaf)
            else:
                valid_leaf.extend(self.get_leaf(leaf))
        
        for leaf in valid_leaf:
            fix_node_seqs.append(leaf.fix_node_seq)
            fix_graph_seqs.append(leaf.fix_graph_seq)
            
        if len(fix_node_seqs) > 0:
            return fix_node_seqs, fix_graph_seqs, True
        else:
            return fix_node_seqs, fix_graph_seqs, False