# -*- coding: utf-8 -*-
"""
A small library of functions for representing, analyzing and processing
acyclic directed mixed graphs for structural causal modelling.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

from queue import Empty
import grapl.util as util

class ADMGNode():
    """An ADMG node (random variable) object. The member variables, parents, children,
       bidirects are Sets of Strings, each string is the name of a random variable (node)
       in the graph. Bidirects are (nodes) connected to this one via a latent edge. The
       properties member is an optional Set of properties.
    """
    def __init__(self, parents={}, bidirects={}, children={}, properties={}):
        self.parents = set(parents)
        self.children = set(children)
        self.bidirects = set(bidirects)
        self.properties = set(properties)

class ADMG():
    """An ADMG graph object. The title member variable (String) is the name of graph
       (mostly for display purposes). The vars member is a Dict of String-ADMGNode pairs,
       representing his graph's random variables (ADMGNodes), indexed by node name (String).
    """
    def __init__(self, title=''):
        self.title = title
        self.vars = dict()

    def settitle(self, title=''):
        """Set graph title (String)."""
        self.title = title

    def addvar(self, name='', parents={}, bidirects={}, children={}, properties={}):
        """Add a new random variable to the graph with string name, with optional parents, bidirects
           children and/or properties. The name must not be empty.
           Each random variable (parents, bidirects, children) is a Set of Strings,
           each string representing a random variable in the graph.
        """
        if name != '':
            self.vars[name] = ADMGNode(parents,bidirects,children,properties)
        else:
            raise ValueError('admg:addvar: Invalid node name')

    def addedges(self, name='', parents={}, bidirects={}, children={}):
        """Add (optional) new edges connecting parents, bidirects and/or child random variables,
           to the random variable given in name (which must not be empty).
           Edges are Sets of Strings, each string representing a random variable in the graph.
        """
        if name != '':
            self.vars[name].parents = self.vars[name].parents.union(parents)
            self.vars[name].children = self.vars[name].children.union(children)
            self.vars[name].bidirects = self.vars[name].bidirects.union(bidirects)
        else:
            raise ValueError('admg:addedges: Invalid node name')

    def connect(self):
        """Ensure all parent/child/bidirect edges are consistent across the graph.
           Explicitly: every parent of a random variable, must have that variable
           as a child, and every random variable connected by a latent (bidrected)
           edge to another random variable, must have a bidirected edge to that
           first random variable.
        """
        for node in set(self.vars):
            for p in self.pa({node}):
                self.vars[p].children = self.ch({p}).union({node})
            for b in self.bi({node}):
                self.vars[b].bidirects = self.bi({b}).union({node})

    def display(self):
        """Pretty-print a text representation of the graph (ADMG) object."""
        nodes = set(self.vars)
        print('Title: ' + self.title)
        print('Vars: ' + util.csepstr(nodes))
        print('Parents:')
        for node in nodes:
            parents = self.pa({node})
            if len(parents) > 0:
                print(node,'<-', util.csepstr(parents))
        print('Children:')
        for node in nodes:
            children = self.ch({node})
            if len(children) > 0:
                print(node,'->', util.csepstr(children))
        print('Bidirects:')
        for node in nodes:
            bidirects = self.bi({node})
            if len(bidirects) > 0:
                print(node,'<-->', util.csepstr(bidirects))

    def nodes(self):
        """Return set of all random variables (nodes) in the graph.
           This is represented as a Set of Strings, each string being the name of
           a random variable.
        """
        return set(self.vars)

    def sub(self, subvars):
        """Find a subgraph of a graph. The subvars parameter is a Set of Strings,
           each string naming a random variable. If there are no such random variables
           in the graph, the result will be an empty graph. Returns an ADMG object.
        """
        sub = ADMG(self.title)
        for node in set(self.vars):
            if node in subvars:
                parents = self.pa({node}).intersection(subvars)
                bidirects = self.bi({node}).intersection(subvars)
                sub.vars[node] = ADMGNode(parents,bidirects,self.vars[node].properties)
        sub.connect()
        return sub
    
    def ch(self, nodes):
        """Find the union of all children of a set of random variables (nodes).
           The set of nodes is given as a Set of Strings, each string being the
           name of a random variable. Each such random variable must exist in the graph.
           Returns a Set of Strings, each string being the name of a random variable.
        """
        try:
            children = set()
            for node in nodes:
                children = children.union(self.vars[node].children)
            return children
        except KeyError as err:
            raise ValueError(f"admg.ch:Var '{node}' does not exist")

    def pa(self, nodes):
        """Find the union of all parents of a set of random variables (nodes).
           The set of nodes is given as a Set of Strings, each string being the
           name of a random variable. Each such random variable must exist in the graph.
           Returns a Set of Strings, each string being the name of a random variable.
        """
        try:
            parents = set()
            for node in nodes:
                parents = parents.union(self.vars[node].parents)
            return parents
        except KeyError as err:
            raise ValueError(f"admg.pa:Var '{node}' does not exist")

    def bi(self, nodes):
        """Find the union of all bidirected (latent) edges of a set of random variables (nodes).
           The set of nodes is given as a Set of Strings, each string being the
           name of a random variable. Each such random variable must exist in the graph.
           Returns a Set of Strings, each string being the name of a random variable.
        """
        try:
            bidirects = set()
            for node in nodes:
                bidirects = bidirects.union(self.vars[node].bidirects)
            return bidirects
        except KeyError as err:
            raise ValueError(f"admg.bi:Var '{node}' does not exist")
    
    def an(self, nodes):
        """Find the union of all ancestors of a given set of random variables (nodes),
           and the nodes themselves. The set of nodes is given as a Set of Strings,
           each string being the name of a random variable. Each such random variable must
           exist in the graph. Returns a Set of Strings, each string being the name
           of a random variable.
        """

        def _an_recurse(node, ancestors):
            parents = self.pa({node})
            ancestors = ancestors.union(parents)
            for parent in parents:
                ancestors = _an_recurse(parent, ancestors)
            return ancestors

        ancestors = set()
        for node in nodes:
            ancestors = ancestors.union(_an_recurse(node, set()))
        return ancestors.union(nodes)

    def de(self, nodes):
        """Find the union of all descendants of a given set of random variables (nodes),
           and the nodes themselves. The set of nodes is given as a Set of Strings,
           each string being the name of a random variable. Each such random variable must
           exist in the graph. Returns a Set of Strings, each string being the name
           of a random variable.
        """

        def _de_recurse(node, descendants):
            children = self.ch({node})
            descendants = descendants.union(children)
            for child in children:
                descendants = _de_recurse(child, descendants)
            return descendants

        descendants = set()
        for node in nodes:
            descendants = descendants.union(_de_recurse(node, set()))
        return descendants.union(nodes)
    
    def nd(self, nodes):
        """Find all non-descendants of a given set of random variables (nodes). The set of nodes
           is given as a Set of Strings, each string being the name of a random variable.
           Each such random variable must exist in the graph. Returns a Set of Strings, each string
           being the name of a random variable.
        """
        return self.nodes().difference(self.de(nodes))

    def dispa(self, node):
        """Find the district parents of a random variable (node) in a graph. The node is given
           as a String, which is the name of the random variable. This random variable must
           exist in the graph. Returns a Set of Strings, each string being the name of a random
           variable.
        """
        district = self.district(node)
        parents = self.pa(district)
        return parents.union(district.difference({node}))

    def mb(self, node):
        """Find the Markov blanket of a random variable (node) in a DAG. The node is given
           as a String, which is the name of the random variable, which must exist in the graph.
           If the graph is a DAG, returns a Set of Strings, each string being the name of a random
           variable, otherwise returns None.
        """
        if not self.isdag():
            return None
        else:
            children = self.ch({node})
            mb_nodes = self.pa({node}).union(children)
            for child in children:
                mb_nodes = mb_nodes.union(self.pa({child}))
            return mb_nodes.difference({node})

    def districts(self):
        """Compute all districts (c-components) of an ADMG. Returns a List of Set of Strings,
           each Set being a district of the graph. The sets are names of the random variables
           in each district.
        """
        components = []
        nodes = set(self.vars)
        while len(nodes) > 0:
            next_node = nodes.pop()
            component = self.district(next_node)
            components.append(component)
            nodes = nodes.difference(component)
        return components
            
    def district(self, node):
        """Find the district (c-components) of a random variable (node). The node is given
           as a String, which is the name of the random variable, which must exist in the graph.
           Returns a Set of Strings, each string being the name of a random variable in the
           district.
        """
        
        def _district_recurse(node, component):
            component = component.union({node})
            bidirects = self.bi({node})
            for bidirect in bidirects:
                if bidirect not in component:
                    component = _district_recurse(bidirect, component)
            return component
        
        return _district_recurse(node, set())

    def fix(self, node):
        """Apply the interventional fixing operation to a random variable (node) in the
           graph. This changes the graph (ADMG) object. The node is given as the string
           name of the random variable, which must exist in the graph. 
        """
        for parent in self.pa({node}):
            self.vars[parent].children = self.vars[parent].children.difference({node})
        self.vars[node].parents = set()
        for bidirect in self.bi({node}):
            self.vars[bidirect].bidirects = self.vars[bidirect].bidirects.difference({node})
        self.vars[node].bidirects = set()

    def isfixable(self, node):
        """Check whether it is possible to apply the interventional fixing operation
           to a random variable (node). The node is given as the string
           name of the random variable, which must exist in the graph. Returns True
           if fixable, False if not.
        """
        district = self.district(node)
        descend = self.de({node}).union({node})
        fix = district.intersection(descend)
        return (fix == {node})

    def isdag(self):
        """Returns True if graph is a DAG (no bidirects => no latent variables)."""
        bidirects = False
        for node in set(self.vars):
            bidirects = bidirects or (len(self.bi({node})) > 0)
        return (not bidirects)
    
    def fixable(self):
        """Returns set of all fixable nodes in the graph. Returns a Set of Strings,
           each string being the name of a fixable random variable in the graph.
        """
        nodes = set(self.vars)
        fix_set = set()
        for node in nodes:
            if self.isfixable(node):
                fix_set = fix_set.union({node})
        return fix_set

    def topsort(self):
        """List of random variables (nodes) sorted topologically. Returns a List of Strings,
           each string being the name of a random variable.
        """

        # Recursive depth-first search
        def _visit(node):
            if node not in visited:
                children = self.ch({node})
                for child in children:
                    _visit(child)
                visited.add(node)
                sorted.insert(0, node)

        visited = set()
        sorted = list()
        nodes = set(self.vars)

        # Initiate depth-first search from all possible starting points
        while True:
            unvisited = nodes.difference(visited)
            if unvisited:
                _visit(unvisited.pop())
            else:
                break

        return sorted

    def isayclic(self):
        """Returns True if the graph is acyclic. Attempt to sort nodes topologically,
           and uses this search to detect if graph is cyclic."""

        # Recursive depth-first search, marking nodes which are part of a cyclic path
        def _visit(node):
            if node in visited:
                return False
            if node in marked:
                return True
            marked.add(node)
            cycle = False
            children = self.ch({node})
            for child in children:
                cycle = _visit(child)
                if cycle:
                    break
            if cycle:
                return True
            marked.discard(node)
            visited.add(node)
            return False

        visited = set()
        marked = set()
        nodes = set(self.vars)

        # Initiate depth-first search from all possible starting points, with early termination
        # if cycle detected
        cyclic = False
        while not cyclic:
            unvisited = nodes.difference(visited)
            if unvisited:
                cyclic = _visit(unvisited.pop())
            else:
                break

        return not cyclic
