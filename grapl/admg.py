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
    """An ADMG node (random variable) object."""
    def __init__(self, parents={}, bidirects={}, children={}, properties={}):
        self.parents = set(parents)
        self.children = set(children)
        self.bidirects = set(bidirects)
        self.properties = set(properties)

class ADMG():
    """An ADMG graph object."""
    def __init__(self, title=''):
        self.title = title
        self.vars = dict()

    def settitle(self, title=''):
        """Set graph title."""
        self.title = title

    def addvar(self, name='', parents={}, bidirects={}, children={}, properties={}):
        """Add a variable to a graph with optional parents, bidirects and/or children."""
        if name != '':
            self.vars[name] = ADMGNode(parents,bidirects,children,properties)
        else:
            raise ValueError('admg:addvar: Invalid node name')

    def addedges(self, name='', parents={}, bidirects={}, children={}):
        """Add parents, bidirects and/or child edges to a graph."""
        if name != '':
            self.vars[name].parents = self.vars[name].parents.union(parents)
            self.vars[name].children = self.vars[name].children.union(children)
            self.vars[name].bidirects = self.vars[name].bidirects.union(bidirects)
        else:
            raise ValueError('admg:addedges: Invalid node name')

    def connect(self):
        """Ensure all parent/child/bidirect edges are consistent across the graph."""
        for node in set(self.vars):
            for p in self.pa({node}):
                self.vars[p].children = self.ch({p}).union({node})
            for b in self.bi({node}):
                self.vars[b].bidirects = self.bi({b}).union({node})

    def display(self):
        """Pretty-print a text representation of the graph."""
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
        """Return set of all nodes (variables)."""
        return set(self.vars)

    def sub(self, subvars):
        """Find a subgraph of a graph."""
        sub = ADMG(self.title)
        for node in set(self.vars):
            if node in subvars:
                parents = self.pa({node}).intersection(subvars)
                bidirects = self.bi({node}).intersection(subvars)
                sub.vars[node] = ADMGNode(parents,bidirects,self.vars[node].properties)
        sub.connect()
        return sub
    
    def ch(self, nodes):
        """Find the union of all children of a set of nodes."""
        try:
            children = set()
            for node in nodes:
                children = children.union(self.vars[node].children)
            return children
        except KeyError as err:
            raise ValueError(f"admg.ch:Var '{node}' does not exist")

    def pa(self, nodes):
        """Find the union of all parents of a set of nodes."""
        try:
            parents = set()
            for node in nodes:
                parents = parents.union(self.vars[node].parents)
            return parents
        except KeyError as err:
            raise ValueError(f"admg.pa:Var '{node}' does not exist")

    def bi(self, nodes):
        """Find the union of all bidirected edges of a set of nodes."""
        try:
            bidirects = set()
            for node in nodes:
                bidirects = bidirects.union(self.vars[node].bidirects)
            return bidirects
        except KeyError as err:
            raise ValueError(f"admg.bi:Var '{node}' does not exist")
    
    def an(self, nodes):
        """Find the union of all ancestors of a given set of nodes, and the nodes themselves."""

        def an_recurse(node, ancestors):
            parents = self.pa({node})
            ancestors = ancestors.union(parents)
            for parent in parents:
                ancestors = an_recurse(parent, ancestors)
            return ancestors

        ancestors = set()
        for node in nodes:
            ancestors = ancestors.union(an_recurse(node, set()))
        return ancestors.union(nodes)

    def de(self, nodes):
        """Find the union of all descendants of a given set of nodes, and the nodes themselves."""

        def de_recurse(node, descendants):
            children = self.ch({node})
            descendants = descendants.union(children)
            for child in children:
                descendants = de_recurse(child, descendants)
            return descendants

        descendants = set()
        for node in nodes:
            descendants = descendants.union(de_recurse(node, set()))
        return descendants.union(nodes)
    
    def nd(self, nodes):
        """Find all non-descendants of a given set of nodes."""
        return self.nodes().difference(self.de(nodes))

    def dispa(self, node):
        """Find the district parents of a node in an ADMG."""
        district = self.district(node)
        parents = self.pa(district)
        return parents.union(district.difference({node}))

    def mb(self, node):
        """Find the Markov blanket of a node in a DAG, returns None if graph is not a DAG."""
        if not self.isdag():
            return None
        else:
            children = self.ch({node})
            mb_nodes = self.pa({node}).union(children)
            for child in children:
                mb_nodes = mb_nodes.union(self.pa({child}))
            return mb_nodes.difference({node})

    def districts(self):
        """Compute the set of all districts (c-components) of an ADMG."""
        components = []
        nodes = set(self.vars)
        while len(nodes) > 0:
            next_node = nodes.pop()
            component = self.district(next_node)
            components.append(component)
            nodes = nodes.difference(component)
        return components
            
    def district(self, node):
        """Find the district (c-components) of a node."""
        
        def district_recurse(node, component):
            component = component.union({node})
            bidirects = self.bi({node})
            for bidirect in bidirects:
                if bidirect not in component:
                    component = district_recurse(bidirect, component)
            return component
        
        return district_recurse(node, set())

    def fix(self, node):
        """Apply the interventional fixing operation to a node."""
        for parent in self.pa({node}):
            self.vars[parent].children = self.vars[parent].children.difference({node})
        self.vars[node].parents = set()
        for bidirect in self.bi({node}):
            self.vars[bidirect].bidirects = self.vars[bidirect].bidirects.difference({node})
        self.vars[node].bidirects = set()

    def isfixable(self, node):
        """Check whether it is possible to apply the interventional fixing operation
           to a node."""
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
        """Returns set of all fixable nodes in the graph."""
        nodes = set(self.vars)
        fix_set = set()
        for node in nodes:
            if self.isfixable(node):
                fix_set = fix_set.union({node})
        return fix_set

    def topsort(self):
        """Return list of nodes sorted topologically."""

        # Recursive depth-first search
        def visit(node):
            if node not in visited:
                children = self.ch({node})
                for child in children:
                    visit(child)
                visited.add(node)
                sorted.insert(0, node)

        visited = set()
        sorted = list()
        nodes = set(self.vars)

        # Initiate depth-first search from all possible starting points
        while True:
            unvisited = nodes.difference(visited)
            if unvisited:
                visit(unvisited.pop())
            else:
                break

        return sorted

    def isayclic(self):
        """Returns True if the graph is acyclic. Attempt to sort nodes topologically,
           and uses this search to detect if graph is cyclic."""

        # Recursive depth-first search, marking nodes which are part of a cyclic path
        def visit(node):
            if node in visited:
                return False
            if node in marked:
                return True
            marked.add(node)
            cycle = False
            children = self.ch({node})
            for child in children:
                cycle = visit(child)
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
                cyclic = visit(unvisited.pop())
            else:
                break

        return not cyclic
