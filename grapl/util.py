# -*- coding: utf-8 -*-
"""
Various utility functions for the GRAPL library.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.admg as admg

def csepstr(nodes=[]):
    """Create a comma-separated string for a given list of nodes, given as a List of Strings,
       each string being the name of a random variable. If empty, then an empty string is
       returned. Returns a String.
    """
    return ','.join(str(s) for s in nodes)

def mrgstr(nodes):
    """Create a Latex marginalization string from a list of nodes, provided as a List of Strings,
       each string being the name of a random variable. If empty, then an empty string is
       returned. Returns a String."""
    if len(nodes) > 0:
        return '\sum_{' + csepstr(nodes) + '}'
    else:
        return ''

def phistr(nodes):
    """Create a Latex fixing operation string from a list of nodes, provided as a List of Strings,
       each string being the name of a random variable. If empty, then an empty string is
       returned. Returns a String.
    """
    if len(nodes) > 0:
        return '\phi_{' + csepstr(nodes) + '}'
    else:
        return ''

def probstr(nodes):
    """Create a Latex distribution string from a list of nodes, given as a List of Strings,
       each string being the name of a random variable. If empty, then an empty string is
       returned. Returns a String.
    """
    if len(nodes) > 0:
        return 'p(' + csepstr(nodes) + ')'
    else:
        return ''

def probcndstr(nodes=[], cnd_nodes=[]):
    """Create a Latex conditional distribution string from a list of nodes and conditional nodes.
       The parameters are Lists of Strings; nodes being the list of conditioned random variables,
       and cnd_nodes being the list of their parents; if the latter is empty then the distribution
       string is not conditional, and if the former is empty then an empty string is returned.
       Returns a String.
    """
    prob_str = ''
    if len(nodes) > 0:
        prob_str = prob_str + 'p(' + csepstr(nodes)
        if len(cnd_nodes) > 0:
            prob_str = prob_str + '|' + csepstr(cnd_nodes)
        prob_str = prob_str + ')'
    return prob_str
