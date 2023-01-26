# -*- coding: utf-8 -*-
"""
Algorithms for representing sets of conditional independences in ADMGs.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.util as util
import copy

class CondInd():
    def __init__(self,X=set(),Y=set(),Z=set()):
        self.X = copy.copy(X)
        self.Y = copy.copy(Y)
        self.Z = copy.copy(Z)

    def _repr_latex_(self):
        return '$$' + self.tolatexstr() + '$$'

    def tolatexstr(self):
        """Convert a conditional independence statement to a Latex expression for display. Returns a Latex
           format String."""
        str = util.csepstr(self.X) + ' \perp\!\!\!\perp ' + util.csepstr(self.Y)
        if self.Z:
            str = str + ' \mid ' + util.csepstr(self.Z)
        return '(' + str + ')'

    def tostr(self):
        """Convert a conditional independence statement to a (Unicode) String statement for display."""
        str = util.csepstr(self.X) + '\u22a5' + util.csepstr(self.Y)
        if self.Z:
            str = str + '|' + util.csepstr(self.Z)
        return '(' + str + ')'

class CondIndSet():
    def __init__(self,ciset=set()):
        self.ciset = copy.copy(ciset)

    def _repr_latex_(self):
        return '$$' + self.tolatexstr() + '$$'
    
    def tolatexstr(self):
        """Convert a set of conditional independence statements to a Latex expression for display.
           Returns a Latex format String.
        """
        strset = set()
        for ci in self.ciset:
            strset.add(ci.tolatexstr())
        return util.csepstr(strset)

    def tostr(self):
        """Convert a set of conditional independence statements to a single (Unicode) String statement
           for display.
        """
        strset = set()
        for ci in self.ciset:
            strset.add(ci.tostr())
        return util.csepstr(strset)
