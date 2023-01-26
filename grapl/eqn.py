# -*- coding: utf-8 -*-
"""
Algorithms for representing equations as pairs of LHS, RHS Expr expression objects, used to
represent various distributions on ADMGs.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.expr as expr

class Eqn():
    def __init__(self,lhs=expr.Expr(),rhs=expr.Expr()):
        self.lhs = lhs
        self.rhs = rhs

    def _repr_latex_(self):
        return '$$' + self.tocondstr() + '$$'

    def tostr(self):
        """Convert an equation to a Latex expression for display. Returns a Latex
           format String."""
        return self.lhs.tostr() + '=' + self.rhs.tostr()

    def tocondstr(self):
        """Convert an equation to conditional form, and then into Latex format. Returns a Latex
           format String.
        """
        return self.lhs.tocondstr() + '=' + self.rhs.tocondstr()
