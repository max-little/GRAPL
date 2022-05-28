# -*- coding: utf-8 -*-
"""
Algorithms for representing and manipulating ratios of distribution
functions, used to construct interventional distributions on ADMGs.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.util as util


class Expr():
    """A distribution expression object."""
    def __init__(self,num=[],den=[],mrg=set()):
        self.num = list(num)
        self.den = list(den)
        self.mrg = set(mrg)

    def addvars(self, num=[], den=[], mrg=set()):
        """Multiply by a set of numerator variables, divide by a set of denominator
           variables, and integrate over a set of variables to marginalize out."""
        if len(num) > 0:
            self.num.append(num)
        if len(den) > 0:
            self.den.append(den)
        if len(mrg) > 0:
            self.mrg.add(mrg)
    
    def subsvar(self, old, new):
        """Substitute the old name for a variable, with a new name."""
        self.num = [{new if (v == old) else v for v in term} for term in self.num]
        self.den = [{new if (v == old) else v for v in term} for term in self.den]
        self.mrg = {new if (v == old) else v for v in self.mrg}
    
    def fix(self, fix, blanket):
        """Apply the interventional fixing operation to a distribution expression."""
        if len(blanket) == 0:
            self.addvars(den={fix})
        else:
            self.addvars(num=blanket, den=blanket.union({fix}))

    def fixmarginal(self, fix):
        """Apply the interventional fixing operation where it can be marginalized."""
        new_var = fix + chr(39)     # Avoid variable name clashes
        self.addvars(mrg = new_var)
        self.subsvar(fix, new_var)

    def cancel(self):
        """An algorithm for cancelling variables in a distribution."""
        num = self.num.copy()
        den = self.den.copy()
        changes = False
        for term in self.num:
            if term in den:
                den.remove(term)
                num.remove(term)
                changes = True
        if changes:
            self.num = num.copy()
            self.den = den.copy()
        return changes

    def marginal(self):
        """An algorithm for marginalizing out variables in a distribution."""
        num = self.num.copy()
        mrg = self.mrg.copy()
        changes = False
        for node in self.mrg:
            num_cnt = sum([sum({1 if (v == node) else 0 for v in term}) for term in self.num])
            den_find = any([any({(v == node) for v in term}) for term in self.den])
            if (num_cnt == 1) and not den_find:
                changes = True
                num = [term.difference({node}) for term in num]
                mrg.remove(node)
        if changes:
            self.num = num.copy()
            self.mrg = mrg.copy()
        return changes

    def simplify(self):
        """An algorithm for simplifying a distribution, by successive cancellation and marginalization
           until a fixed point is reached."""
        changes = True
        simplified = False
        while changes:
            change1 = self.cancel()
            change2 = self.marginal()
            changes = change1 or change2
            simplified = simplified or changes
        return simplified

    def tostr(self):
        """Convert a distribution to a Latex expression for display."""
        expr_str = ''
        if len(self.mrg) > 0:
            expr_str = expr_str + util.mrgstr(self.mrg)
        if len(self.num) > 1:
            expr_str = expr_str + '['
        for term in self.num:
            expr_str = expr_str + util.probstr(term)
        if len(self.den) > 0:
            expr_str = expr_str + '/'
            if len(self.den) > 1:
                expr_str = expr_str + '{'
            for term in self.den:
                expr_str = expr_str + util.probstr(term)
            if len(self.den) > 1:
                expr_str = expr_str + '}'
        if (len(self.num) > 1):
            expr_str = expr_str + ']'
        return expr_str

    def tocondstr(self):
        """Convert a distribution to conditionals, and then into Latex format."""
        expr_str = ''
        if len(self.mrg) > 0:
            expr_str = expr_str + util.mrgstr(self.mrg)
        if len(self.num) > 1:
            expr_str = expr_str + '['
        num_ord = sorted(self.num,key=len,reverse=True)
        den_ord = sorted(self.den,key=len,reverse=True)
        while len(num_ord) > 0:
            next_num = num_ord[0]
            num_ord = num_ord[1:]
            assign_term = False
            den_list = den_ord.copy()
            while (len(den_list) > 0) and (not assign_term):
                next_den = den_list[0]
                den_list = den_list[1:]
                if next_den.issubset(next_num):
                    den_ord.remove(next_den)
                    assign_term = True
                    expr_str = expr_str + util.probcndstr(next_num.difference(next_den), next_den)
            if not assign_term:
                expr_str = expr_str + util.probstr(next_num)
        if len(den_ord) > 0:
            expr_str = expr_str + '/'
            for den in den_ord:
                expr_str = expr_str + util.probstr(den)
        if (len(self.num) > 1):
            expr_str = expr_str + ']'
        return expr_str
