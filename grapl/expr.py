# -*- coding: utf-8 -*-
"""
Algorithms for representing and manipulating ratios of distribution
functions, used to construct interventional distributions on ADMGs.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.util as util


class Expr():
    """A distribution expression object, as a ratio of products of marginalized distributions.
       The member variables are: num (numerator), a List of Sets of Strings, each set being
       a distribution, consisting of a set of random variables (node names); den (denominator)
       same as with num; and mrg (marginals), a Set of Strings, each string being the name of a
       random variable which is marginalized out over the whole ratio distribution.
    """
    def __init__(self,num=[],den=[],mrg=set()):
        self.num = list(num)
        self.den = list(den)
        self.mrg = set(mrg)

    def addvars(self, num=[], den=[], mrg=set()):
        """Adds a new set of numerator, denominator and marginal variables to an expression object.
           The parameters are: num (numerator), a List of Sets of Strings, each set being
           a distribution, consisting of a set of random variables (node names); den (denominator)
           same as with num; and mrg (marginals), a Set of Strings, each string being the name of a
           random variable which is marginalized out.
        """
        if len(num) > 0:
            self.num.append(num)
        if len(den) > 0:
            self.den.append(den)
        if len(mrg) > 0:
            self.mrg.add(mrg)
    
    def subsvar(self, old, new):
        """Substitute the old name for a random variable, with a new name. The parameters are:
           old (old variable) String name for an existing random variable; and new: String
           name to which the old variable will be renamed.
        """
        self.num = [{new if (v == old) else v for v in term} for term in self.num]
        self.den = [{new if (v == old) else v for v in term} for term in self.den]
        self.mrg = {new if (v == old) else v for v in self.mrg}
    
    def fix(self, fix, blanket):
        """Apply the interventional fixing operation to a distribution expression. Parameters
           are fix: a String naming the random variable to fix; and blanket: a Set of Strings,
           giving the effective parents of the variable in fix.
        """
        if len(blanket) == 0:
            self.addvars(den={fix})
        else:
            self.addvars(num=blanket, den=blanket.union({fix}))

    def fixmarginal(self, fix):
        """Apply the interventional fixing operation where it can be marginalized. The fix
           parameter is a string giving the name of the random variable to fix. To avoid
           naming clashes, a prime character is appended to the fix variable, and any
           random variables with the name in fix, are renamed accordingly.
        """
        new_var = fix + chr(39)     # Avoid variable name clashes
        self.addvars(mrg = new_var)
        self.subsvar(fix, new_var)

    def cancel(self):
        """An algorithm for cancelling variables in a distribution expression. This seeks to greedily
           match and remove terms appearing in both numerator and denominator of an expression.
           Returns True if any changes to the expression occurred as a result, and False
           otherwise.
        """
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
        """An algorithm for marginalizing out variables in a distribution expression. Greedily
           removes variables appearing in both the numerator and the set of marginal variables.
           Returns True if any changes to the expression occurred as a result, and False
           otherwise.
        """
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
           until a fixed point is reached. Returns True if simplifications were possible, and
           False otherwise."""
        changes = True
        simplified = False
        while changes:
            change1 = self.cancel()
            change2 = self.marginal()
            changes = change1 or change2
            simplified = simplified or changes
        return simplified

    def tostr(self):
        """Convert a distribution expression to a Latex expression for display. Returns a Latex
           format String."""
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
        """Convert a distribution expression to conditional form, and then into Latex format. This
           greedily matches terms in the numerator with corresponding terms in the denominator,
           so as to represent ratios of terms as conditionals in the numerator. Returns a Latex
           format String.
        """
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
