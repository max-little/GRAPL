# -*- coding: utf-8 -*-
"""
Implements a lexer and parser of a simple domain-specific language
(GRAPL -- GRAph Processing Language) for describing acyclic directed mixed graphs.
A quasi-EBNF description of the language goes as follows, where 'grapl' is the top
level rule:

        grapl = TITLE EOC graph | graph
        graph = graph command | command
        command = NODE EOC | NODE DIR_EDGE NODE EOC | NODE BI_EDGE NODE EOC

The terminal identifiers are:

        TITLE: any string of characters enclosed in quotes
        NODE: any non-empty string of alphabetical characters including underline,
              followed optionally by any natural number
        EOC: ';', the end of command character
        DIR_EDGE: '->', indicating a directed edge
        BI_EDGE: '<->', indicating a bidirected edge

It is recommended not to call functions without docstrings in this module, externally.

(CC BY-SA 4.0) 2021. If you use this code, please cite: M.A. Little, R. Badawy,
2019, "Causal bootstrapping", arXiv:1910.09648
"""

import grapl.ply.lex as lex
import grapl.ply.yacc as yacc
import grapl.admg as admg
import copy


class GraplDSL():
    """A GRAPL lexer and parser object, implemented using PLY. See PLY documentation
       in this package, for further details.
    """
    def __init__(self):
        self.lexer = lex.lex(module=self)
        self.parser = yacc.yacc(module=self)

    tokens = (
        'NODE',
        'DIR_EDGE',
        'BI_EDGE',
        'EOC',
        'TITLE'
    )

    t_DIR_EDGE = '->'
    t_BI_EDGE = '<->'
    t_EOC = ';'
    t_ignore = ' \t'

    def t_TITLE(self,t):
        r'\"(.)*\"'
        return t

    def t_NODE(self,t):
        r'([a-zA-Z_])+(\d)*'
        return t

    def t_newline(self,t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    def t_error(self,t):
        raise SyntaxError("Grapl lexical error at line: " + str(t.lineno) + ", near: " + t.value)

    def p_grapl(self,p):
        '''grapl : TITLE EOC graph
                 | graph'''
        if len(p) == 4:
            self.admg.settitle(p[1])

    def p_graph(self,p):
        '''graph : graph command
                | command'''
        pass

    def p_command(self,p):
        '''command : NODE EOC
                | NODE DIR_EDGE NODE EOC
                | NODE BI_EDGE NODE EOC'''
        if len(p) == 3:
            self.admg.addvar(p[1])
        elif p[2] == self.t_DIR_EDGE:
            if p[1] not in self.admg.vars:
                print("Grapl error: node " + p[1] + " not already defined at line: " + str(p.lineno(1)))
            elif p[3] not in self.admg.vars:
                print("Grapl error: node " + p[3] + " not already defined at line: " + str(p.lineno(3)))
            else:
                self.admg.addedges(p[3],parents={p[1]})
        elif p[2] == self.t_BI_EDGE:
            if p[1] not in self.admg.vars:
                print("Grapl error: node " + p[1] + " not already defined at line: " + str(p.lineno(1)))
            elif p[3] not in self.admg.vars:
                print("Grapl error: node " + p[3] + " not already defined at line: " + str(p.lineno(3)))
            else:
                self.admg.addedges(p[3],bidirects={p[1]})
        
    def p_error(self,p):
        if p is not None:
            raise SyntaxError("Grapl syntax error at line: " + str(p.lineno) + ", near: " + p.value)
        else:
            raise SyntaxError("Grapl syntax error: unexpected end of input")

    def readgrapl(self, scm_string=''):
        """Read a GRAPL syntax string describing an ADMG, and create a graph from it. Returns
           an ADMG object; see the description of the ADMG class, for further details. Raises
           a SyntaxError exception if the GRAPL string is not well formed.
        """
        self.admg = admg.ADMG()
        try:
            self.parser.parse(scm_string)
        except SyntaxError as ex:
            print(ex)
            return None
        self.admg.connect()
        return copy.deepcopy(self.admg)
    
    def savegrapl(self, admg):
        """Given an ADMG object as parameter, outputs a corresponding GRAPL syntax string.
           See the description of the ADMG class for further details about this structure.
           Returns a String.
        """
        nodes = set(admg.vars)
        scm_string = admg.title + self.t_EOC + ' '

        for node in nodes:
            scm_string = scm_string + node + self.t_EOC + ' '
        
        for node in nodes:
            parents = admg.pa({node})
            for parent in parents:
                scm_string = scm_string + parent + self.t_DIR_EDGE + node + self.t_EOC + ' '
        
        for node in nodes:
            bidirects = admg.bi({node})
            for bidirect in bidirects:
                scm_string = scm_string + node + self.t_BI_EDGE + bidirect + self.t_EOC + ' '

        return scm_string

    def displaygrapltokens(self, scm_string=''):
        """For debugging, print out the sequence of lexical tokens in a GRAPL syntax string."""
        self.admg = admg.ADMG()
        self.lexer.input(scm_string)
        while True:
            tok = self.lexer.token()
            if not tok: 
                break
            print(tok)
