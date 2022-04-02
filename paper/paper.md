---
title: 'GRAPL: A computational library for nonparametric structural causal modelling, analysis and inference'
tags:
  - Python
  - structural causal models
  - causal inference
  - DAGs
  - ADMGs
authors:
  - name: Max A. Little
    orcid: 0000-0002-1507-3822
    affiliation: "1, 2"
affiliations:
 - name: School of Computer Science, University of Birmingham, UK
   index: 1
 - name: MIT, Cambridge, MA, USA
   index: 2
date: 28 March 2022
bibliography: paper.bib

---

# Summary

*Structural causal models* (SCMs) provide a probabilistic language for describing directed
relationships between random variables. SCMs are widely used in science, engineering
and statistical modelling to capture causal relationships between quantitative, measured
phenomena in the real world. Two SCM formalisms, *directed acyclic graphs* (DAGs) and
*acyclic directed mixed graphs* (ADMGs) have been extensively studied. In these formalisms,
the conditions under which causal dependence between variables occurs is well understood.
Furthermore, analytical techniques have been developed which allow manipulation of the model
so as to perform *nonparametric causal adjustment*, that is, the isolation of desired causal
relationships from the SCM. The `GRAPL` library described in this paper brings together the
most important and useful of such algorithms in one convenient Python package. Using this
library it is possible to represent, analyze and manipulate DAGs and ADMGs of arbitrary complexity.

# Statement of need

There exist a large number of techniques for statistical estimation of *causal* rather than merely
*associational relationships* from data [@Hernan:2020]. Under special causal relationships such as the existence of
*observed confounders*, classical *causal inference* methods such as stratification, propensity scoring and IPW are widely
used to *adjust* for these spurious associations. However, there are countless other physical scenarios where e.g. the
confounders are *unobserved* or there are multiple, more intricate, cause-effect relationships between measured variables.
In these circumstances, the classical techniques are not valid and more generalized *adjustment* methods are required [@Pearl:2009].

Where all measured variables are observed (so that the SCM is a DAG), it can be shown that it is always
possible to compute the adjusted nonparametric cause-effect (*interventional*) distribution which can then
be used to derive a suitable adjustment technique [@Bareinboim:2020]. Doing these computations by hand for complex DAGs can be
laborious and error-prone; this package automates the construction of such distributions, additionally enabling the
determination of important DAG properties such as all *local Markov conditional dependence* relationships between
variables [@Koller:2009]. These relationships encode for meaningful statistical dependence relationships which can, for instance,
be tested against real data to validate any proposed DAG model.

When there are non-observed variables and the SCM is e.g. an ADMG (such as the existence of *unmeasured common causes*
of pairs of variables) then the above guarantee no longer applies and in some cases it will not be possible,
even in principle, to compute an interventional distribution in order to adjust for these spurious associations [@Richardson:2012].
The computation of exactly when adjustment is possible, and the derivation of a nonparametric, interventional
adjustment distribution when it is possible, is addressed by several recently developed, generic algorithms
[@Shpitser:2008; @Richardson:2012]. These computations are complex and too tedious and/or difficult to apply by hand, and require automation to be practical. However, the relevant algorithms are fairly complex and (to the author's knowledge) there are no publicly-available implementations written in widely used languages such as Python which are accessible to non-statisticians.
(At the time of writing, the one possible exception to this, the R package `causaleffect`, does implement one of these algorithms,
but does not provide the constituent, fine-grained analytical methods required to understand the topological causal structure of a given SCM).

The design of `GRAPL` was intended to bring together all the most general algorithms for handling all associational
and causal distribution computations in DAGs and ADMGs, of arbitrary complexity, into one convenient package.
The library is structured so as to expose all layers of computations such as the necessary fine-grained analysis of
SCM graph topology and symbolic computations with arbitrary nonparametric marginalized conditional distributions on
those graphs. Distributions are output in Latex format for convenient inclusion in publications, and SCMs can be
represented using the built-in domain-specific language for specifying directed relationships between variables.
It is chiefly aimed at researchers wanting to carry out causal inference analysis for a wide array of
disciplines including machine learning, AI, data science, mathematical statistics and also in practical application
areas such as bioinformatics or medical statistics, but can be used by any interested researcher wanting to understand
the full topological, causal structural implications of specific SCMs, for both theoretical and practical reasons.

It is already in regular use as a tool for developing specialized causal inference methods, in particular
*causal bootstrapping* techniques for general machine learning applications [@Little:2019], and for constructing
and analyzing SCMs representing the complex causal processes behind neurodegeneration [@Sturchio:2020].

# Acknowledgements

This work partially funded by NIH grant UR-Udall Center, award number P50 NS108676.

# References