title: libdirectional Tutorial
author: Amir Vaxman
date: 15 June 2017
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libdirectional tutorial notes


# Table of contents

* [Chapter 1: Discretization and Representation](#chapter1:discandrep)
    * [Discretization](#discretization)
    * [Representation](#representation)
    * [101 Principal Matching](#principalmatching)
    * [102 Sampling](#sampling)
    * [103 Combing](#combing)
* [Chapter 2: Cartesian Representations](#chapter2:cartesian)
    * [Cartesian Fields](#cartesian)
    * [201 Globally Optimal Fields](#globallyoptimal)
    * [202 PolyVectors](#polyvectors)
    * [203 Integrable PolyVectors](#integrablepolyvectors)
    * [204 Conjugate Fields](#conjugatefields)
* [Chapter 3: Polar Representations](#chapter3:polar)
    * [301 Trivial Connection](#trivialconnection)
* [Chapter 4: Visualization](#chapter4:visualization)
    * [401 Streamline Tracing](#streamlinetracing)
* [Outlook for continuing development](#future)

# Chapter 1 [#chapter1:discandrep]

## [Discretization](#discretization)[discretization]


# Outlook for continuing development [future]

libdirectional is a budding project, and there are many algorithms that we look forward to implement. Prominent examples are:

1. Face-based polar representation, and mixed-integer directional algorithms.

2. Support for 3D *Octahedral* fields, both in tet meshes and with the boundary-element method.

3. Discrete exterior calculus.

4. Differential operators and Hodge decomposition.

5. Cutting, integration, and parameterization.

6. Support for tensor fields.

[#bommes_2009]: David Bommes, Henrik Zimmer, Leif Kobbelt.
  [Mixed-integer
  quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf),
  2009.
[#bouaziz_2012]: Sofien Bouaziz, Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly
  [Shape-Up: Shaping Discrete Geometry with
  Projections](http://lgg.epfl.ch/publications/2012/shapeup.pdf), 2012
[#diamanti_2014]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga
  Sorkine-Hornung. [Designing N-PolyVector Fields with Complex
  Polynomials](http://igl.ethz.ch/projects/complex-roots/), 2014
[#diamanti_2015]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga
  Sorkine-Hornung. [Integrable PolyVector Fields](http://igl.ethz.ch/projects/integrable/), 2015
[#knoppel_2013]: Felix Knöppel, Keenan Crane, Ulrich Pinkall, and Peter
  Schröder. [Globally Optimal Direction
  Fields](http://www.cs.columbia.edu/~keenan/Projects/GloballyOptimalDirectionFields/paper.pdf),
  2013.
[#levy_2008]: Nicolas Ray, Bruno Vallet, Wan Chiu Li, Bruno Lévy.
  [N-Symmetry Direction Field
  Design](http://alice.loria.fr/publications/papers/2008/DGF/NSDFD-TOG.pdf),
  2008.
[#liu_2011]: Yang Liu, Weiwei Xu, Jun Wang, Lifeng Zhu, Baining Guo, Falai Chen, Guoping
  Wang.  [General Planar Quadrilateral Mesh Design Using Conjugate Direction
  Field](http://research.microsoft.com/en-us/um/people/yangliu/publication/cdf.pdf),
  2008.
[#panozzo_2014]: Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga
  Sorkine-Hornung.  [Frame Fields: Anisotropic and Non-Orthogonal Cross
  Fields](http://cs.nyu.edu/~panozzo/papers/frame-fields-2014.pdf),
  2014.
[#vaxman_2016]: Amir Vaxman, Marcel Campen, Olga Diamanti, Daniele Panozzo,
  David Bommes, Klaus Hildebrandt, Mirela Ben-Chen. [Directional Field
  Synthesis, Design, and
  Processing](https://www.google.com/search?q=Directional+Field+Synthesis+Design+and+Processing),
  2016

