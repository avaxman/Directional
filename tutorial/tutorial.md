title: libdirectional Tutorial
author: Amir Vaxman
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libdirectional Tutorial Notes


# Table of contents

* [Introduction](#chapter0:introduction)
* [Chapter 1: I/O and Visualization](#chapter1:iovis)
    * [101 Glyph rendering](#glyphrendering)
    * [102 Picking and editing](#pickingediting)
    * [103 Streamline tracing](#streamlinetracing)
* [Chapter 2: Discretization and representation](#chapter2:discandrep)
    * [Discretization](#discretization)
    * [Representation](#representation)
    * [201 Principal matching](#principalmatching)
    * [202 Sampling](#sampling)
    * [203 Combing](#combing)
* [Chapter 3: Cartesian methods](#chapter3:cartesian)
    * [Cartesian fields](#cartesian)
    * [301 Power fields](#powerfields)
    * [302 PolyVectors](#polyvectors)
    * [303 Integrable PolyVectors](#integrablePVs)
    * [304 Conjugate Fields](#conjugatefields)
* [Chapter 4: Polar Methods](#chapter4:polar)
    * [Polar Fields](#polar)
    * [401 Index Prescription](#indexprescription)
* [Outlook for continuing development](#future)
* [References](#references)


# <a name="chapter0:introduction">Introduction</a>
libdirectional is a C++ geometry processing library written as an extension to [libigl](http://libigl.github.io/libigl/), with a  speciality in directional fields. The functionality is based on the definitions and taxonomy surveyed theoretically in [#vaxman_2016], and through it by much of the relevant papers in the literature. It contains tools to edit, analyze, and visualize directional fields of various degrees and symmetries. 

The underlying structure extends the general philosophy of [libigl](http://libigl.github.io/libigl/): the library is header only, where each header contains a set (often only one) of functions closely related (for instance, the precomputation and computation of some directional quantity over a mesh). The data structures are, for the most part, simple matrices in [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), and the library avoids complicated and nested structures, instead directly working with calling to standalone functions. The visualization is done using the libigl viewer with some extended options that allow the rendering of directional fields.

The header files contain documentation of the parameters to each function and their required composition; in this tutorial we will mostly tie the functionality to the theoretical concepts of directional fields and their capabilities.

###Installing the tutorial examples

This tutorial comprises an exhaustive set of examples that demonstrate the capabilities of libdirectional, where every subchapter relates to a single example. The tutorial code can be installed by going into the `Tutorial` folder, and following the following instructions:


```cpp
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

This will build all tutorial examples in the `build` folder. To build in windows, use the `cmake-gui .` options instead of the last two, and creates the project using Visual Studio, with the proper tutorial chapter as the "startup project". 

To access a single example, say ``202_Sampling``, go to the ``build`` subfolder, and the executable will be there. Command-line arguments are never required; the data is read from the ``shared`` folder directly for each example. 

Most examples contain a component of user interaction; the instructions of what to do will be given in the command-line output upon execution.

# <a name="chapter1:iovis">Chapter 1: I/O and Visualization</a>

## <a name="glyphrendering">101 Glyph Rendering</a>

The most basic operation on directional fields is reading them from a file and drawing them in the most explicit way. In [Example 101](101_GlyphRendering/main.cpp), a field is read from a file as follows:

```cpp
directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singPositions, singIndices);
```

The field is read in <i>raw</i> format, which is detailed in [Chapter 2: Discretization and Representation](#chapter2:discandrep). The field is <i>face-based</i>, and the singularities are consequently <i>vertex-based</i>.

The field is drawn on the mesh as follows:

```cpp
if (drawSingularities)
    directional::singularity_spheres(V, F, singPositions, singIndices, positiveIndexColors, negativeIndexColors, false, true, fullV, fullF, fullC);
directional::glyph_lines_raw(V, F, rawField, rawGlyphColor, false, true, fullV, fullF, fullC);
```

These two operations do not produce any active drawing; they create meshes that extend the original geometry, and then get passed to libigl viewer.

``directional::singularity_spheres()`` creates small spheres on vertices, where the size of the sphere is devised automatically (but can be configured using the extended version of this function). The spheres are only created where the index is different than $0$.

``directional::glyph_lines_raw()`` creates lines on the faces that constitute the <i>glyph drawing</i>: simply drawing the vectors upon the faces in their $\left(x,y,z\right)$ coordinates, starting from the face barycenter. There are several ways to give colors to these vectors, which can be individual or global; check the documentation to the function in the header.
  
Vectors are drawn in their given magnitudes, up to a global scale. By default, this scale is set to be related to the average edge length---it can be manually set by the extended version of the function. 

![([Example 101](101_GlyphRendering/main.cpp)) Glyph Rendering on a mesh, with singularities visible.](images/101_GlyphRendering.png) 

## [102 Picking and editing](#pickingediting)[pickingediting]

This is a simple tutorial that demonstrates how libdirectional uses libigl picking to make directional field editing possible. A face and a vector within the face are chosen, and clicking on a new direction for the vector changes it. Note the ability to set different colors for glyphs, via the following code in [Example 102](101_PickingEditing/main.cpp).

```cpp
Eigen::MatrixXd fullGlyphColor(F.rows(),3*N);
  for (int i=0;i<F.rows();i++){
    for (int j=0;j<N;j++){
      if (i==currF)
        fullGlyphColor.block(i,3*j,1,3)<<(j==currVec ? selectedVectorGlyphColor : selectedFaceGlyphColor);
      else
        fullGlyphColor.block(i,3*j,1,3)<<defaultGlyphColor;
    }
  }
```

The size of ``fullGlyphColor`` can either be one color per vector, per face, or uniform for all vetors in the entire mesh; the intent will automatically be devised from the size of the matrix.

![([Example 102](102_PickingEditing/main.cpp)) Editing several vectors on a single face.](images/102_PickingEditing.png)

## [103 Streamline Tracing](#streamlinetracing)[streamlinetracing]

Vector fields on surfaces are commonly visualized by tracing [streamlines] (https://en.wikipedia.org/wiki/Streamlines,_streaklines,_and_pathlines). libdirectional supports the seeding and tracing of streamlines, for all type of directionals. The seeds for the streamlines are initialized using `streamlines_init`, and the lines are traced using `streamlines_next`. Each call to `streamlines_next` extends each line by one triangle, allowing interactive rendering of the traced lines, as demonstrated in [Example 103](103_StreamlineTracing/main.cpp).

![([Example 103](103_StreamlineTracing/main.cpp)) Interactive streamlines tracing.](images/103_StreamlineTracing.png)

# Chapter 2: Discretization and Representation [chapter2:discandrep]

## [discretization](#discretization)[discretization]

There are several ways to discretize directional fields, although not all of them fit for all representations. They roughly divide into face-based, edge-based, and vertex-basex discretization (see [#degoes_2016] for an in-depth analysis). The only discretization currently supported by libdirectional is face-based fields, where the discrete tangent plane is considered as the supporting plane to each (naturally flat) face. We also use a local basis (mostly provide by `igl::local_basis()`) to parameterize each tangent plane. Notions like connection, parallel transport, and consequently smoothness, are measured through dual edges between adjacent faces.

## [Representation](#representation)[Representation]

libdirectional uses several different representations to describe directional fields. We denote the number of faces in the mesh as $|F|$, the set of inner edges (adjacent to two triangles) as $E_I$, and the degree of the field (must be fixed for the entire field) as $N$. The supported  representations are as follows, where the taxonomy is based on that of the directional field course [#vaxman_2016]:

1. **Raw** - A $|F|\times 3N$ double matrix, representing an $1^N$-vector field (a directional with $N$ independent vectors in each face) in the form $X_1, Y_1, Z_1, X_2, Y_2, Z_2, \cdots X_N, Y_N, Z_N$ per face. Vectors are assumed to be ordered in counterclockwise order in most libdirectional functions that process raw fields.
2. **Representative**. A $|F| \times 3$ double matrix that represents a rotationally symmetric $N$-vector field, known as an $N$-RoSy. The single vector is an arbitrary "first" vector in the face, and the rest of the vectors are deduced by rotations of $\frac{2\cdot\pi}{N}$
3. **Rotation Angles**. A $|E_I|$-sized double vector representing the rotation angle between two directions (without magnitude information) on two neighbouring triangles. The rotation represents the deviation from the Levi-Civita parallel transport [#ray_2008], [#crane_2010]. This representation may only encode $N$-direction fields. Note that the <i>effort</i> (sum of all rotations) is then $N$ times rotation angles. Since this is a differential quantity, an extra global rotation needs to be given to uniquely create the full face-based field.
4. **Power Field** - An $|F|$-sized *complex* vector, representing an $N$-RoSy as a single complex number $y=u^N$, where the $N$-RoSy is the set of roots $u$. The magnitude is also encoded this way, though it may be neglected in some applications.
5. **PolyVector** - A $|F| \times N$ complex double matrix, encoding the coefficients $a$ of a complex polynomial $f(z)=\sum_{i=0}^{N-1}{a_iz^i}+z^N$, which roots $u$ are an $1^N$-vector field. Every row is encoded as $a_{0},\cdots, a_{N-1}$, where $a_0$ is the free coefficient. In case where the field is an $N$-RoSy, all coefficients but $a_0$ are zero.

libdirectional provides a number of conversion functions to switch between different representations. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the representation names in the above list. e.g., `rotation_to_representative()` and `polyvector_to_raw()`. Some possible combinations are given by composing two functions in sequence.

However, note that every conversion is possible; for instance, it is not possible to convert between PolyVectors and rotation angles, as they do not possess the same generality of directional fields (with current state-of-the-art...). For $N$-RoSy fields, for instance, you will most likely work primarily with the power field, representative, or rotation-angle representation. converting into the more explicit raw representation is done mostly for I/O and visualization.

In the following sections, we show some effects of working with different representations and converting between them.

## [201 Principal Matching](#principalmatching)[principalmatching]

One of the fundamental operations in directional-field processing is <i>matching</i>. That is, defining which vectors in face $f_i$ correspond to those in adjacent face $f_j$. In libdirectional, we only work with order-preserving matchings: if vector $k$ in face $f_i$ is matched to vector $m$ in face $f_j$, then for any $l$, $k+l$ is matched to $m+l$ (modulu $N$) in the respective faces. Suppose that the orientation of the dual edge is $f_i \rightarrow f_j$, then the matching is encoded as $m-k$. Some representations, like rotation angles, already encode the matching explicitly, but others do not. Therefore, it needs to be devised from the field.

Given a raw field (in assumed CCW order in every face), it is possible to devise the rotation angles $\delta_{ij}$ by the process of *principal matching* [#diamanti_2014]. Principal matching creates the matching that minimizes the effort of the matching, always putting it within the range of $[-\pi, \pi)$ (and therefore denoted as "principal"). It corresponds to the "smallest angle" matching for $N$-RoSy fields.

principal matching is done through this function (from in [Example 201](201_PrincipalMatching/main.cpp)):

```cpp
directional::principal_matching(V, F,EV, EF, FE, rawField, matching, effort);
directional::effort_to_indices(V,F,EV, EF, effort,N,prinIndices);
```
`directional::effort_to_indices()` computes the <i>index</i> of each vertex from the effort around it. The index of a vertex is the amount of rotations a directional undergoes along a cycle around the vertex. a directional must return to itself after a cycle, and therefore the index is an integer $I$ when a vector $m$ in the face ended up in vector $m+I$ (modulu $N$). Note that this can also include multiple full rotations, and therefore the index is unbounded. The <i>fractional</i> part of the index is encoded by the matching; however, matching alone cannot encode <i>integral</i> indices (for instance, a single vector field has trivial (Zero) matching anywhere, but can have singularities).

![([Example 201](201_PrincipalMatching/main.cpp)) A Field is shown with singularities, and a single face is shown with the principal matching to its neighbors (in multiple colors).](images/201_PrincipalMatching.png)


## [202 Sampling](#sampling)[sampling]

This is an educatory example that demonstrates the loss of information when moving between a polar (in this case, rotation angle) representation, to a Cartesian representation, where the matching between vectors in adjacent faces is done with principal matching. There are three modes seen in the example:

1. In the polar mode, the user can control the index of a singularity directly. As such, the rotation angles between faces become arbitrarily large, and appear as noise in the low valence cycles. 

2. In the principal matching mode, the rotations are computed from the field, without prior knowledge of the polar prescribed rotations from the previous mode. The noise of the field then gives rise to a "singularity party". 

3. In the Cartesian mode, the field is interpolated on the white faces, keeping the red band fixed from the polar mode. We see a field that in smooth in the Cartesian sense, with more uniformly dispersed singularities.

![([Example 202](202_Sampling/main.cpp)) Left to right: the polar mode, the principal matching mode, and the Cartesian mode.](images/202_Sampling.png)

## [203 Combing](#combing)[combing]

Given a matching (in this case, principal matching), it is possible to "comb" the field. That is, re-index each face (keeping the CCW order), so that the vector indexing aligns perfectly with the matching to the neighbors---then, the new matching on the dual edges becomes trivially zero. This operation is important in order to preapre a directional field for integration, for instance. In the presence of singularities, the field can only be combed up to a set of connected paths that connect between singularities, also known as <i>cuts</i>. Note that such paths do not cut the mesh to a simply-connected patch, but only connects subgroups with indices adding up to an integer; as a trivial example, a 1-vector field is trivially combed to begin with, even in the presence of its (integral) singularities, and nothing happens. The combing is done through the function `directional::principal_combing()`.

![([Example 203](203_Combing/main.cpp)) Vector indices inside each face are colored. They are uncombed in the left image, and combed, with the cut showing, in the right image.](images/203_Combing.png) 

# Chapter 3: Cartesian Representations [chapter3:cartesian]

## [Cartesian Fields](#cartesian)[cartesian]

The Cartesian representation is a meta-category for representation of vectors in explicit coordinates, either $\left(x,y\right)$ in some local $2D$ basis on a tangent plane, or $\left(x,y,z\right)$ in the ambient coordinates of space. The raw, representative (of an $N$-RoSy), power field, and PolyVector representations are all such examples. Cartesian fields often do not automatically contain information about the matching, or rotation, of a field between one face and the next, and it needs to be computed using principal matching. This chapter focuses on computing fields with this representation. 

## [301 Power Fields](#powerfields)[powerfields]

This representation is offered in [#knoppel_2013], but they did not give it a specific name---we use the name given in [#azencot_2017].

A power field representation uses a complex basis in each tangent plane (face in our implementation), and represents an $N$-RoSy using a <i>power vector</i>---a single complex number $y$ per face so that its root set $y=u^N$ comprise the $N$-RoSy. 

By prescribing constraints $y_B$ on a set of faces $B$, the algorithm interpolates the field to the rest of the faces $y_I$ by minimizing the face-based Dirichlet energy: $$y_I=\text{argmin}\sum_{e=(f,g)\in F \times F}{\left|y_fe_f^N - y_ge_g^N\right|^2},$$
where $e_f$ is the representation of the vector of edge $e$ in the basis of $f$, and similarly for $g$. The field is computed through the function `directional::power_field()`.

For fixed set $B$ and changing $y_B$, It is possible to speed up computations by precomputing the solver (sparse Cholsky for the positive-definite matrix) used to compute the power field. This is done by using the function `directional::power_field_precompute()`, coupled with the appropriate version of `directional::power_field()`. Note that field can be converted to representative and raw forms using the appropriate `power_to_X` functions.

![([Example 301](301_GloballyOptimal/main.cpp)) Setting up a small subset of constraints (red faces), and interpolating (and normalizing in magnitude) the power field to the rest of the mes. Note the singularities that are discovered through principal matching.](images/301_GloballyOptimal.png) 


## [302 PolyVectors](#polyvectors)[polyvectors]

A Polyvector field [#diamanti_2014] is a generalization of power fields that allows to represent independent vectors in each tangent plane. The representation is as the coefficient set $a_{0 \cdots N-1}$ of a complex polynomial in the local compex basis:
$$P(z) = a_0 + a_1z + \ldots + a_{N-1} z^{N-1} + z^N$$
where the roots $P(z)=0$ are the vectors represented as complex numbers. The Dirichlet energy is as before, but there is a term for each $a_i$, with the right power $i$. Note that an $N$-RoSy is represented as a polynomial where all $a$ are zero except $a_0$. Principal matching, combing, and effort are well-defined on PolyVectors as well.

The example allows a user to set individual vectors within each face, and see the interpolated result. The responsible function is `directional::polyvector_field()`. In the case as well, the solver can be prefactored in advance using `directional::polyvector_precompute()`.

![([Example 302](302_PolyVectors/main.cpp)) Vectors are constrained individually in the constrained faces (red), and interpolated to the rest of the faces](images/302_PolyVectors.png) 

## [303 Integrable PolyVectors](#integrablePVs)[integrablePVs]

<span style="color:red">
This tutorial is a direct migration of the libigl version, as the functionality would reside in libdirectional. It is fully functional, but the syntax is not libdirectional-compatible as of yet.
</span>

Vector-field guided surface parameterization is based on the idea of designing the gradients
of the parameterization functions (which are tangent vector fields on the surface) instead of the functions themselves. Thus, vector-set fields (N-Rosy, frame fields, and polyvector fields) that are to be used for parameterization (and subsequent remeshing) need to be integrable: it must be possible to break them down into individual vector fields that are gradients of scalar functions. Fields obtained by most smoothness-based design methods (eg. [#ray_2008][], [#knoppel_2013][], [#diamanti_2014][], [#bommes_2009][], [#panozzo_2014][]) do not have this property. In [#diamanti_2015][], a method for creating integrable polyvector fields was introduced. This method takes as input a given field and improves its integrability by removing the vector field curl, thus turning it into a gradient of a function ([Example 303](303_IntegrablePVs/main.cpp)).

![Integration error is removed from a frame field to produce a field aligned parameterization free of triangle flips.](images/303_IntegrablePVs.png)

This method retains much of the core principles of the polyvector framework - it expresses the condition for zero discrete curl condition (which typically requires integers for the vector matchings) into a condition involving continuous variables only. This is done using coefficients of appropriately defined polynomials. The parameterizations generated by the resulting fields are exactly aligned to the field directions and contain no inverted triangles.

## [304 Conjugate Fields](#conjugatefields)[conjugatefields]

<span style="color:red">
This tutorial is a direct migration of the libigl version, as the functionality would reside in libdirectional. It is fully functional, but the syntax is not libdirectional-compatible as of yet.
</span>

Two tangent vectors lying on a face of a triangle mesh are conjugate if

\\[ k_1 (u^T d_1)(v^T d_1) + k_2(u^T d_2)(v^T d_2) = 0. \\]

This condition is very important in architectural geometry: The faces of an
infinitely dense quad mesh whose edges are aligned with a conjugate field are
planar. Thus, a quad mesh whose edges follow a conjugate field  are easier to
planarize [#liu_2011].

Finding a conjugate vector field that satisfies given directional constraints
is a standard problem in architectural geometry, which can be tackled by
deforming a Poly-Vector field to the closest conjugate field.

This algorithm [#diamanti_2014] alternates a global step, which enforces
smoothness, with a local step, that projects the field on every face to the
closest conjugate field ([Example 304](304_ConjugateField/main.cpp)).

![A smooth 4-PolyVector field (left) is deformed to become a conjugate field
(right).](images/304_ConjugateFields.png)


# Chapter 4: Polar Representations [chapter4:polar]

## [Polar Fields](#polar)[polar]

Polar fields are represented using angles. These angles may encode the rotation from some known basis on a tangent plane (and so it is a "logarithmic" representation, when compared to Cartesian methods), or an angle difference between two neighboring tangent planes (in the sense of deviation from parallel transport). The former usually requires integer variables for directional field design. The latter does not, but state-of-the-art methods require the prescription of indices around independent dual cycles in the mesh. Currently, libdirectional supports the latter.

## [401 Index Prescription](#indexprescription)[indexprescription]

The notation of encoding rotation angles on dual edges, as means to encode deviation from parallel transport between adjacent tangent planes, appeared in several formats in the literature [#ray_2008], [#crane_2010]. The formulation and notation we use in libdirectional is that of Trivial Connections [#crane_2010]. Trivial connection solves for a single rotation angle $\delta_{ij}$ per (dual) edge $e_{ij}$ between two faces $f_i,f_j$, encoding the deviation from parallel transport between them. The algorithm first computes a spanning set of basis cycles (see next section), around which the sum of $\delta_{ij}$ has to be prescribed. The summation is defined as matrix $H$. Every such cycle (row in the matrix) has a curvature, defined as a discrete angle defect, and the prescribed index defines an alternative curvature. The algorithm solves for the smoothest field, in the 2-norm least squares sense, as follows:

$$
\delta = \text{argmin}\ |\delta_{ij}|^2\ s.t.\ H\delta = -K_0 + K.
$$

$H$ is the matrix that defines the basis-cycles sum, $K_0$ is the original curvature of the basis cycle, and $K$ is the prescribed curvature. $K$ defines singularities: for regular cycles, we prescribe $K=0$, and for a singular cycle with prescribed singularity index $\frac{1}{N}$, we set $K=\frac{2\pi}{N}$. the sum of $K$ has to conform to the Poincar&eacute; index theorem, except handle cycles which can have unbounded index. See [#crane_2010] for exact details. If the input obeys the sum, the result obeys the prescribed indices around the cycles, and nothing else. As the representation is differential, there is still a global degree of freedom in setting a single direction in a single arbitrary face.

Note that the correct definition for "cycle curvature" corresponds to the so-called "cycle holonomy", only up to integer multiples of $2\pi$. However, in the discrete setting, the curvature should theoretically be computed as the exact angle defect, in which for inner vertices we use $2\pi-\sum{\alpha}$, and for boundary vertices we use $\pi - \sum{\alpha}$ ($\alpha$ are the angles at the corners of a vertex). For a cycle aggregating many vertices, such as a boundary loop cycle, we add up all the defects. That is requires for exact discrete Poincar&eacute; index consistency.

### Basis Cycles

The basis cycles form the cycles around which curvatures (and singularities) are prescribed on the mesh. The sum on basis cycles is described in a sparse matrix $H$ of size $|cycles|\times |E_I|$, where $E_I$ is the number of inner edges in the mesh. Each row in the matrix describes the sum over one cycle, and contains 1 or -1 values depending on the (arbitrary) orientation of the dual edge participating in the cycle to the respective face. There are three types of cycles, so ordered in the rows of $H$: 

1. $1$-ring dual cycles around each inner vertex, on which vertex-based singularities can be encoded (the relevant part of $H$ is basically $\left(d_0\right)^T$ in discrete exterior calculus, restricted to inner edges).
2. Cycles around mesh boundary loops. 
3. Cycles around topological generators (independent handles).

The method `directional::dual_cycles()` computes the proper basis cycles and matrix $H$. To be able to intuitively prescribe singularities to inner vertices, the method also returns a conversion vector ``vertex2cycle``, and the list of indices of inner edges from the list of edges.

The singularity indices that are prescribed contain the singularity index corresponding to each basis cycle. A value of $k \in \mathbb{Z}$ represents an $\frac{2\pi k}{N}$ rotation around the respective cycle.

If the prescribed indices do not conform to correct sum, a result will still be computed by least squares, but it will be unpredictable.

The algorithm is performed through the function ``directional::index_prescription()``, which also accepts a solver for precomputation.


![([Example 401](401_IndexPrescription/main.cpp)) Indices are prescribed on three vertex singularities, and on the boundary loop, to match the index theorem. The computed field is smooth and obeys these indices exactly.](images/401_IndexPrescription.png) 


# Outlook for continuing development [future]

libdirectional is a budding project, and there are many algorithms in the state-of-the-art that we look forward to implement, with the help of volunteer researchers and practitioners from the field. Prominent examples of desired implementations are:

1. Face-based polar representation, and mixed-integer directional algorithms.

2. Support for 3D *Octahedral* fields [#solomon_2017], both in tet meshes and with the boundary-element method.

3. A discrete exterior calculus framework.

4. Differential operators and Hodge decomposition.

5. Cutting, integration, and parameterization. Note the libigl has this capacity that could be called from libdirectional, but they are not entirely compatible.

6. Support for tensor fields.

7. Advanced and better visualization techniques.

# References [references]
[#azencot_2017]: Omri Azencot, Etienne Corman, Mirela Ben-Chen, Maks Ovsjanikov, [Consistent Functional Cross Field Design for Mesh Quadrangulation](http://www.cs.technion.ac.il/~mirela/publications/cfc.pdf), 2017.
[#bommes_2009]: David Bommes, Henrik Zimmer, Leif Kobbelt, [Mixed-integer quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf), 2009.
[#bouaziz_2012]: Sofien Bouaziz, Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly, [Shape-Up: Shaping Discrete Geometry with Projections](http://lgg.epfl.ch/publications/2012/shapeup.pdf), 2012.
[#crane_2010]: Keenan Crane, Mathieu Desbrun, Peter Schr&ouml;der, [Trivial Connections on Discrete Surfaces](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/), 2010.
[#diamanti_2014]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung, [Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/), 2014.
[#diamanti_2015]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung, [Integrable PolyVector Fields](http://igl.ethz.ch/projects/integrable/), 2015.
[#degoes_2016]: Fernando de Goes, Mathieu Desbrun, Yiying Tong, [Vector Field Processing on Triangle Meshes](http://geometry.caltech.edu/pubs/dGDT16.pdf), 2016.
[#knoppel_2013]: Felix Knöppel, Keenan Crane, Ulrich Pinkall, and Peter Schr&ouml;der, [Globally Optimal Direction Fields](http://www.cs.columbia.edu/~keenan/Projects/GloballyOptimalDirectionFields/paper.pdf), 2013.
[#ray_2008]: Nicolas Ray, Bruno Vallet, Wan Chiu Li, Bruno Lévy, [N-Symmetry Direction Field Design](http://alice.loria.fr/publications/papers/2008/DGF/NSDFD-TOG.pdf), 2008.
[#liu_2011]: Yang Liu, Weiwei Xu, Jun Wang, Lifeng Zhu, Baining Guo, Falai Chen, Guoping Wang, [General Planar Quadrilateral Mesh Design Using Conjugate Direction Field](http://research.microsoft.com/en-us/um/people/yangliu/publication/cdf.pdf), 2008.
[#solomon_2017]: Justin Solomon, Amir Vaxman, David Bommes, [Boundary Element Octahedral Fields in Volumes](http://www.staff.science.uu.nl/~vaxma001/frames3d.pdf), 2017.
[#panozzo_2014]: Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung,  [Frame Fields: Anisotropic and Non-Orthogonal Cross Fields](http://cs.nyu.edu/~panozzo/papers/frame-fields-2014.pdf), 2014.
[#vaxman_2016]: Amir Vaxman, Marcel Campen, Olga Diamanti, Daniele Panozzo, David Bommes, Klaus Hildebrandt, Mirela Ben-Chen, [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis), 2016.

