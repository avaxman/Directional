title: libdirectional Tutorial
author: Amir Vaxman
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libdirectional Tutorial Notes


# Table of contents

* [Introduction] 
(#chapter0:introduction)
* [Chapter 1: I/O and Visualization](#chapter1:iovis)
    * [101 Basic Glyph Rendering](#glyphrendering)
    * [102 Picking and editing](#pickingediting)
    * [103 Streamline Tracing](#streamlinetracing)
* [Chapter 2: Discretization and Representation](#chapter2:discandrep)
    * [Discretization](#discretization)
    * [Representation](#representation)
    * [201 Principal Matching](#principalmatching)
    * [202 Sampling](#sampling)
    * [203 Combing](#combing)
* [Chapter 3: Cartesian Methods](#chapter3:cartesian)
    * [Cartesian Fields](#cartesian)
    * [301 Globally Optimal Fields](#globallyoptimal)
    * [302 PolyVectors](#polyvectors)
    * [303 Integrable PolyVectors](#integrablePVs)
    * [304 Conjugate Fields](#conjugatefields)
* [Chapter 4: Polar Methods](#chapter4:polar)
    * [Polar Fields](#polar)
    * [401 Index Prescription](#indexprescription)
* [Outlook for continuing development](#future)


# Introduction [chapter0:introduction]
libdirectional is a C++ geometry processing library written as an extension to [libigl](http://libigl.github.io/libigl/), with a  speciality in directional fields. The functionality is based on the definitions and arrangement surveyed theoretically in [#vaxman_2016], and through it by much of the relevant papers in the literature. It contains tools to edit, analyze, and visualize directional fields of various degrees and symmetries. 

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

To access a single example, say ``202_Sampling``, go to the ``build`` and the executable will be there. Command-line arguments are never required; the data is read from the ``shared`` folder automatically in each example. 

Most examples contain a component of user interaction; the instruction on what to do will be given in the commandline output immediately after executing the example.

# Chapter 1: I/O and Visualization [chapter1:iovis]

## [101 Glyph Rendering](#glyphrendering)[glyphrendering]

The most basic operation on directional fields is reading them from a file and drawing them in the most explicit way. In [Example 101](101_GlyphRendering/main.cpp), a field is read from a file as follows:

```cpp
directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singPositions, singIndices);
```

The field is read in <i>raw</i> format, which is detailed in [Chapter 2: Discretization and Representation](#chapter2:discandrep). The field is <i>face-based</i>, and the singularities are then <i>vertex-based</i>. 

The field is then drawn on the mesh as follows:

```cpp
if (drawSingularities)
    directional::singularity_spheres(V, F, singPositions, singIndices, positiveIndexColors, negativeIndexColors, false, true, fullV, fullF, fullC);
  
  directional::glyph_lines_raw(V, F, rawField, rawGlyphColor, false, true, fullV, fullF, fullC);
```

These two operations in fact do not produce any drawing; they create meshes that extend the original geometry, and then get passed to libigl viewer. 

``directional::singularity_spheres`` creates small spheres on vertices, where the size of the sphere is devised automatically (but can be configured using the extended version of this function. The spheres are only created where the index is different than $0$.

``directional::glyph_lines_raw`` creates lines on the faces that constitute the simple <i>glyph drawing</i>: simply drawing the vectors upon the faces as they are. There are several ways to give colors to these vectors, which can be individual or global; check the documentation to the function in the header.
  
By default, the size of each vector is set to be related to the average edge length, keeping the ratios between the lengths of the actual vectors intact. The base length and with can be manually set by the extended version of the function. 

![([Example 101](101_GlyphRendering/main.cpp)) Glyph Rendering on a mesh.](images/101_GlyphRendering.png)

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

The size of ``fullGlyphColor`` can either be one color per vertex, per face, or per the entire mesh, and the intent will be automatically devised from the size.

![([Example 102](102_PickingEditing/main.cpp)) Editing several vectors on a single face.](images/102_PickingEditing.png)
## [103 Streamline Tracing](#streamlinetracing)[streamlinetracing]

Vector fields on surfaces are commonly visualized by tracing [streamlines] (https://en.wikipedia.org/wiki/Streamlines,_streaklines,_and_pathlines). libdirectional supports the seeding and tracing of streamlines, for all type of directionals. The seeds for the streamlines are initialized using `streamlines_init`, and the lines are traced using `streamlines_next`. Each call to `streamlines_next` extends each line by one triangle, allowing interactive rendering of the traced lines, as demonstrated in [Example 103](103_StreamlineTracing/main.cpp).

![([Example 103](103_StreamlineTracing/main.cpp)) Interactive streamlines tracing.](images/103_StreamlineTracing.png)



# Chapter 2: Discretization and Representation [chapter2:discandrep]

## [Representation](#representation)[Representation]

The only discretization currently supported by libdirectional is face-based fields, where the discrete tangent plane is considered as the supporting plane to each (naturally flat) face. Nevertheless, libdirectional uses several different representations to describe directional fields. We denote $|F|$ as the number of faces in the mesh, $E_I$ as the set of inner edges (adjacent to two triangles), and $N$ as the degree of the field (must be fixed for the entire field). The supported  representations are as follows, where the taxonomy is based on that of the directional field course [#vaxman_2016]:

1. **Raw** - A $|F|\times 3N$ double matrix, representing an $1^N$-vector field (a directional with $N$ independent vectors in each face) in the form $X_1, Y_1, Z_1, X_2, Y_2, Z_2, \cdots X_N, Y_N, Z_N$ per face. Vectors are assumed to be ordered in counter clockwise order in most libdirectional functions that process raw fields.
2. **Representative**. A $|F| \times 3$ double matrix represents an $N$-vector field, known as an $N$-RoSy. The single vector is the "first" vector in the face, and the rest of the vectors are deduced by rotations of $\frac{2\cdot\pi}{N}$
3. **Rotation Angles**. A $|E_I|$-sized double vector representing the rotation angle between two directions (without magnitude information) on two neighbouring triangles. The rotation represents the deviation from the Levi-Civita parallel transport [#levy_2008], [#crane_2010]. The type may only encode $N$-direction fields. Note that the <i>effort</i> (sum of all rotations) is then $N$ times rotation angles. Since this is a differential quantity, an extra global rotation needs to be given to uniquely create the full raw field.
4. **Power Field** - A $|F|$-sized *complex* vector, representing an $N$ rosy as a single complex number $y=u^N$, where all the possible roots $u$ comprise an $N$-RoSy. The magnitude is also encoded this way, though it may be neglected in some applications.
5. **PolyVector** - A $|F| \times N$ complex double matrix, encoding the coefficients $X$ of a complex polynomial $f(z)=\sum_{i=0}^{N-1}{X_iz^i}$, which roots $u$ are an $1^N$-vector field. Every row is encoded as $X_{0},\cdots, X_{N-1}$, where $X_0$ is the free coefficient. In case where the field is an $N$-RoSy, all coefficients but $X_0$ are zero.

libdirectional provides a number of conversion functions to switch between different forms of representation. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the representation names in the above list. e.g. `rotation_to_representative()` and `polyvector_to_raw()`. Some possible combinations are given by composing two functions in sequence.

However, note that every conversion is possible; for instance, it is not possible to convert between PolyVectors and rotation angles, as they do not represent the same generality of directional fields (with current state-of-the-art...). For $N$-RoSy fields, for instance, you will most likely work primarily with the power field, representative, and rotation-angle representation. converting into the more explicit raw representation is done mostly for I\O and visualization.

In the following subchapter, we show some effects of working with different representations and converting between them.

## [201 Principal Matching](#principalmatching)[principalmatching]

One of the fundamental operations in directional-field processing is <i>matching</i>. That is, defining which vectors in face $f_i$ correspond to those in adjacent face $f_j$. In libdirectional we only treat order-preserving matching: if vector $k$ in face $f_i$ is matched to vector $m$ in face $f_j$, then for any $l$, $k+l$ is matched to $m+l$ (modulu $N$). Suppose that the orientation of the dual edge is $f_i \rightarrow f_j$, then the matching is encoded as $m-k$. Some representations, like rotation angles, already encode the matching explicitly, but others do not. Therefore, it needs to be devised from the field.

Given a raw field (in assumed CCW order in every face), it is possible to devise the rotation angles $\delta_{ij}$ by the process of *principal matching* [#diamanti_2014]. Principal matching creates the matching that minimizes the effort of the matching, always putting it within the range of $[-\pi, \pi)$ (and therefore denoted as "principal"). It corresponds to the "smallest angle" matching.

principal matching is done through this function (from in [Example 201](201_PrincipalMatching/main.cpp):

```cpp
directional::principal_matching(V, F,EV, EF, FE, rawField, matching, effort);
directional::effort_to_indices(V,F,EV, EF, effort,N,prinIndices);
```

The second line devises the <i>indices</i> of each vertex from the effort. The index of a vertex is the amount of rotations a directional undergoes along a cycle around the vertex. It must return to itself, and therefore the index is an integer $I$ when a vector $m$ in the face ended up in vector $m+I$ (modulu $N$). Note that this can also include multiple full rotation, and therefore the index is unbounded. The <i>fractional</i> part of the index is encoded by the matching; however, matching alone cannot encode indices (for instance, a single vector field has trivial matching anywhere, but can have singularities).

![([Example 201](201_PrincipalMatching/main.cpp)) Field is shown with singularities, and a single face with principal matching to its neighbors (in multiple colors). Manually follow the (natural) principal matching to see that the singularities are indeed correct.](images/201_PrincipalMatching.png)


## [202 Sampling](#sampling)[sampling]

## [203 Combing](#combing)[combing]

# Chapter 3: Cartesian Representations [chapter3:cartesian]

## [Cartesian Fields](#cartesian)[cartesian]
## [301 Globally Optimal Fields](#globallyoptimal)[globallyoptimal]

The algorithm libdirectionl uses to compute such fields is "Globally Optimal"<sup>[6](#fn6)</sup>. The field is defined as a single complex number $y$ (relative to local basis), that represents an $N$-RoSy by its roots. By prescribing constraints $y_B$ on a set of faces $B$, the algorithm interpolates to the rest of the faces $y_I$ by minimizing the Dirichlet energy $y_I=\text{argmin}\ | \nabla y |^2$. Note that $y$ also encodes magnitude in general; however, it is possible to normalize the field  after the computation, as long as it is not identically zero.


### Defining Constraints
Constraints are defined as a pair of matrices of equal hight, refered to as `soft_ids` and `soft_values`. `Soft_ids` is a 1 wide integer matrix containing the ids of the faces (index in the F matrix) on which the constraints are placed. Meanwhile `soft_values` contains the x, y and z values of the matching vector, representing it in the same way as the representative vectors. Constraints do not need to be of unit length, but their size does matter for how they affect the field.

### Precomputing Solver
It is possible to speed up computations by precomputing the solver used to compute the complex field. This solver can then be reused to compute changes in the field as long as the mesh and **the ids of the constrained faces** remain the same.

### Examples
The *complex_field* example contains a small program which allows setting the soft constraints dynamically and see how it affects the field.

The below code creates a field of degree 3 and sets the first face so that its first vector aligns with the first edge. The other vectors are equally spaced to create a 3-rosy. V and F are the Vertices and Faces of the mesh. To see an example that alligns all vectors on the first face with an edge see the polyvector field example.

Without precalculations:
```cpp
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1));

// Matrix containing the field
Eigen::MatrixXcd complex;

//Calculate the field
directional::complex_field(V, F, soft_ids, soft_values, N, complex);
```

With precalculations:
```cpp
//Degree of the field (number of vectors within each directional)
int N = 3;

Eigen::MatrixXi TT;
igl::triangle_triangle_adjacency(F, TT);
Eigen::MatrixXd B1, B2, x;
igl::local_basis(V, F, B1, B2, x);
    
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1));
    
// Prepare the solver, must be recalculated whenever soft_ids changes (optional)
Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> solver;
Eigen::SparseMatrix<std::complex<double>> energy;
Complex_field_prepare_solver(V, F, TT, B1, B2, soft_id, N, solver, energy);

// Calculate the field
complex_field(B1, B2, soft_id, soft_value, solver, energy, N, complex);
```
## [302 PolyVectors](#polyvectors)[polyvectors]

## Polyvector Field
The Polyvector field is a generalisation of the standard complex field method, which allows defining each vector in a directional individually in both direction and length. Besides that it works largely the same way as the complex field.<sup>[5](#fn5)</sup> 

### Defining Constraints
Just like the Complex Field, the Polyvector Field takes a `soft_ids` matrix defining the face indices and a matching `soft_values` matrix defining the directionals on the faces, however the polyvector `soft_values` matrix is 3\*N wide, containing the X, Y, and Z values for each individual vector.

### Precomputing Solvers
It is possible to precompute the solvers for the Polyvector Field. To precompute the solvers one should pass an empty vector of SimplicialLDLT pointers (`std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*>`) into the `poly_vector_prepare_solvers()` function before passing them along with the energy matrix to the `poly_vector()` method.

The Solvers can be reused as long as the `soft_ids` remain the same, and must be properly `deleted` afterwards.

### Examples
The *poly_vectors* example shows the polyvector field in action, allowing the user to set constraints for each vector on each face individually. 

The below code creates a field of degree 3 and sets the first face so that its vectors each align with one edge of the triangle.  V and F are the Vertices and Faces of the mesh.

Without precalculations:
```cpp
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1)), V.row(F(0,1)) - V.row(F(0,2)), V.row(F(0,2)) - V.row(F(0,0));

// Matrix containing the field
Eigen::MatrixXcd complex;

//Calculate the field
directional::poly_vector(V, F, soft_ids, soft_values, N, complex);
```

With precalculations:
```cpp
//Degree of the field (number of vectors within each directional)
int N = 3;

Eigen::MatrixXi TT;
igl::triangle_triangle_adjacency(F, TT);
Eigen::MatrixXd B1, B2, x;
igl::local_basis(V, F, B1, B2, x);
    
// Set constraints
Eigen::VectorXi soft_ids(1);
Eigen::MatrixXd soft_values(1, N*3);
// Set to all faces that should be constrained
soft_ids(0) = 0;
// Set each matching row to the N vectors on that face
soft_values << V.row(F(0,0)) - V.row(F(0,1)), V.row(F(0,1)) - V.row(F(0,2)), V.row(F(0,2)) - V.row(F(0,0));
    
// Prepare the solvers, must be recalculated whenever soft_ids changes (optional)
std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*> solvers;
std::vector<Eigen::SparseMatrix<std::complex<double>>> energy;
poly_vector_prepare_solvers(V, F, TT, B1, B2, soft_ids, N, solvers, energy);

// Matrix containing the field
Eigen::MatrixXcd poly;

// Calculate the field
poly_vector(B1, B2, soft_ids, soft_values, solvers, energy, N, poly);

...

// Make sure to properly dispose of all solvers
for (std::vector< Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>* >::iterator it = solvers.begin(); it != solvers.end(); ++it)
{
    delete (*it);
}
```

## [303 Integrable PolyVectors](#integrablePVs)[integrablePVs]

<span style="color:red">
This tutorial is a direct migration of the libigl version, as the functionality would reside in libdirectional. It is fully functional, but the syntax is not libdirectional-compatible as of yet.
</span>

Vector-field guided surface parameterization is based on the idea of designing the gradients
of the parameterization functions (which are tangent vector fields on the surface) instead of the functions themselves. Thus, vector-set fields (N-Rosy, frame fields, and polyvector fields) that are to be used for parameterization (and subsequent remeshing) need to be integrable: it must be possible to break them down into individual vector fields that are gradients of scalar functions. Fields obtained by most smoothness-based design methods (eg. [#levy_2008][], [#knoppel_2013][], [#diamanti_2014][], [#bommes_2009][], [#panozzo_2014][]) do not have this property. In [#diamanti_2015][], a method for creating integrable polyvector fields was introduced. This method takes as input a given field and improves its integrability by removing the vector field curl, thus turning it into a gradient of a function ([Example 303](303_IntegrablePVs/main.cpp)).

![Integration error is removed from a frame field to produce a field aligned parameterization free of triangle flips.](images/303_IntegrablePVs.png)

This method retains much of the core principles of the polyvector framework - it expresses the condition for zero discrete curl condition (which typically requires integers for the vector matchings) into a condition involving continuous variables only. This is done using coefficients of appropriately defined polynomials. The parameterizations generated by the resulting fields are exactly aligned to the field directions and contain no inverted triangles.

## [304 Conjugate Fields](#conjugatefields)[conjugatefields]

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

## [401 Trivial Connection](#trivialconnection)[trivialconnection]

The notation of encoding rotation angles appeared in several formats in the literature (see STAR for details). The formulation and notation we use here is the of Trivial Connections<sup>[4](#fn4)</sup>. Trivial connection solve for a single rotation angle $\delta_{ij}$ per (dual) edge $e_{ij}$ between two faces $f_i,f_j$ that encodes the deviation from parallel transport. The algorithm defines a spanning set of basis cycles (see next section), around which the sum of $\delta_{ij}$ is prescribed. The sum defines the *curvature* of the cycle. The algorithm solves for the smoothest field, or the least-squares 2-norm $\delta$:

$$
\delta = \text{argmin}\ |\delta_{ij}|^2\ s.t.\ H\delta = -K_0 + K.
$$

$H$ is the matrix that defines the basis-cycles sum, $-K_0$ is the original curvature of the basis cycle, and $K$ is the prescribed curvature. $K$ defines singularities: for a regular cycles, we set $K=0$, and for a singular cycle with singularity index $\frac{1}{N}$, we set $K=\frac{2\pi}{N}$. the sum of $K$ has to conform to the Poincar&eacute; index theorem, except handle cycles which can have unbounded index. See paper for exact details. As a consequence, the singularities of the field, and the indices of other basis cycles, are the input. The output is not yet an $N$-RoSy field: there is a global degree of freedom in setting a single direction in a single face.

The algorithm is defined with "cycle holonomy" that is the curvature "modulu $\pi$". They are the same the same for small cycles, but not so for big cycles. However, curvature (specifically, angle defect) is the correct measure needed for the consistency of trivial connection, and we use it.

### Basis Cycles

The basis cycles form the cycles around which curvatures (and singularities) are described on the mesh. The sum on basis cycles is described in a sparse matrix $H$ of size $|cycles|\times |E|$. Each row in the matrix describes the sum over one cycle, and contains a 1 or -1 depending on the orientation of the dual edge participating in the cycle. There are three types of cycles, according to their order in the rows of $H$: 

1. $1$-ring dual cycles around each vertex, on which vertex-based singularities can be encoded (the relevant part of $H$ is basically $d_0^T$ in discrete exterior calculus).
2. Cycles around mesh boundary loops. 
3. Cycles around topological generator (independent handles).

The method `dual_cycles()` computes the proper basis cycles and matrix $H$.

### Singularities

The singularity indices are prescribed as an `Eigen::VectorXi` object, containing the singularity index corresponding to each basis cycle. A value of $k \in \mathbb{Z}$ represents an $\frac{2\pi k}{N}$ rotation around the respective cycle. In order to create a smooth field it is required that the indices of all singularities add up to $\sum{k}=N *\chi$, where $\chi$ is the Euler characteristic of the mesh. Otherwise, a result will still be computed by least squares, but it will be unpredictable.

###Trivial connection demo

Below is a contraction of code from `examples\trivial_connection` that computes and draws a field from prescribed singularities:

```cpp
igl::edge_topology(meshV, meshF, EV, FE, EF);
directional::dual_cycles(meshF, EV, EF, cycles);
igl::boundary_loop(meshF, boundaryLoops);

directional::trivial_connection(meshV, meshF, cycles, indices, N, rotationField, e);
directional::rotation_to_representative(meshV, meshF, EV, EF, rotationField, N, globalRotation, representative);

int sum = round(indices.head(indices.size() - generators).sum());
if (euler*N != sum){
    std::cout << "Warning: All non-generator singularities should add up to N * the Euler characteristic."<<std::endl;
    std::cout << "Total indices: " << sum << std::endl;
    std::cout << "Expected: " << euler *N << std::endl;
}

// Turn the field into a drawable mesh
directional::drawable_field(meshV, meshF, representative, Eigen::RowVector3d(0, 0, 1), N, directional::field_draw_flags::NONE, fieldV, fieldF, fieldC);

```


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
[#crane_2010]: Keenan Crane, Mathieu Desbrun, Peter Schr&ouml;der, [Trivial Connections on Discrete Surfaces](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/), 2010.<br>
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

