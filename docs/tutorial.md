
<h1>Directional tutorial notes</h1>


## Introduction

Directional is a C++ geometry processing library written as an extension library to [libigl](http://libigl.github.io/libigl/), with a speciality in directional-field processing. The functionality is based on the definitions and taxonomy surveyed theoretically in [^vaxman_2016], and through it by much of the relevant papers in the literature. It contains tools to edit, analyze, and visualize directional fields of various degrees and symmetries.

The underlying structure extends the general philosophy of [libigl](http://libigl.github.io/libigl/): the library is header only, where each header contains a set (often only one) of functions closely related (for instance, the precomputation and computation of some directional quantity over a mesh). The data structures are, for the most part, simple matrices in [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), and the library avoids complicated and nested structures, instead directly working with calling to standalone functions. The visualization is done using the libigl viewer with some extended options that allow the rendering of directional fields.

The header files contain documentation of the parameters to each function and their required composition; in this tutorial we will mostly tie the functionality to the theoretical concepts of directional fields and their capabilities.

### Installing the tutorial examples

This tutorial comprises an exhaustive set of examples that demonstrate the capabilities of Directional, where every subchapter relates to a single example. The tutorial code can be installed by going into the `Tutorial` folder, and following the following instructions:


```cpp
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

This will build all tutorial examples in the `build` folder. To build in windows, use the `cmake-gui .` options instead of the last two, and creates the project using Visual Studio, with the proper tutorial chapter as the "startup project".

To access a single example, say ``202_Sampling``, go to the ``build`` subfolder, and the executable will be there. Command-line arguments are never required; the data is read from the ``shared`` folder directly for each example.

Most examples contain a component of user interaction; the instructions of what to do will be given in the command-line output upon execution.

### Discretization

There are several ways to encode directional fields on discrete triangle meshes, although not all of them compatible with the different information that the fields conveys. Mainstream discretizations roughly divide into face-based, edge-based, and vertex-basex discretization (see [^degoes_2016] for an in-depth analysis). The only discretization currently supported by Directional is that of face-based fields, where the discrete tangent plane is considered as the supporting plane to each (naturally flat) face. We also use a local basis (mostly provide by `igl::local_basis()`) to parameterize each tangent plane. Notions like connection, parallel transport, and consequently smoothness, are measured through *dual edges* between adjacent faces.

### Representation

The representation of a directional field is the way it is encoded in each discrete tangent plane, and cross adjacent tangent planes. Directional uses several different representations to describe directional fields. We denote the number of faces in the mesh as $|F|$, the set of inner edges (adjacent to two triangles) as $E_I$, and the degree of the field (must be fixed for the entire field) as $N$. The supported  representations are as follows, where the taxonomy is based on that of the directional field course [^vaxman_2016]:

1. **Raw** - A $|F|\times 3N$ double matrix, representing an $1^N$-vector field (a directional with $N$ independent vectors in each face) in the form $X_1, Y_1, Z_1, X_2, Y_2, Z_2, \cdots X_N, Y_N, Z_N$ per face. Vectors are assumed to be ordered in counterclockwise order in most Directional functions that process raw fields.
2. **Representative**. A $|F| \times 3$ double matrix that represents a rotationally symmetric $N$-vector field, known as an $N$-RoSy. The single vector is an arbitrary "first" vector in the face, and the rest of the vectors are deduced by rotations of $\frac{2\cdot\pi}{N}$. That means the all functions should also accept $N$ as input to decode the full field.
3. **Rotation Angles**. A $|E_I|$-sized double vector representing the rotation angle between two directions (without magnitude information) on two adjacent triangles. The rotation represents the deviation from the Levi-Civita parallel transport [^ray_2008], [^crane_2010]. This representation may only encode $N$-direction fields. Note that the *effort* (sum of all rotations) is then $N$ times rotation angles. Since this is a differential quantity, an extra global rotation needs to be given to uniquely decode the full face-based field.
4. **Power Field** - An $|F|$-sized *complex* vector, representing an $N$-RoSy as a single complex number $y=u^N$, where the $N$-RoSy is the set of roots $u$. The magnitude is also encoded this way, though it may be neglected in some applications. The representation depends on a local $2D$ basis, such as one that could be obtained from ```igl::local_basis()```.
5. **PolyVector** - A $|F| \times N$ complex double matrix, encoding the coefficients $a$ of a monic complex polynomial $f(z)=z^N+\sum_{i=0}^{N-1}{a_iz^i}$, which roots $u$ are an $1^N$-vector field. Every row is encoded as $a_{0},\cdots, a_{N-1}$, where $a_0$ is the free coefficient. In case where the field is an $N$-RoSy, all coefficients but $a_0$ are zero.

Directional provides a number of conversion functions to switch between different representations. Each of the functions is of the form ```rep1_to_rep2```, where ```rep1``` and ```rep2``` are the representation names in the above list. e.g., ```rotation_to_representative()``` and  ```polyvector_to_raw()```. Some possible combinations are given by composing two functions in sequence.

However, note that every conversion is possible; for instance, it is not possible to convert between PolyVectors and rotation angles, as they do not possess the same generality of directional fields (with current state-of-the-art...). For $N$-RoSy fields, for instance, you will most likely work primarily with the power field, representative, or rotation-angle representation. converting into the more explicit raw representation is often needed for I/O and visualization.


## Chapter 1: I/O and Visualization

### Visualization paradigm

Directional uses the libigl viewer (although a viewer is not necessary for the core functionality), and augments it with auxuliary functionality and colors schemes that pertain to directional fields specifically. The drawing paradigm in Directional is that the visualization functions create actual meshes (vertices, faces, and per-face colors) for the different components of the field: box-like meshes for per-face glyphs representing vectors, sphere meshes for representing singularities, piecewise-cylinder curves to represent streamlines, and more. These visualization meshes can then be stored in libigl viewer as multiple meshes and visualized as needed.

The color schemes and directional visualization mesh creation functionality are mostly given in [visualization_schemes.h]({{ repo_url }}/include/directional/visualization_schemes.h), and can be called to create a homogeneous look for visualizing fields throughout applications. Most notable are ```indexed_glyph_colors()``` that gives back glyph meshes colored by index per face, and ```default_singularity_colors()``` that produces singularity meshes according to index. Other visualization functionality is detailed below in the context of the different examples.

### 101 Glyph Rendering

The most basic operation on directional fields is reading them from a file and drawing them in the most explicit way. In [Example 101]({{ repo_url }}/tutorial/101_GlyphRendering/main.cpp), a field is read from a file as follows:

```cpp
directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singVertices, singIndices);
```

The field is read in *raw* format (see File Formats), which is detailed in the [Introduction](#introduction). The field is *face-based*, and the singularities are consequently *vertex-based*, where ```singVertices``` are the singular vertices, and ```singIndices``` are the corresponding integer indices, so that the actual fractional indices are $\frac{singIndices}{N}$.

The visualization meshes for the glyphs and the singularities are obtained as follows

```cpp
directional::glyph_lines_raw(VMesh, FMesh, rawField, directional::default_glyph_color(), VField, FField, CField);
directional::singularity_spheres(VMesh, FMesh, N, singVertices, singIndices, VSings, FSings, CSings);
```

These two operations do not produce any active drawing; they create meshes that extend the original geometry, and then get passed to libigl viewer.

``directional::glyph_lines_raw()`` creates box-like meshes on the faces that constitute the *glyph drawing*: simply drawing the vectors upon the faces in their $\left(x,y,z\right)$ coordinates, starting from the face barycenter. There are several ways to give colors to these vectors, which can be individual or global; check the documentation to the function in the header. In this case, we give the ```default_glyph_color()``` from the Directional visualization schemes. Vectors are drawn in their given magnitudes, up to a global scale. By default, this scale is set to be related to the average edge length---it can be manually set by the extended version of the function.

``directional::singularity_spheres()`` creates a mesh of small spheres on vertices, where the size of the sphere is devised automatically (but can be configured using the extended version of this function). The spheres are only created where the index is different than $0$.


![[Example 101]({{ repo_url }}/tutorial/101_GlyphRendering/main.cpp) Glyph Rendering on a mesh, with singularities visible.](images/101_GlyphRendering.png)


### 102 Picking and editing

This tutorial demonstrates the editing paradigm in Directional, using libigl picking to make directional field editing possible. A face and a vector within the face are chosen, and clicking on a new direction for the vector changes it. Note the different setting of colors for glyphs: for the selected face and for the selected vector in the face particularly, via the following code in [Example 102]({{ repo_url }}/101_PickingEditing/main.cpp).

```cpp
Eigen::MatrixXd glyphColors=directional::default_glyph_color().replicate(FMesh.rows(),N);
glyphColors.row(currF)=directional::selected_face_glyph_color().replicate(1,N);
glyphColors.block(currF,3*currVec,1,3)=directional::selected_vector_glyph_color();

directional::glyph_lines_raw(VMesh, FMesh, rawField, glyphColors, VField, FField, CField);
```

The selected face coloring is done as follows:

```cpp
CMesh=directional::default_mesh_color().replicate(FMesh.rows(),1);
CMesh.row(currF)=directional::selected_face_color();
```


![([Example 102]({{ repo_url }}/102_PickingEditing/main.cpp)) Editing several vectors on a single face.](images/102_PickingEditing.png)

### 103 Streamline Tracing

Vector fields on surfaces are commonly visualized by tracing [streamlines] (https://en.wikipedia.org/wiki/Streamlines,_streaklines,_and_pathlines). Directional supports the seeding and tracing of streamlines, for all type of directionals. The seeds for the streamlines are initialized using `streamlines_init`, and the lines are traced using `streamlines_next`. Each call to `streamlines_next` extends each line by one triangle, allowing interactive rendering of the traced lines, as demonstrated in [Example 103]({{ repo_url }}/103_StreamlineTracing/main.cpp). The streamline are compiled into meshes with the following:

```cpp
  directional::line_cylinders(sl_state.start_point, sl_state.end_point, 0.0005, color.replicate(sl_state.start_point.rows(),1), 4, VFieldNew, FFieldNew, CFieldNew);
```
creating meshes from many cylinders to have the appearance of continuous curves.

![([Example 103]({{ repo_url }}/103_StreamlineTracing/main.cpp)) Interactive streamlines tracing.](images/103_StreamlineTracing.gif)

### 104 Dynamic Visualization

This tutorial simulates the behavior of particles within a n-degree field and is based on the simulation of the earth elements as seen on https://earth.nullschool.net/. This was implemented in the dynamic_visualization.h file. There are two main functions. The first is the initialize function, which creates the initial noodles and colors the mesh. By default the mesh is colored based on the average effort (average of the effort per edge) per face it takes for a noodle to move from one face to the next. This can also be done by giving it an custom user defined function as an extra argument. The initialization without the user function for the coloring looks like this:

```cpp
  directional::initialize_noodles(n_data, VMesh, CMesh, FMesh, streamLengths, degree, MaxLifespan, 1.0);
```


The second function is the update_noodles function. This function uses the same method as used in [103 Streamline Tracing] to make the noodles move. Instead of showing the entire path, it only shows part of it. It also shades the tail of the noodle to indicate the travel direction. The update function is called like so:

```cpp
  directional::update_noodles(n_data, VMesh, FMesh);
```

As a result we get a mesh over which the noodles are traveling which are representing the flow in the vector field.
![([Example 104]({{ repo_url }}/104_DynamicVisualization/main.cpp)) Dynamic Visualization](images/104_DynamicVisualization.gif)


## Chapter 2: Discretization and Representation


In the following sections, we show some effects of working with different representations and converting between them.

### 201 Principal Matching

One of the fundamental operations in directional-field processing is *matching*. That is, defining which vectors in face $f_i$ correspond to those in adjacent face $f_j$. In Directional, we only work with order-preserving matchings: if vector $k$ in face $f_i$ is matched to vector $m$ in face $f_j$, then for any $l$, $k+l$ is matched to $m+l$ (modulu $N$) in the respective faces. Suppose that the orientation of the dual edge is $f_i \rightarrow f_j$, then the matching is encoded as $m-k$. Some representations, like rotation angles, already encode the matching explicitly, but others do not. Therefore, it needs to be devised from the field.

Given a raw field (in assumed CCW order in every face), it is possible to devise the rotation angles $\delta_{ij}$ by the process of *principal matching* [^diamanti_2014]. Principal matching creates the matching that minimizes the effort of the matching, always putting it within the range of $[-\pi, \pi)$ (and therefore denoted as "principal"). It corresponds to the "smallest angle" matching for $N$-RoSy fields.

principal matching is done through the function ```principal_matching()```  (from in [Example 201]({{ repo_url }}/201_PrincipalMatching/main.cpp)) as follow:

```cpp
directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawField, matching, effort);
directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching,N,singVertices,singIndices);
```
`directional::effort_to_indices()` computes the <i>index</i> of each vertex from the effort around it. The index of a vertex is the amount of rotations a directional object undergoes along a cycle around the vertex. a directional must return to itself after a cycle, and therefore the index is an integer $I$ when a vector $m$ in the face ended up in vector $m+I$. Note that this can also include multiple full rotations(so: this is *not* taken modulu $N$), where the index is unbounded. The *fractional* part of the index is encoded by the matching; however, matching alone cannot encode *integral* indices (for instance, a single vector field has trivial (Zero) matching anywhere, but can have singularities). ```singVertices``` and ```singIndices``` only enumerate the singular vertices.

![([Example 201]({{ repo_url }}/201_PrincipalMatching/main.cpp)) A Field is shown with singularities, and a single face is shown with the principal matching to its neighbors (in multiple colors).](images/201_PrincipalMatching.png)


### 202 Sampling

This is an educatory example that demonstrates the loss of information when moving between a polar (in this case, rotation angle) representation, to a Cartesian representation, where the matching between vectors in adjacent faces is done with principal matching. In that case, low valence cycles and undersampling cause aliasing in the perceived field. There are three modes seen in the example:

1. In the polar mode, the user can prescrube the index of a singularity directly. With this, the rotation angles between adjacent faces become arbitrarily large, and appear as noise in the low valence cycles.

2. In the principal matching mode, the rotations are recoconstructed from the field, without prior knowledge of the polar-prescribed rotations from the previous mode. The noise of the field then gives rise to a "singularity party".

3. In the Cartesian mode, the field is interpolated on the free faces (white) from the constrained faces (red), keeping the red band fixed from the polar mode. We see a field that in smooth in the Cartesian sense, with more uniformly dispersed singularities.

![([Example 202]({{ repo_url }}/202_Sampling/main.cpp)) Animating: the polar mode, the principal matching mode, and the Cartesian mode.](images/202_Sampling.gif)

### 203 Combing

Given a matching (in this case, principal matching), it is possible to "comb" the field. That is, re-index each face (keeping the CCW order), so that the vector indexing aligns perfectly with the matching to the neighbors---then, the new matching on the dual edges becomes trivially zero. This operation is important in order to prepare a directional field for integration, for instance. In the presence of singularities, the field can only be combed up to a set of connected paths that connect between singularities, also known as *seams*. Note that such paths do not necessarily the mesh to a simply-connected patch, but only connects subgroups with indices adding up to an integer; as a trivial example, a 1-vector field is trivially combed to begin with, even in the presence of its (integral) singularities, and the set of seams is zero. The combing is done through the function ```directional::combing()``` as follows, taken from [Example 203]({{ repo_url }}/203_Combing/main.cpp)

```cpp
directional::combing(VMesh,FMesh, EV, EF, FE, rawField, matching, combedField);
```

where ```combedField``` is the re-indexed ```rawField```, done according to the matching. We can re-do the matching on the combed field to retrieve the seams:

```cpp
directional::principal_matching(VMesh, FMesh,EV, EF, FE, combedField, combedMatching, combedEffort);
```

![([Example 203]({{ repo_url }}/203_Combing/main.cpp)) Vector indices inside each face are colored, combed (with seams) and uncombed).](images/203_Combing.gif)

## Chapter 3: Cartesian Methods

### Cartesian Fields

The Cartesian representation is a meta-category for representation of vectors in explicit coordinates, either $\left(x,y\right)$ in some local $2D$ basis on a tangent plane, or $\left(x,y,z\right)$ in the ambient coordinates of space. The raw, representative (of an $N$-RoSy), power field, and PolyVector representations are all such examples. Cartesian fields often do not automatically contain information about the matching, or rotation, of a field between one face and the next, and it needs to be computed using principal matching. This chapter focuses on computing fields with this representation.

### 301 Power Fields

This representation is offered in [^knoppel_2013], but they did not give it a specific name (the method in general is called "globally optimal"). We use the name "power fields" given in [^azencot_2017].

A power field representation uses a complex basis in each tangent plane (face in our implementation), and represents an $N$-RoSy using a *power vector*---a single complex number $y$ per face so that its root set $y=u^N$ comprises the vectors of the $N$-RoSy.

By prescribing constraints $y_B$ on a set of faces $B$, the algorithm interpolates the field to the rest of the faces $y_I$ by minimizing the face-based Dirichlet energy:

$$y_I=\text{argmin}\sum_{e=(f,g)\in F \times F}{\left|y_fe_f^N - y_ge_g^N\right|^2},$$

where $e_f$ is the representation of the vector of edge $e$ in the basis of $f$, and $e_g$ is for $g$ respectively. The field is computed through the function `directional::power_field()`.

For fixed set $B$ and changing $y_B$, It is possible to speed up computations by precomputing the solver (sparse Cholsky for the positive-definite matrix) used to compute the power field. This is done by using the function `directional::power_field_precompute()`, coupled with the appropriate version of `directional::power_field()`. Note that field can be converted to representative and raw forms using the appropriate `power_to_X` functions.

If the set $B$ is empty, then the computed field is the first Eigenvector of the Dirichlet energy.

![([Example 301]({{ repo_url }}/301_PowerFields/main.cpp)) Setting up a small subset of constraints (red faces), and interpolating (and normalizing in magnitude) the power field to the rest of the mesh. Note the singularities that are discovered through principal matching.](images/301_PowerFields.png)


### 302 PolyVectors

A Polyvector field [^diamanti_2014] is a generalization of power fields that allows to represent independent vectors in each tangent plane. The representation is as the coefficient set $a_{0 \cdots N-1}$ of a monic complex polynomial in the local compex basis:

$$P(z) = a_0 + a_1z + \ldots + a_{N-1} z^{N-1} + z^N,$$

where the roots $P(z)=0$ are the vectors of the face-based directional object, represented as complex numbers in the local basis. The Dirichlet energy is as for power fields, except with a term for each $a_i$, with the appropriate power $i$. Note that an $N$-RoSy is represented as a polynomial where all $a$ are zero except $a_0$. Principal matching, combing, and effort are well-defined on PolyVectors as well.

[Example 302]({{ repo_url }}/302_PolyVectors/main.cpp) allows a user to set individual vectors within each face, and see the interpolated result. The responsible function is `directional::polyvector_field()`. In this case as well, the solver can be prefactored in advance using `directional::polyvector_precompute()`.

![([Example 302](302_PolyVectors/main.cpp)) Vectors are constrained individually in the constrained faces (red), and interpolated to the rest of the faces](images/302_PolyVectors.png)

### 303 PolyCurl Reduction

Vector-field guided surface parameterization is based on the idea of designing the *candidate* gradients
of the parameterization functions (which are tangent vector fields on the surface) instead of the functions themselves. Thus, vector-set fields ($N$-Rosy, frame fields, and polyvector fields) that are to be used for parameterization (and subsequent remeshing) need to be *integrable*: it must be possible to locally break them down into individual vector fields that are gradients of scalar functions. Fields obtained by "as-smooth-as-possible" design methods (eg. [^ray_2008], [^knoppel_2013], [^diamanti_2014], [^bommes_2009], [^panozzo_2014]) do not have this property in general. In [^diamanti_2015], a method for creating integrable polyvector fields was introduced by the process of *curl-reduction*. This method takes as input a given field and improves its integrability by iteratively reducing the *PolyCurl* of the field; that is, the coefficients of a dual-edge-based polynomial, whose roots are the curl of the matched vectors of both respective adjacent faces. By working with PolyCurl instead of matching, the optimization can be done on the POlyVector itself, allowing for singularities to naturally move around. However, the optimization is nonlinear---it reduces the PolyCurl iteratively, while preserving the CCW order of the vectors for a bijective parameterization, and keeping it as smooth and orthogonal as possible for quality results.

A field that has zero PolyCurl everywhere is locally (away from singularities) integrable into $N$ different scalar functions; globally, it is integrable into a *rotationally-seamless* multi-branched function, which we further demonstrate in [Chapter 5](#chapter-5-seamless-parameterization). In [Example 303]({{ repo_url }}/303_PolyCurlReduction/main.cpp) we demonstrate the optimization for curl.

![([Example 303](303_PolyCurlReduction/main.cpp)) PolyCurl is iterative reduced from an initial PolyVector field, Top: fields for iteration 0 (original), 10, and 50. Bottom: PolyCurl plots. The color is the norm of the vector of the roots of the PolyCurl. Its maximal value (infinity norm) appears below.](images/303_PolyCurlReduction.png)

### 304 Conjugate Fields


Two tangent vectors $u$ and $v$ are conjugate if

$$k_1 (u^T d_1)(v^T d_1) + k_2(u^T d_2)(v^T d_2) = 0. $$,
where $k_1$ and $k_2$ are the principal curvatures and $d_1$ and $d_2$ are the respective principal directions.

Conjugate vector fields are very important in computational geometry: their integral lines form, informally speaking, an infinitesimal planar quad mesh. As such, the finite quad mesh that result out of integrating these fields is considered as a good candidate for consequent planarity parameterization [^liu_2011].

Finding a conjugate vector field that satisfies given directional constraints is a standard problem in architectural geometry, which can be tackled by
deforming a $2^2$ PolyVector field to the closest conjugate field. Such an algorithm was presented in [^diamanti_2014], which alternates beweetn a global smoothness and orthogonality step, and a local step that projects the field on every face to the closest conjugate field ([Example 304]({{ repo_url }}/304_ConjugateField/main.cpp)).


![([Example 304]({{ repo_url }}/304_ConjugateFields/main.cpp)) A smooth $2^2$-PolyVector field (left) is deformed to become a conjugate field (right). Top: fields Bottom: conjugacy plots.](images/304_ConjugateFields.png)


## Chapter 4: Polar Methods

### Polar Fields

Polar fields are represented using angles. These angles may encode the rotation from some known basis on a tangent plane (and so it is a "logarithmic" representation, when compared to Cartesian methods), or an angle difference between two neighboring tangent planes (in the sense of deviation from parallel transport). The former usually requires integer variables for directional field design. The latter does not, but state-of-the-art methods require the prescription of indices around independent dual cycles in the mesh. Currently, Directional supports the latter.

### 401 Index Prescription

The notion of encoding rotation angles on dual edges, as means to encode deviation from parallel transport between adjacent tangent planes, appeared in several formats in the literature [^ray_2008], [^crane_2010]. The formulation and notation we use in Directional is that of Trivial Connections [^crane_2010]. Trivial connection solves for a single rotation angle $\delta_{ij}$ per (dual) edge $e_{ij}$ between two faces $f_i$ and $f_j$, encoding the deviation from parallel transport between them. The algorithm first computes a spanning set of *basis cycles*, around all of which the sum of $\delta_{ij}$ has to be prescribed. The summation is defined as a matrix $H$. Every such cycle (row in the matrix) has a curvature, defined as a discrete angle defect, and the prescribed index defines an alternative curvature. The algorithm solves for the smoothest field, in the 2-norm least squares sense, as follows:

$$
\delta = \text{argmin}\ |\delta_{ij}|^2\ s.t.\ H\delta = -K_0 + K.
$$

$H$ is the matrix that defines the basis-cycles sum, $K_0$ is a vector of the original curvatures of every basis cycle, and $K$ is the prescribed curvature. $K$ defines the prescribed singularity indices: for regular cycles, we prescribe $K=0$, and for a singular cycle with prescribed singularity index $\frac{i}{N}$, we set $K=\frac{2\pi i}{N}$. the sum of $K$ has to conform to the Poincar&eacute; index theorem, except generator (handle) cycles that admit unbounded indices. See [^crane_2010] for exact details. If the input obeys the sum, the result obeys the prescribed indices around the cycles everywhere. The representation is *differential*, and there is then a global degree of freedom in setting a single direction in a single arbitrary face.

Note that the correct definition for "cycle curvature" corresponds to the so-called "cycle holonomy", only up to integer multiples of $2\pi$. However, in the discrete setting, the curvature should theoretically be computed as the exact discrete angle defect, in which for inner vertices we use $2\pi-\sum{\alpha}$, and for boundary vertices we use $\pi - \sum{\alpha}$ ($\alpha$ are the angles at the corners of a vertex). For a cycle aggregating many vertices, such as a boundary-loop cycle, we add up all the defects. That is required for exact discrete Poincar&eacute; index consistency. Note that the boundary indices define how many rotations of the vector field the boundary loop "sees". As an example, a constant field on a simple disc in the plane has all indices $0$ inside, but the boundary index is in fact $1$---This obeys the total index sum $\chi = 2-2g-b = 2-0-1=1$ ($g$ stands for genus and $b$ for number of boundary loops)

#### Basis Cycles

The basis cycles form the cycles around which curvatures (and singularities) have to be prescribed on the mesh. The sum on basis cycles is described in a sparse matrix $H$ of size $|cycles|\times |E_I|$, where $E_I$ is the number of inner edges in the mesh. Each row in the matrix describes the sum over one cycle, and contains $1$ or $-1$ values depending on the (arbitrary) orientation of the dual edge participating in the cycle to the respective face. There are three types of cycles, so ordered in the rows of $H$:

1. $1$-ring dual cycles around each inner vertex, on which vertex-based singularities can be encoded (the relevant part of $H$ is basically $\left(d_0\right)^T$ in discrete exterior calculus, restricted to inner edges).
2. Cycles around mesh boundary loops.
3. Cycles around the $2g$ topological generators (independent handles).

The method `directional::dual_cycles()` computes the proper basis cycles and matrix $H$. To be able to intuitively prescribe singularities to inner vertices, the method also returns a conversion vector ``vertex2cycle``, and the list of indices of inner edges from the list of edges.

The singularity indices that are prescribed contain the singularity index corresponding to each basis cycle. A value of $k \in \mathbb{Z}$ represents an $\frac{2\pi k}{N}$ rotation around the respective cycle. If the prescribed indices do not conform to the Poincar&eacute; index theorem, a result will still be computed by least squares, but it will be unpredictable. The algorithm is performed through the function ``directional::index_prescription()``, which can also accept a solver for precomputation, for the purpose of prefactoring $H$ only once.


![([Example 401]({{ repo_url }}/401_IndexPrescription/main.cpp)) Indices are prescribed on several vertex singularities, and on a generator loop, to match the index theorem.](images/401_IndexPrescription.png)

## Chapter 5: Seamless Parameterization

Directional fields are commonly used to create seamless parameterizations [^bommes_2009],[^Kaelberer_2007],[^Myles_2014]. Recall that [combing](#203-combing) trivializes the matching everywhere but a sparse set of seams. We augment these seams so that the mesh is cut into disc topology. Then, we treat a combed $N$-directional $\left\{v_0,\cdots,v_{N-1}\right\}$ as a set of candidate gradients for $N$ vertex-based functions $\left\{F_0,\cdots,F_{N-1}\right\}$ on the cut mesh. On the cut mesh, we then solve the Poisson problem:

$$ F = argmin{\sum_{i=0}^{N-1}{\left|\nabla F_i - v_i\right|^2}} $$.

Consider a seam edge $e_{ij}$ between original vertices $v_i$ and $v_j$, and between adjacent faces $f_k$ and $f_l$. The two vertices are then cut into four corners $v_{i,k},v_{j,k},v_{i,l},v_{j,l}$. Note that some corners might be identical, if the seam edge is at a singularity. across the seam edge, we enforce the (linear) seamless conditions:

$$F_{i,k}= \pi_e \cdot F_{i,l}  + T_e,$$

where $\pi_e:N \times N$ is a permutation matrix attached to the (dual) edge $e$, matching values in the integrated function $F$ as it did for the directional field $v$. and $T_e:N \times 1$ is a *translational jump* (also: period jump), that encodes the discontinuity in $F$ across the seam. For quick intuition, this encodes the integration of the function over a loop around the mesh beginning and ending with the seam edge: in a quad-mesh parameterization, it is the number of quads in such a loop.

Every parameterization that obeys the seamless constraint is seamless; it can be easily shown [^Kaelberer_2007] that the translational jump $T_e$ is in fact uniform across seam curves between singularities. Thus, the amount of such jumps is the number of seam curves in the mesh.

### 501 Seamless Parameterization

In [Example 501]({{ repo_url }}/501_SeamlessParameterization/main.cpp) we demonstrate the computation of such a parameterization. The core functionality is in these lines:

```cpp
directional::ParameterizationData pd;
directional::cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
  ...
directional::setup_parameterization(N, VMeshWhole, FMeshWhole,  EV, EF, FE, combedMatching, singVertices, pd, VMeshCut, FMeshCut);
 double lengthRatio=0.01;
 bool isInteger = false;  //do not do translational seamless.
 std::cout<<"Solving parameterization"<<std::endl;
 directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, cutUV);
```

```directional::cut_mesh_with_singularities()``` encodes the the seams in ```pd.face2cut```. ```directional::setup_parameterization()``` creates the Poisson system and the constraints, and creates the actual cut mesh in ```VMeshCut``` and ```FMeshCut```.  ```directional::parameterize()``` solves the parameterization, and puts the result in ```cutUV```. The functionality is currently limited only to $2^2$ fields ($N=4$) with symmetry assumed. ```lengthRatio``` encodes a global scale for the Poisson problem (scaling the fields uniformly), where the ratio is measured against the bounding box diagonal.

The variable ```isInteger``` refers set to ```false``` which means that the parameterization would only be *rotationally-seamless*; the appearance across the seams would align only in direction. The next example considers the stronger option.


![([Example 501]({{ repo_url }}/501_SeamlessParameterization/main.cpp)) Left: directional field. Right: rotationally-seamless parameterization. Note that the direction of the texture aligns across seams.](images/501_SeamlessParameterization.png)

### 502 Mixed-Integer Parametrization

To use seamless parameterizations for the purpose of quad meshing, we require the texture to be *fully-seamless* across seams. With this, the seams are virtually invisible on the cut mesh. This can be done by setting all $T_e$ to values in $2\mathbb{Z}$ (for why double integers are necessary, see [^Kaelberer_2007]). The solving is done by simple iterative rounding [^Bommes_2009]: we choose the value in $T_e$ which is the closest to a double integer, round and fix it, and repeat until all $T_e$ are rounded. Note that CoMISo [^Bommes_2012] has a more sophisticated Gauss-Seidel rounding algorithm that is more efficient; however, we include an implementation in Directional to avoid the dependency, and since this will be used for general $N$-function parameterization in future versions of Directional.

Mixed-integer parameterization is demonstrated in [Example 502]({{ repo_url }}/502_MixedIntegerParameterization/main.cpp). The essential difference from [Example 501]({{ repo_url }}/501_SeamlessParameterization/main.cpp) is by simply setting ```isInteger=true```, when passed to ```directional::parameterize()```.

*Note:* The input field ```horsers-cf.rawfield``` is computed according to [Example 303](#303-polycurl-reduction) to have negligible PolyCurl. As such, the rotationally-seamless parameterization has a very low error ($L_\infty$ of $1.25778\cdot 10^{-5})$). The rounding iterations incur some error, but it is rather low as well (after rounding $116$ variables it climbs to $L_\infty=0.657496$); we therefore recommend to warm-start a mixed-integer parameterization with a curl-reduced directional field.


![([Example 502]({{ repo_url }}/502_MixedIntegerParameterization/main.cpp)) Top Left: directional field. Top right: rotationally-seamless parameterization. Bottom: fully-seamless parameterization.](images/502_MixedIntegerParameterization.png)


## Outlook for continuing development

Directional is a budding project, and there are many algorithms in the state-of-the-art that we look forward to implement, with the help of volunteer researchers and practitioners from the field. Prominent examples of desired implementations are:

1. Face-based polar representation, and mixed-integer directional algorithms.

2. Support for 3D *Octahedral* fields [^solomon_2017], both in tet meshes and with the boundary-element method.

3. A discrete exterior calculus framework.

4. Differential operators and Hodge decomposition.

5. Cutting, integration, and parameterization. Note the libigl has this capacity that could be called from Directional, but they are not entirely compatible.

6. Support for tensor fields.

7. Advanced and better visualization techniques.

## References
[^azencot_2017]: Omri Azencot, Etienne Corman, Mirela Ben-Chen, Maks Ovsjanikov, [Consistent Functional Cross Field Design for Mesh Quadrangulation](http://www.cs.technion.ac.il/~mirela/publications/cfc.pdf), 2017.
[^bommes_2009]: David Bommes, Henrik Zimmer, Leif Kobbelt, [Mixed-integer quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf), 2009.
[^Bommes_2012]: David Bommes, Henrik Zimmer, Leif Kobbelt, [Practical Mixed-Integer Optimization for Geometry Processing](https://www.graphics.rwth-aachen.de/publication/0319/), 2012.
[^bouaziz_2012]: Sofien Bouaziz, Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly, [Shape-Up: Shaping Discrete Geometry with Projections](http://lgg.epfl.ch/publications/2012/shapeup.pdf), 2012.
[^crane_2010]: Keenan Crane, Mathieu Desbrun, Peter Schr&ouml;der, [Trivial Connections on Discrete Surfaces](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/), 2010.
[^diamanti_2014]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung, [Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/), 2014.
[^diamanti_2015]: Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung, [Integrable PolyVector Fields](http://igl.ethz.ch/projects/integrable/), 2015.
[^degoes_2016]: Fernando de Goes, Mathieu Desbrun, Yiying Tong, [Vector Field Processing on Triangle Meshes](http://geometry.caltech.edu/pubs/dGDT16.pdf), 2016.
[^Kaelberer_2007]: Felix K&auml;lberer, Matthias Nieser, Konrad Polthier, [QuadCover - Surface Parameterization using Branched Coverings](http://www.mi.fu-berlin.de/en/math/groups/ag-geom/publications/ressources/db/KNP07-QuadCover.pdf), 2007
[^knoppel_2013]: Felix Kn&ouml;ppel, Keenan Crane, Ulrich Pinkall, and Peter Schr&ouml;der, [Globally Optimal Direction Fields](http://www.cs.columbia.edu/~keenan/Projects/GloballyOptimalDirectionFields/paper.pdf), 2013.
[^ray_2008]: Nicolas Ray, Bruno Vallet, Wan Chiu Li, Bruno Lévy, [N-Symmetry Direction Field Design](http://alice.loria.fr/publications/papers/2008/DGF/NSDFD-TOG.pdf), 2008.
[^liu_2011]: Yang Liu, Weiwei Xu, Jun Wang, Lifeng Zhu, Baining Guo, Falai Chen, Guoping Wang, [General Planar Quadrilateral Mesh Design Using Conjugate Direction Field](http://research.microsoft.com/en-us/um/people/yangliu/publication/cdf.pdf), 2008.
[^Myles_2014]: Ashish Myles, Nico Pietroni, Denis Zorin, [Robust Field-aligned Global Parametrization](http://vcg.isti.cnr.it/Publications/2014/MPZ14/), 2014.
[^solomon_2017]: Justin Solomon, Amir Vaxman, David Bommes, [Boundary Element Octahedral Fields in Volumes](http://www.staff.science.uu.nl/~vaxma001/frames3d.pdf), 2017.
[^panozzo_2014]: Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung,  [Frame Fields: Anisotropic and Non-Orthogonal Cross Fields](http://cs.nyu.edu/~panozzo/papers/frame-fields-2014.pdf), 2014.
[^vaxman_2016]: Amir Vaxman, Marcel Campen, Olga Diamanti, Daniele Panozzo, David Bommes, Klaus Hildebrandt, Mirela Ben-Chen, [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis), 2016.
