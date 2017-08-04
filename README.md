# libdirectional
libdirectional is a C++ library for creating, manipulating and drawing directional fields on 3D meshes, build on [libigl](www.github.com/libigl/libigl)<sup>[1](#fn1)</sup> and, in turn, [Eigen](http://eigen.tuxfamily.org/). The  libdirectional represents directional objects, that is sets of vectors, per face of the mesh. The content and the notations are based on the Eurographics 2016 star, and consequently SIGGRAPH Asia 2016/SIGGRAPH 2017 course on [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis)<sup>[2](#fn2)</sup>. Some visualization code is borrowed from the [libhedra](https://github.com/avaxman/libhedra) library<sup>[3](#fn3)</sup>. 



## Installation
libdirectional is a header-only library where each file generally includes one function. To use the library, simply add the _include_ directory to your include path and make sure libigl and its prerequisites are set up properly. After that you can include any files you need normally, using for example `#include <directional/trivial_connection.h>`.

To use the examples simply clone the repository using:
```git
git clone --recursive https://github.com/avaxman/libdirectional.git
```

Then open a shell in the folder containing the example you wish to run, in the folder `examples' and call:
```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

This should properly set up the example including all dependencies upon which you can build it using your favourite compiler. For windows, you should use `cmake-gui ..` and follow the instructions to create a compilable Visual Studio file.




## Representations

libdirectional uses several different representations to describe directional fields. However, it currently only supports tangent face-based representations, where each face is a  single discrete tangent plane. $|F|$ is the number of faces in the mesh, where for each face one of the following is given. $E_I$ is the set of inner edges (adjacent to two triangles). We further use $N$ to describe the degree of the field (must be fixed for the entire field).

| Method            | Representation                                                                                                            |
|-------------------|---------------------------------------------------------------------------------------------------------------------------|
| **Raw**               | A $|F|\times 3N$ double matrix, representing an $1^N$-vector field (a directional with $N$ independent vectors in each face) in the form $X_1, Y_1, Z_1, X_2, Y_2, Z_2, \cdots X_N, Y_N, Z_N$. Vectors are assumed to be ordered in counter clockwise order in most libdirectional functions.|
| **Representative**    | A $|F| \times 3$ double matrix representing the first vector in a directional. Represents an $N$-vector field, known as an $N$-RoSy.|
| **Rotation Angles** | A $|E_I|$-sized double vector representing the rotation angle between directions (without magnitude information) on two neighbouring triangles, used in combination with a global rotation to uniquely define the field. The type may only encode $N$-direction fields.|
| **Complex Cartesian**           | A $|F|$-sized *complex* vector, representing an $N$ rosy as a single complex number $y=u^N$, where all the possible roots $u$ comprise an $N$-RoSy. The magnitude is also encoded this way, though it may be neglected.
| **PolyVector**        | A $|F| \times N$ complex double matrix, encoding the coefficients $X$ of a complex polynomial $f(z)=\sum_{i=0}^{N-1}{X_iz^i}$, which roots $u$ are an $1^N$-vector field. Every row is encoded as $X_{N-1},\cdots, X_0$, where $X_0$ is the free coefficient. In case where the field is an $N$-RoSy, all coefficients but $X_0$ are zero.

libdirectional provides a number of conversion functions to switch between different forms of representation. Each of the functions is of the form \<method 1>\_to\_\<method 2>, where \<method 1> and \<method 2> are the bold parts of the method name in the above table in lower cacse. e.g. `rotation_to_representative()`. Some possible combinations are given by composing two functions in sequence.

Not every combination is possible; for instance, it is not possible to convert between PolyVectors and rotation angles, which do not represent the same generality of directional fields. For $N$-RoSy fields, for instance, you will most likely work primarily with the Complex Cartesian, Representative, and Rotation Angles representation. Using the more explicit raw representation is done mostly for output (or drawing). 

## Trivial Connection

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

## Principal Matching and Singularities
When given a field in representative form, it is possible to devise the rotation angles $\delta_{ij}$ by the process of *principal matching*. Principal matching finds which vectors in face $f_i$ match those of face $f_j$, in order-preserving manner. The sum of rotation angles that each matching induces is known as the *effort* of the matching. Assuming only 

This is done by passing a representative or raw field into the `matching_energy()` function.

With the given rotation angles, it is also possible to devise the singularities on $1$-rings, using the `singularities()` method. Singularities can be calculated using either the principal matching  or the adjustment angles. Singularity calculation from principal matching is subject to aliasing, so unless calculated using the original adjustment angles you will most likely not obtain the same singularities as used to create the field. 

An illustration of the sampling problem is visualized in the *Singularities* example, which allows you to toggle between prescribed singularities, on which a trivial connection is computed, and the calculated singularities from a (blind) principal matching of the resulting field. It is possible to save fields generated with the *trivial_connection* example and load them in the *singularities* example.


## Complex Cartesian Fields

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

## Drawing fields
To draw a vector field first it must be converted into either its representative or raw representation. Then it can be passed to the `drawable_field()` method which will create a list of vertices, colors and faces to represent the field, finally those can be merged with the mesh data and passed to the viewer. If the viewer is set to take colors per vertex the `COLOR_PER_VERTEX` flag should be set.

There are several ways you can set the field color:<br>
The easiest way is to simply pass one color as a single row matrix. This will create one uniform ly colored field.<br>
The second option is to pass a matrix with \|F\| rows, each ro representing the color for the matching face.<br>
It is also possible to pass N\*\|F\| colors to give the color for every single vector on every single directional. In this case colors should be ordered to first give the color for every first vector on each face and only after that the color for every second vertex. e.g. F<sub>1</sub>1, F<sub>2</sub>1 ... F<sub>1</sub>2, F<sub>2</sub>2 ... F<sub>1</sub>N, F<sub>2</sub>N.<br>
Finally it is possible to pass one color per vertexafter setting the `PER_VECTOR_COLOR` flag.

By default the size of each vector is set to be related to the average edge length, as well as the length of the actual vector. The base length and with can be manually set if needed. If you want all vectors to be equal in size you can scale them by normalizing each vector in the field matrix. 

##Future Plans

The following functionality is still needed for libdirectional:
 
Face-based polar representation, and consequent mixed-integer algorithms.
Integrable PolyVector fields
Support for 3D octahedral fields.
Other discretizations: vertex-based, DEC.
Discrete vector calculus: operators and Hodge decomposition.
Poisson equations and parameterization.

If you would like to suggest further topics, would like to collaborate in implementation, complain about bugs or ask questions, please address [Amir Vaxman] (avaxman@gmail.com) (or open an issue in the repository)

##Acknowledge libdirectional

If you use libhedra in your academic projects, please cite the implemented papers directly, and/or the EG STAR 2016 when appropriate. To cite the library in general, you could use this BibTeX entry:

```
@misc{libdirectional,
  title = {{libdirectional}: directional field synthesis, design, and processing,
  author = {Amir Vaxman and Sam de Rijkenheid and others},
  note = {https://github.com/avaxman/libdirectional},
  year = {2017},
}
```

## References
<a name="fn1">1</a>: A. Jacobson and D. Panozzo and others, [libigl: A simple C++ geometry processing library](http://libigl.github.io/libigl/), 2016.<br>
<a name="fn2">2</a>: A. Vaxman et al., [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis), 2016.<br>
<a name="fn3">3</a>: A. Vaxman et al., [libhedra](https://github.com/avaxman/libhedra), 2016.<br>
<a name="fn4">4</a>: K. Crane and M. Desbrun and P. Schr&ouml;der, [Trivial Connections on Discrete Surfaces](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/), 2010.<br>
<a name="fn5">5</a>: O. Diamanti et al., [Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/n-polyvector-fields.pdf), 2014.<br>
<a name="fn6">6</a>: F. Knöppel et al., [Globally optimal direction fields](https://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/paper.pdf), 2013.




