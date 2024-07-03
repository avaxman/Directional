Directional file formats
===================

### Raw Field (```.rawfield```)

This  file format describe a raw field given in uncompressed xyz form as follows:

    [Degree]  [#tangent_spaces]
    [x] [y] [z] [x][y][z]  .... N*Degree columns
     ... #faces rows
     
The convention is that the vectors are assumed to be ordered by the given index within each tangent space. Mostly expected to be CCW ordering around the normal for several functions. Altering this might result in unexpected results. The indices are **not** assumed to be matching across neighboring tangent spaces, and are independently indexed per space. Note that this file format does not assume anything about the tangent bundle itself or any of its implementations---it can store, for instance, either face- or vertex-based fields.

### Singularities (```.sing```)

This file format describes prescribed singularities as follows:

    [Degree]  [#singularities]
    [vertex_i] [index_i]
    ... #singularities rows

Where ```index_i``` are integers, and the actual fractional index for vertex $i$ is $\frac{index_i}{Degree}$. The indices can be negative and are considered unbounded if need be. These singularities should match the local cycles on the designated tangent bundle.

### Matching (```.matching```)

This file format described the matching between each two neighboring tangent spaces $i$ and $j$ (across an edge on the tangent bundle: ([vertex_e], [vertex_f])) as follows:

    [Degree]  [#dual_edges]
    [face_i] [face_j] [vertex_e] [vertex_f] [matching_k]
    ... #dual_edges rows
    
  That means the vector $k$ in face $i$ is matched to vector $k+matching_i$ (modulu $Degree$) in face $j$, and the rest in an order-preserving manner. See tutorial for details.
    
### Data structures used from libigl

- [.dmat](./dmat) uncompressed ASCII/binary files for dense matrices
- [.off](https://wias-berlin.de/software/tetgen/fformats.off.html) Geomview's polyhedral file format
- [.obj](https://en.wikipedia.org/wiki/Wavefront_.obj_file#File_format) Wavefront object file format. Usually unsafe to assume anything more than vertex positions and triangle indices are supported
- [.png](https://en.wikipedia.org/wiki/Portable_Network_Graphics) Portable Network Graphics image file. IGLLIB (in the libiglpng extra) supports png image files via the [yimg](https://github.com/yig/yimg) library. Alpha channels and compression are supported.


