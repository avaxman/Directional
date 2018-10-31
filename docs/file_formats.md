Directional file formats
===================

### Raw Field (```.rawfield```)

This  file format describe a face-based field given in uncompressed xyz form as follows:

    [Degree]  [#faces]
    [x] [y] [z] [x][y][z]  .... N*Degree columns
     ... #faces rows
     
The convention is that the vectors are assumed to be ordered by the given index within each face. Mostly expected to be CCW ordering around the normal for several functions. Altering this might result in unexpected results. The indices are **not** assumed to be matching across neighboring faces, and are independently indexed per face.

### Singularities (```.sing```)

This file format describes prescribed singularities as follows:

    [Degree]  [#singularities]
    [vertex_i] [index_i]
    ... #singularities rows

Where ```index_i``` are integers, and the actual fractional index for vertex $i$ is $\frac{index_i}{Degree}$. The indices can be negative and are considered unbounded if need be.

### Matching (```.matching```)

This file format described the matching between each two neighboring faces $i$ and $j$ (across a dual edge) as follows:

    [Degree]  [#dual_edges]
    [face_i] [face_j] [matching_i]
    ... #dual_edges rows
    
  The order is generally compatible with the result of the field ```EF``` in ```igl::edge_topology()``` for the same mesh, but that is not a guarantee. That means the vector $k$ in face $i$ is matched to vector $k+matching_i$ (modulu $Degree$) in face $j$, and the rest in an order-preserving manner. See tutorial for details.
    
    



### Singularities and Matching

### Data structures used from libigl:

- [.dmat](./dmat) uncompressed ASCII/binary files for dense matrices
- [.off](http://wias-berlin.de/software/tetgen/fformats.off.html) Geomview's polyhedral file format
- [.obj](http://en.wikipedia.org/wiki/Wavefront_.obj_file#File_format) Wavefront object file format. Usually unsafe to assume anything more than vertex positions and triangle indices are supported
- [.png](https://en.wikipedia.org/wiki/Portable_Network_Graphics) Portable Network Graphics image file. IGLLIB (in the libiglpng extra) supports png image files via the [yimg](https://github.com/yig/yimg) library. Alpha channels and compression are supported.


