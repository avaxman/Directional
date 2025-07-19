# Directional - A Directional-Field Processing Library

[Github repository](https://github.com/avaxman/Directional/)

Directional is a C++ library for creating, manipulating, and visualizing directional fields on 3D meshes. It is based on [Eigen](http://eigen.tuxfamily.org/). Directional represents directional fields:  discrete sets of vectors on meshes. Much of the content and the notations are based on the Eurographics 2016 star (adapted subsequently to SIGGRAPH Asia 2016/SIGGRAPH 2017 courses) on [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis). Some visualization code is borrowed from the [libhedra](https://github.com/avaxman/libhedra) library.

## New: version 3.0

This version marks the biggest change from previous versions, in which it stops depending on both [Libigl](https://libigl.github.io/) and [CGAL](https://www.cgal.org/), making it much more independent, lightweight, and easy to maintain. Directional now only depends directly on [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) and [PolyScope](https://polyscope.run/). Both are bundled in. There is an optional dependency on [GMP](https://gmplib.org/) just for the Mesher package, but the Mesher can work without it (albeit slower).

There are other considerable changes; check the [Release history](https://avaxman.github.io/Directional/RELEASE_HISTORY/).

## Installation
Directional is a header-only library where each file generally includes one function. To use the library, simply add the _include_ directory to your include path and make sure Directional and its prerequisites are set up properly. After that you can include any files you need normally, using for example `#include <directional/index_prescription.h>`.

To get the library, simply clone the repository using:
```git
git clone --recursive https://github.com/avaxman/Directional.git
```

## Features
The current version is 2.0, which represents a big paradigm shift from previous versions. Namely, the library uses polymorphic classes to represent abstract discrete tangent bundles, and contains their impelmentations of these interfaces by finite elements on meshes. The library further comprises the following features

1. Visualization of fields using (sparsified) glyphs and streamline tracing.
2. Principal and curl matching, and combing for $N$-directional fields.
3. Computation of power fields of $N$-RoSy fields.
4. PolyVector fields.
5. Optimization for curl reduction.
6. Prescription of singularity, generator, and boundary indices.
7. Rotationally- and fully-seamless integration of branched $N$-functions.
8. Meshing of the arrangement of $N$-functions isolines into polygonal meshes.
9. Finite-element fields and Discrete Exterior calculus forms.

Directional is **a header-only library**. You do not need to compile anything to use,
just include directional headers (e.g. `#include <directional/index_prescription.h>`) and run.  Each
header file contains a single function (e.g. `igl/index_prescription.h` contains
`igl::index_prescription()`). The library is only dependent on [Eigen]

## Tutorial
A [Tutorial](https://avaxman.github.io/Directional/tutorial/) that walks through the entire functionality of Directional is available. To compile it, go to the `tutorial` folder, open a shell and call:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```
This should properly set up the tutorial project, with the individual chapters as subprojects, and create project makefiles upon which you can build it using your favourite compiler. For windows, we recommend to use `cmake-gui ..` and follow the instructions to create a compilable Visual Studio file.

## How to Contribute

If you are interested in joining development, please fork the repository and submit a [pull request](https://help.github.com/articles/using-pull-requests/) with your changes.

## License
Directional is primarily [MPL2](https://www.mozilla.org/MPL/2.0/) licensed ([FAQ](https://www.mozilla.org/MPL/2.0/FAQ/)). Some files contain third-party code under other licenses.

## Attribution

If you use Directional in your academic projects, please cite the implemented papers directly, and/or the EG STAR 2016 when appropriate. To cite the library in general, you could use this BibTeX entry:

```
@misc{Directional,
author       = {Amir Vaxman and others},
title        = {Directional: A library for Directional Field Synthesis, Design, and Processing},
doi          = {10.5281/zenodo.3338174},
url          = {https://doi.org/10.5281/zenodo.3338174}
}
```

## Contact

Directional is led by [Amir Vaxman](https://avaxman.github.io/). Please [contact me](mailto:avaxman@gmail.com) if you have questions or comments. For troubleshooting, please post an [issue](https://github.com/avaxman/Directional/issues) on github.

If you're using Directional in your projects, [drop me a note](mailto:avaxman@gmail.com). Tell me who you are and what you're using it for. This helps justify spending time maintaining this library!

## Copyright
2017 Amir Vaxman, Sam de Redelijkheid, Daniele Panozzo, Olga Diamanti, Olga Sorkine-Hornung, Bram Custers and others.

Please see individual files for appropriate copyright notices.
