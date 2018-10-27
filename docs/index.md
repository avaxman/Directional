# libdirectional - A Directional-Field Processing Library

<!---[![Build Status](https://travis-ci.org/libigl/libigl.svg?branch=master)](https://travis-ci.org/libigl/libigl)
[![Build status](https://ci.appveyor.com/api/projects/status/mf3t9rnhco0vhly8/branch/master?svg=true)](https://ci.appveyor.com/project/danielepanozzo/libigl-6hjk1/branch/master)
![](libigl-teaser.png)!--->

<https://github.com/avaxman/libdirectional/>

libdirectional is a C++ library for creating, manipulating and visualizing directional fields on 3D meshes, built as an extension to [libigl](www.github.com/libigl/libigl)<sup>[1](#fn1)</sup> on the foundations of [Eigen](http://eigen.tuxfamily.org/). libdirectional represents directional fields:  discrete sets of vectors on meshes. The content and the notations are based on the Eurographics 2016 star (adapted subsequently to SIGGRAPH Asia 2016/SIGGRAPH 2017 courses) on [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis)<sup>[2](#fn2)</sup>. Some visualization code is borrowed from the [libhedra](https://github.com/avaxman/libhedra) library. 

## Installation
libdirectional is a header-only library where each file generally includes one function. To use the library, simply add the _include_ directory to your include path and make sure libdirectional and its prerequisites are set up properly. After that you can include any files you need normally, using for example `#include <directional/index_prescription.h>`.

To get the library, simply clone the repository using:
```git
git clone --recursive https://github.com/avaxman/libdirectional.git
```

## Features
The current version is 1.5, entailing the following funcionality:

1. Representation of per-face directional fields of any given degree and symmetry.
2. Visualization using glyphs and streamline tracing.
3. Principal and curl-matching and combing for $N$-directional fields.
4. Computation of power fields of $N$-RoSy fields.
5. PolyVector fields.
6. Optimization for curl reduction.
7. Conjugate fields.
8. Prescription of singularity and generator indices.
9. Rotationally- and fully-seamless parameterization.

libdirectional is **a header-only library**. You do not need to compile anything to use,
just include directional headers (e.g. `#include <directional/index_prescription.h>`) and run.  Each
header file contains a single function (e.g. `igl/index_prescription.h` contains
`igl::index_prescription()`). 


## Tutorial
A [Tutorial](https://avaxman.github.io/libdirectional/tutorial/tutorial.html) that walks through the entire functionality of libdirectional is available. To compile it, go to the `tutorial` folder, open a shell and call:

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```
This should properly set up the tutorial project, with the individual chapters as subprojects, and create project makefiles upon which you can build it using your favourite compiler. For windows, you should use `cmake-gui ..` and follow the instructions to create a compilable Visual Studio file.

## Coding Guidelines and Tips

libdirectional inherits and follows the strict coding guidelines of libigl: please take a look [here](http://libigl.github.io/libigl/style-guidelines) before submitting your pull requests.


## How to Contribute

If you are interested in joining development, please fork the repository and submit a [pull request](https://help.github.com/articles/using-pull-requests/) with your changes.

## License
libdirectional is primarily [MPL2](http://www.mozilla.org/MPL/2.0/) licensed ([FAQ](http://www.mozilla.org/MPL/2.0/FAQ.html)). Some files contain third-party code under other licenses. 

## Attribution

If you use libdirectional in your academic projects, please cite the implemented papers directly, and/or the EG STAR 2016 when appropriate. To cite the library in general, you could use this BibTeX entry:

```
@misc{libdirectional,
title = {{libdirectional}: directional field synthesis, design, and processing,
author = {Amir Vaxman and others},
note = {https://github.com/avaxman/libdirectional},
year = {2017},
}
```

## Contact

libdirectional is led by [Amir Vaxman](http://www.staff.science.uu.nl/~vaxma001/). Please [contact me](mailto:avaxman@gmail.com) if you have questions or comments. For troubleshooting, please post an [issue](https://github.com/avaxman/libdirectional/issues) on github.

If you're using libirectional in your projects, quickly [drop me a note](mailto:avaxman@gmail.com). Tell me who you are and what you're using it for. This helps justify spending time maintaining this library!

## Future Plans

The following functionality is still in workds for libdirectional:

1. Other discretizations: discrete exterior calculus, vector-based fields.
2. 3D fields.
3. Discrete vector calculus: operators and Hodge decomposition.
4. Line-integral convolution visualization.
5. N-Symmetry seamless parameterization (for hexagonal and triangle remeshing).
6. Subdivision fields.

If you would like to suggest further topics, would like to collaborate in implementation, complain about bugs or ask questions, please address [Amir Vaxman](avaxman@gmail.com) (or open an issue in the repository).

## Copyright
2017 Amir Vaxman, Sam de Redelijkheid, Daniele Panozzo, Olga Diamanti, Olga Sorkine-Hornung, and others.

Please see individual files for appropriate copyright notices.
