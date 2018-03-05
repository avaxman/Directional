# libdirectional
libdirectional is a C++ library for creating, manipulating and drawing directional fields on 3D meshes, build on [libigl](www.github.com/libigl/libigl)<sup>[1](#fn1)</sup> and, in turn, [Eigen](http://eigen.tuxfamily.org/). The  libdirectional represents directional objects, that is sets of vectors, per face of the mesh. The content and the notations are based on the Eurographics 2016 star, and consequently SIGGRAPH Asia 2016/SIGGRAPH 2017 course on [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis)<sup>[2](#fn2)</sup>. Some visualization code is borrowed from the [libhedra](https://github.com/avaxman/libhedra) library<sup>[3](#fn3)</sup>. 



## Installation
libdirectional is a header-only library where each file generally includes one function. To use the library, simply add the _include_ directory to your include path and make sure libigl and its prerequisites are set up properly. After that you can include any files you need normally, using for example `#include <directional/trivial_connection.h>`.

To use the examples simply clone the repository using:
```git
git clone --recursive https://github.com/avaxman/libdirectional.git
```

## Tutorial

A tutorial that walks through the entire functionality of libdirectional is forthcoming---real soon! 

<!---
Then open a shell in the folder containing the example you wish to run, in the folder `examples' and call:
```shell
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

This should properly set up the example including all dependencies upon which you can build it using your favourite compiler. For windows, you should use `cmake-gui ..` and follow the instructions to create a compilable Visual Studio file.
-->

##Future Plans

The following functionality is still needed for libdirectional:
 
Face-based polar representation, and consequent mixed-integer algorithms.
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
  author = {Amir Vaxman and others},
  note = {https://github.com/avaxman/libdirectional},
  year = {2017},
}
```

## Copyright
2017 Amir Vaxman, Sam de Redelijkheid, Daniele Panozzo, Olga Diamanti, Olga Sorkine-Hornung, and others.

## References
<a name="fn1">1</a>: A. Jacobson and D. Panozzo and others, [libigl: A simple C++ geometry processing library](http://libigl.github.io/libigl/), 2016.<br>
<a name="fn2">2</a>: A. Vaxman et al., [Directional Field Synthesis, Design, and Processing](https://github.com/avaxman/DirectionalFieldSynthesis), 2016.<br>
<a name="fn3">3</a>: A. Vaxman et al., [libhedra](https://github.com/avaxman/libhedra), 2016.<br>
<a name="fn4">4</a>: K. Crane and M. Desbrun and P. Schr&ouml;der, [Trivial Connections on Discrete Surfaces](https://www.cs.cmu.edu/~kmcrane/Projects/TrivialConnections/), 2010.<br>
<a name="fn5">5</a>: O. Diamanti et al., [Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/n-polyvector-fields.pdf), 2014.<br>
<a name="fn6">6</a>: F. Knöppel et al., [Globally optimal direction fields](https://www.cs.cmu.edu/~kmcrane/Projects/GloballyOptimalDirectionFields/paper.pdf), 2013.




