# Directional Version Log

## Version 3.0 Changes

- Libigl and CGAL are not longer dependencies.
- The viewer is now based on PolyScope.
- Added finite elements fields and discrete exterior calculus forms. This are special examples of the cochain complexes module that can incorporate rich structures.
- The Mesher is now independent, but can use GMP to speed up, if installed.
- Added PolyVector customize algorithms, with the possibility to iterate implicit steps and projections. Exemplified by Ginzburg-Landau fields and integrable fields.
- Phased out (temporarily) subdivision fields and conjugate fields, to modernize the code to the new version.


## Version 2.0 Changes

- Paradigm shift of the entire library to work with polymoprhic discrete tangent bundles. All algorithms (if possible) converted to work with this paradigm directly.

## Version 1.8.0 Changes

- Extension of the PolyVector and power field implementation to include soft and hard partial alignment, and rotational symmetry objectives.

## Version 1.7.0 Changes

- Major overhaul of the viewer, where most operation are encapsulated in a new  ```DIrectionalViewer``` class, allowing properties such as isolines, vertex, face, and edge data to be associated with a mesh. User does not (need to) directly control shapes and colors, but rather only provide the geometric quantities.
- Many changes in the actual visualizations - thin arrows for glyphs, edge bars for seams and $N$-function isolines, and nicer spheres for singularities.
- Several bug fixes.

## Version 1.6.0 Changes

- The seamless parameterization package was entirely replaced with a seamless integration and meshing package.
- Several bug fixes.

## Version 1.5.0 Changes

The major changes are:

- Change in name from "libdirectional" to "Directional".
- Full compatibility to the libigl 2.0 release, including similar website build and tutorial cmake paradigms.
- Introducing seamless parameterization.
- Rendering the different component of field visualization (mesh, field, streamlines, glyphs, singularities) as separate meshes, to facilitate coding. As a result, the "appending" options of ``directional::line_cylinders()`` etc. were removed.
- Visualization was canonicalized using ```directional\visualization_schemes.h``` for a more homogeneous look.
- The tutorial portions imported from libigl (Conjugate fields, curl reduction, streamline tracing) are now fully compatible with Directional representation and data structures.
- Combing is now dependent on a given matching, and therefore there is no separate ```curl_combing()``` or  ```principal_combing()```.
- Functionality was renamed to better reflect its general application and target, rather then the name assigned by the relevant papers.

Former Name | New Name
--------|----------------------------------------------------------------------
Integrable PolyVectors   | PolyCurl Reduction
Trivial Connection | Index prescription
Globally Optimal | Power FIelds


## Version 1.0 Changes
Alpha version of Directional (then called "libdirectional"). Introducing the following functionality:

1. Glyph Drawing with singularities.
2. Trivial connection.
3. Globally optimal fields.
4. Polyvectors + integrable polyvectors.
5. Principal matching and combing.


