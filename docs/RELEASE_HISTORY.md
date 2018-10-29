# Directional version tracking


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


