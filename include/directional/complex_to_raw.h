#ifndef COMPLEX_TO_RAW_H
#define COMPLEX_TO_RAW_H
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/complex_to_representative.h>
#include <igl/local_basis.h>


namespace directional
{
	// Computes the raw vector field given a complex field.
	// Inputs:
	//  B1, B2:
	//  B3: #F normals for each face/B3 from igl::local_base.
	//  N: the degree of the field.
	//  complex: Representation of the field as complex double
	// Outputs:
	//  raw: #F by 3*N matrix with all N explicit vectors of each directional in the order X,Y,Z,X,Y,Z, ...
	IGL_INLINE void complex_to_raw(const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::MatrixXd& B3,
		const Eigen::MatrixXcd& complex,
		int N,
		Eigen::MatrixXd& raw)
	{
		Eigen::MatrixXd representative;
		complex_to_representative(B1, B2, complex, N, representative);
		representative_to_raw(B3, representative, N, raw);
	}

	// Computes the raw vector field given a complex field.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  adjustAngles: #E angles that encode deviation from parallel transport.
	//  complex: Representation of the field as complex doubles
	//  N: the degree of the field.
	// Outputs:
	//  raw: #F by 3*N matrix with all N explicit vectors of each directional in the order X,Y,Z,X,Y,Z, ...
	IGL_INLINE void complex_to_raw(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXcd& complex,
		int N,
		Eigen::MatrixXd& raw)
	{
		Eigen::MatrixXd B1, B2, B3;
		igl::local_basis(V, F, B1, B2, B3);
		complex_to_raw(B1, B2, B3, complex, N, raw);
	}
}

#endif
