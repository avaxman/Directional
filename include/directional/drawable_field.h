#ifndef DRAWABLE_FIELD
#define DRAWABLE_FIELD
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <directional/line_cylinders.h>
#include <Eigen/Core>
#include <igl/avg_edge_length.h>


namespace directional
{
	enum field_draw_flags
	{
		NONE = 0,
		// Set if the viewer is expected to get one color per vertex instead of per face.
		COLOR_PER_VERTEX = 1,
		// If set color will be interpreted as a list of colors for each vector in a directional instead.
		PER_VECTOR_COLOR = 2
	};

	// Returns a list of faces, vertices and colour values that can be used to draw a vector field 
	// on a mesh.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  field: Either a representative or a raw vector field.
	//  color: An array of either 1 by 3 color values for each vector, #F by 3 colors for each 
	//         individual directional or #F*N by 3 colours for each individual vector, ordered 
	//         by #F times vector 1, followed by #F times vector 2 etc.
	//         If PER_VECTOR_COLOR is set the colors will be interpreted as a list of colors for
	//         each vector in a directional instead.
	//  N: The degree of the field.
	//  width: The width of each vector in the vector field.
	//  length: The width of each vector in the vector field, respective to each vector's length.
	//  flags: Options:
	//         NONE: Colors are defined per face and vector colors are either defined for the full 
	//               mesh or per face.
	//         COLOR_PER_VERTEX: Set if the viewer is expected to get one color per vertex instead 
	//                           of per face.
	//         PER_VECTOR_COLOR: If set color will be interpreted as a list of colors for each 
	//                           vector in a directional instead.
	// Outputs:
	//  fieldV: The vertices of the field.
	//  fieldF: The faces of the field.
	//  fieldC: The colors of  the field.
	void IGL_INLINE drawable_field(const Eigen::MatrixXd &V, 
		const Eigen::MatrixXi &F, 
		const Eigen::MatrixXd &field,
		const Eigen::MatrixXd &color, 
		int N, 
		double width,
		double length,
		int flags,
		Eigen::MatrixXd &fieldV, 
		Eigen::MatrixXi &fieldF, 
		Eigen::MatrixXd &fieldC)
	{
		Eigen::MatrixXd normals;
		igl::per_face_normals(V, F, normals);

		Eigen::MatrixXd rawField, barycenters, vectorColors, P1, P2;
		igl::barycenter(V, F, barycenters);

		P1.resize(F.rows() * N, 3);
		P2.resize(F.rows() * N, 3);
		vectorColors.resize(F.rows() * N, 3);
		
		if (field.cols() == 3)
			representative_to_raw(normals, field, N, rawField);
		else
			rawField = field;

		normals.array() *= width;
		barycenters += normals;
		P1 = barycenters.replicate(N, 1);

		for (int i = 0; i < N; i++)
			P2.middleRows(F.rows()*i, F.rows()) = rawField.middleCols(3*i, 3);

		P2.array() *= length;
		P2 += P1;
		
		// Duplicate colors so each cylinder gets the proper color
		if (flags & PER_VECTOR_COLOR)
		{
			vectorColors.resize(F.rows()*N, color.cols());

			for (int n = 0; n < N; n++)
			{
				if (n >= color.rows())
				{
					vectorColors.block(F.rows()*n, 0, F.rows()*N - color.rows(), vectorColors.cols()) = Eigen::MatrixXd::Zero(F.rows()*N - color.rows(), vectorColors.cols());
					break;
				}
				vectorColors.block(F.rows()*n, 0, F.rows(), vectorColors.cols()) = color.row(n).replicate(0, F.rows());
			}
		}
		else
		{
			if (color.rows() == 1)
				vectorColors = color.replicate(P1.rows(), 1);
			else if (color.rows() == F.rows())
				vectorColors = color.replicate(N, 1);
			else
				vectorColors = color;
		}

		Eigen::MatrixXd Vc, Cc, Vs, Cs;
		Eigen::MatrixXi Fc, Fs;
		// Draw cylinders
		line_cylinders(P1, P2, width, vectorColors, 6, flags & COLOR_PER_VERTEX, Vc, Fc, Cc);
		
		// Draw sphere over intersection
		point_spheres(barycenters, width*2, vectorColors.topRows(barycenters.rows()), 10, flags & COLOR_PER_VERTEX, Vs, Fs, Cs);

		Fs.array() += Fc.rows();

		// Merge
		fieldV.resize(Vc.rows() + Vs.rows(), Vc.cols());
		fieldC.resize(Cc.rows() + Cs.rows(), Cc.cols());
		fieldF.resize(Fc.rows() + Fs.rows(), Fc.cols());
		fieldV << Vc, Vs;
		fieldC << Cc, Cs;
		fieldF << Fc, Fs;
	}

	// Returns a list of faces, vertices and colour values that can be used to draw a vector field
	// on a mesh.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  field: Either a representative or a raw vector field.
	//  color: An array of either 1 by 3 color values for each vector, #F by 3 colors for each 
	//         individual directional or #F*N by 3 colours for each individual vector, ordered 
	//         by #F times vector 1, followed by #F times vector 2 etc.
	//  N: The degree of the field.
	//  flags: Options, can be combined using the bitwise | or operator:
	//         NONE: Colors are defined per face and vector colors are either defined for the full 
	//               mesh or per face.
	//         COLOR_PER_VERTEX: Set if the viewer is expected to get one color per vertex instead 
	//                           of per face.
	//         PER_VECTOR_COLOR: If set color will be interpreted as a list of colors for each 
	//                           vector in a directional instead.
	// Outputs:
	//  fieldV: The vertices of the field.
	//  fieldF: The faces of the field.
	//  fieldC: The colors of  the field.
	void IGL_INLINE drawable_field(const Eigen::MatrixXd &V,
		const Eigen::MatrixXi &F,
		const Eigen::MatrixXd &field,
		const Eigen::MatrixXd &color,
		int N,
		int flags,
		Eigen::MatrixXd &fieldV,
		Eigen::MatrixXi &fieldF,
		Eigen::MatrixXd &fieldC)
	{
		double l = igl::avg_edge_length(V, F);
		drawable_field(V, F, field, color, N, l/50, l/5, flags, fieldV, fieldF, fieldC);
	}

}

#endif