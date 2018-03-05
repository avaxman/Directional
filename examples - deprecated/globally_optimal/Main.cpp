#include <iostream>
#include <directional/drawable_field.h>
#include <directional/complex_field.h>
#include <directional/complex_to_representative.h>
#include <directional/complex_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/dual_cycles.h>
#include <directional/get_indices.h>
#include <directional/draw_singularities.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include "Main.h"


Eigen::VectorXi cIDs;
Eigen::MatrixXi F, fieldF, meshF, singF;
Eigen::MatrixXd V, C, meshV, meshC, fieldV, fieldC, singV, singC, representative, cValues;
Eigen::MatrixXcd complex;
Eigen::SparseMatrix<double, Eigen::RowMajor> cycles;
igl::viewer::Viewer viewer;

Eigen::MatrixXd positiveIndices(4, 3),
negativeIndices(4, 3);

int N = 4;
bool drag = false;
bool normalized = false;
bool showSing = false;
bool onePressed = false;

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
	V.resize(VA.rows() + VB.rows(), VA.cols());
	V << VA, VB;
	F.resize(FA.rows() + FB.rows(), FA.cols());
	F << FA, (FB.array() + VA.rows());
}

void draw_field()
{
	// Calculate the field
	directional::complex_field(meshV, meshF, cIDs, cValues, N, complex);

	// Convert it so it can be drawn
	directional::complex_to_representative(meshV, meshF, complex, N, representative);

	// Normalize if wanted
	if (normalized)
		representative.rowwise().normalize();

	// Calculate the vectors, faces and colors of the field representation
	directional::drawable_field(meshV, meshF, representative, Eigen::RowVector3d(0, 0, 1), N, 0, fieldV, fieldF, fieldC);
	meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);

	for (int i = 0; i < cIDs.rows(); i++)
		meshC.row(cIDs(i)) = Eigen::RowVector3d(1, 0, 0);

	Eigen::VectorXi indices;
	Eigen::VectorXd effort;
	directional::principal_matching(meshV, meshF, representative, N, effort);

	directional::get_indices(meshV, meshF, cycles, effort, N, indices);

	directional::draw_singularities(meshV, indices, positiveIndices, negativeIndices, .007, singV, singF, singC);

	Eigen::MatrixXd a;
	Eigen::MatrixXi b;
	ConcatMeshes(meshV, meshF, fieldV, fieldF, a, b);
	if (showSing && singF.size())
	{
		ConcatMeshes(a, b, singV, singF, V, F);
		C.resize(F.rows(), 3);
		C << meshC, fieldC, singC;
	}
	else
	{
		V = a;
		F = b;
		C.resize(F.rows(), 3);
		C << meshC, fieldC;
	}

	// Update the viewer
	viewer.data.clear();
	viewer.data.set_face_based(true);
	viewer.data.set_mesh(V, F);
	viewer.data.set_colors(C);
}

bool key_up(igl::viewer::Viewer& viewer, int key, int modifiers)
{
    switch (key)
    {
        case '1': onePressed=false; break;
    }
    return true;
}


bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	int borders;
	switch (key)
	{
		// Toggle field drawing for easier rotation
            
        case '1': onePressed=true; break;
	case 'D':
		drag = !drag;
		break;

		// Toggle singularities
	case 'S':
		showSing = !showSing;
		draw_field();
		break;

	// Reset the constraints
	case 'R':
		cIDs.resize(0);
		cValues.resize(0, 6);
		draw_field();
		break;

	// Toggle normalization
	case 'N':
		normalized = !normalized;
		draw_field();
		break;
	}
	return true;
}

//Select vertices using the mouse
bool mouse_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
	if (drag || (key != 0 && key != 2) || !onePressed)
		return false;
	int fid;
	Eigen::Vector3d bc;

	// Cast a ray in the view direction starting from the mouse position
	double x = viewer.current_mouse_x;
	double y = viewer.core.viewport(3) - viewer.current_mouse_y;
	if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
		viewer.core.proj, viewer.core.viewport, meshV, meshF, fid, bc))
	{
		//Remove constraint
		if (key == 2)
		{
			int i;
			for (i = 0; i < cIDs.rows(); i++)
				if (cIDs(i) == fid)
					break;
			if (i == cIDs.rows())
				return false;
			cIDs(i) = cIDs(cIDs.size() - 1);
			cIDs.conservativeResize(cIDs.rows() - 1);
			cValues.row(i) = cValues.row(cValues.rows() - 1);
			cValues.conservativeResize(cValues.rows() - 1, 3);
			draw_field();
			return true;
		}

		int i;
		for (i = 0; i < cIDs.rows(); i++)
			if (cIDs(i) == fid)
				break;
		if (i == cIDs.rows())
		{
			cIDs.conservativeResize(cIDs.rows() + 1);
			cIDs(i) = fid;
			cValues.conservativeResize(cValues.rows() + 1, 3);
		}

		// Calculate direction from the center of the face to the mouse
		cValues.row(i) = 
			 (meshV.row(meshF(fid, 0)) * bc(0) + 
			 meshV.row(meshF(fid, 1)) * bc(1) + 
			 meshV.row(meshF(fid, 2)) * bc(2) - 
			(meshV.row(meshF(fid, 0)) + 
				meshV.row(meshF(fid, 1)) + 
				meshV.row(meshF(fid, 2))) / 3).normalized();
		draw_field();
		return true;
	}
	return false;
};

int main()
{
	viewer.callback_key_down = &key_down;
    viewer.callback_key_up = &key_up;
	viewer.callback_mouse_down = &mouse_down;
	std::cout <<
		"  R       Reset the constraints" << std::endl <<
		"  N       Toggle field normalization" << std::endl <<
		"  1+L-bttn  Place constraint pointing from the center of face to the cursor" << std::endl <<
		"  1+R-bttn  Remove constraint" << std::endl <<
		"  D       Enable/disable constraint placement" << std::endl <<
		"  S       Toggle singularities" << std::endl;

	// Set colors for Singularities
	positiveIndices << .25, 0, 0,
		.5, 0, 0,
		.75, 0, 0,
		1, 0, 0;

	negativeIndices << 0, .25, 0,
		0, .5, 0,
		0, .75, 0,
		0, 1, 0;

	// Load mesh
	igl::readOBJ("../../data/rocker-arm2500.obj", meshV, meshF);

	cIDs.resize(0);
	cValues.resize(0, 3);


	directional::dual_cycles(meshF, cycles);


	draw_field();
	viewer.launch();
}
