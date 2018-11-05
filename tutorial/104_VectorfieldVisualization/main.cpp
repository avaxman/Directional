#include <igl/barycenter.h>
#include <igl/edge_topology.h>
#include <igl/local_basis.h>
#include <igl/parula.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/sort_vectors_ccw.h>
#include <directional/streamlines.h>
#include <igl/opengl/glfw/Viewer.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/line_cylinders.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>

#include <directional/dynamic_visualization.h>

Eigen::MatrixXd VMesh, CMesh;
Eigen::MatrixXi FMesh;

Eigen::MatrixXcd powerField;
Eigen::MatrixXd raw;

directional::StreamlineData sl_data;		
directional::StreamlineState sl_state;		//The actual state of the streamlines

directional::StreamlineState sl_state0;		//The streamline state when the streamline start


int degree = 1;								// degree of the vector field

int streamLengths = 5;						//The number of segments a streamline consists off
int currentSegment = 0;						//The last segment of the streamline which is currently up to be replaced by the front runner

int MaxLifespan = 20;						//The lifespan of a streamline before it respawns
Eigen::VectorXi currentLifespan;			//The number off itterations each streamline have been alive 

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	if (!viewer.core.is_animating)
		return false;

	dynamic_visualization::update(viewer, sl_data, sl_state, sl_state0, VMesh, FMesh, streamLengths, currentSegment, currentLifespan, MaxLifespan);

	return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	if (key == ' ')
	{
		viewer.core.is_animating = !viewer.core.is_animating;
		return true;
	}
	return false;
}
		
int main(int argc, char *argv[])
{
	using namespace Eigen;
	using namespace std;

	igl::opengl::glfw::Viewer viewer;

	// Load a mesh in OFF format
	igl::readOFF(TUTORIAL_SHARED_PATH "/lion.off", VMesh, FMesh);

	// Create a Vector Field
	Eigen::VectorXi b;
	Eigen::MatrixXd bc;

	b.resize(1);
	b << 0;
	bc.resize(1, 3);
	bc << 1, 1, 1;

	directional::power_field(VMesh, FMesh, b, bc, degree, powerField);

	// Convert it to raw field
	directional::power_to_raw(VMesh, FMesh, powerField, degree, raw, true);
	directional::streamlines_init(VMesh, FMesh, raw, sl_data, sl_state, 0.9);

	//Create a color mask for the imported mesh
	dynamic_visualization::create_mask(raw, degree, CMesh);
	
	//Create the imported mesh in the viewer
	viewer.data().set_mesh(VMesh, FMesh);
	viewer.data().set_colors(CMesh);
	viewer.data().show_lines = false;

	// Viewer Settings
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;


	// Save the spawn points off the seeds
	sl_state0 = sl_state;

	//Initialize Streamlines
	dynamic_visualization::initialize(viewer, streamLengths, sl_data, sl_state, VMesh, FMesh, currentLifespan, MaxLifespan);

	cout <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}
