#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <directional/dynamic_visualization.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>

// Mesh
Eigen::MatrixXd VMesh, CMesh;
Eigen::MatrixXi FMesh;

// Noodles
directional::noodleData n_data;		//Contains all noodle data
int streamLengths = 5;				//The number of segments a noodle consists of
int MaxLifespan = 10;				//The lifespan of a noodle before it respawns

// Vector field
int N = 2;							//Degree of the vector field

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	if (!viewer.core.is_animating)
		return false;

	directional::update_noodles(n_data, VMesh, FMesh);

	viewer.data().set_mesh(n_data.VNoodles, n_data.FNoodles);
	viewer.data().set_colors(n_data.CNoodles);

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
	Eigen::MatrixXcd powerField;
	Eigen::MatrixXd raw;
	Eigen::VectorXi b;
	Eigen::MatrixXd bc;

	b.resize(0);
	bc.resize(0, 3);

	directional::power_field(VMesh, FMesh, b, bc, N, powerField);

	// Convert it to raw field
	directional::power_to_raw(VMesh, FMesh, powerField, N, raw, true);

	//Initialize noodles
	directional::initialize_noodles(n_data, VMesh, CMesh, FMesh, raw, streamLengths, N, MaxLifespan, 4.0);

	// Viewer Settings
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;

	//Create the imported mesh in the viewer
	viewer.data().set_mesh(VMesh, FMesh);
	viewer.data().set_colors(CMesh);
	viewer.data().show_lines = false;

	//Create the noodles
	viewer.append_mesh();
	viewer.data().set_mesh(n_data.VNoodles, n_data.FNoodles);
	viewer.data().set_colors(n_data.CNoodles);
	viewer.data().show_lines = false;

	cout <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}
