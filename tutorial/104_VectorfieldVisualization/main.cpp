#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <directional/dynamic_visualization.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


Eigen::MatrixXd VMesh;						//Vertices of the imported mesh
Eigen::MatrixXi FMesh;						//Faces of the imported mesh

dynamic_visualization::noodleData n_data;

int degree = 3;								//Degree of the vector field
int streamLengths = 5;						//The number of segments a noodle consists off
int MaxLifespan = 20;						//The lifespan of a noodle before it respawns

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	if (!viewer.core.is_animating)
		return false;

	dynamic_visualization::update(viewer, n_data, VMesh, FMesh);

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

	//Initialize noodles
  dynamic_visualization::initialize(viewer, n_data, VMesh, FMesh, streamLengths, degree, MaxLifespan, 1.0);

	// Viewer Settings
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;

	cout <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}
