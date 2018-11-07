#include <igl/readOFF.h>
#include <directional/streamlines.h>
#include <igl/opengl/glfw/Viewer.h>
#include <directional/line_cylinders.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>

#include <directional/dynamic_visualization.h>

Eigen::MatrixXd VMesh;						//Vertices of the imported mesh
Eigen::MatrixXi FMesh;						//Faces of the imported mesh


directional::StreamlineData sl_data;		//The data of the vectorfield through which the
directional::StreamlineState sl_state;		//The actual state of the noodle
directional::StreamlineState sl_state0;		//The state of the noodles when they 


int degree = 1;								//Degree of the vector field
int streamLengths = 5;						//The number of segments a noodle consists off
int currentSegment = 0;						//The last segment of the noodle which is currently up to be replaced by the front runner
int MaxLifespan = 20;						//The lifespan of a noodle before it respawns
Eigen::VectorXi currentLifespan;			//The number off itterations each noodle have been alive 

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

	//Initialize noodles
	dynamic_visualization::initialize(viewer, streamLengths, degree, sl_data, sl_state, sl_state0, VMesh, FMesh, currentLifespan, MaxLifespan);

	// Viewer Settings
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;

	cout <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}
