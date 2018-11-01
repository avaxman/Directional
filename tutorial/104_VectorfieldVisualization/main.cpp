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

#include <directional/Visualization.h>

// Mesh
Eigen::MatrixXd VMesh, CMesh;
Eigen::MatrixXi FMesh;

Eigen::MatrixXcd powerField;
Eigen::MatrixXd raw;

directional::StreamlineData sl_data;
directional::StreamlineState sl_state;

directional::StreamlineState sl_state0;		//The starting points off the streamlines

int degree = 1;								// degree of the vector field

int streamLengths = 5;						//The number of segments a streamline consists off
int currentSegment = 0;						//The last segment of the streamline which is currently up to be replaced by the front runner

int MaxLifespan = 20;						//The lifespan of a streamline before it respawns
Eigen::VectorXi currentLifespan;			//The number off itterations past since the streamlines spawned

bool pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	using namespace Eigen;
	using namespace std;

	if (!viewer.core.is_animating)
		return false;

	directional::streamlines_next(VMesh, FMesh, sl_data, sl_state);

	RowVector3d color(1.0, 1.0, 1.0);

	MatrixXd VFieldNew, CFieldNew;
	MatrixXi FFieldNew;
	directional::line_cylinders(sl_state.start_point, sl_state.end_point, 0.0005, color.replicate(sl_state.start_point.rows(), 1) /*Eigen::MatrixXd::Constant(sl_state.start_point.rows(),3,1.0)*/, 4, VFieldNew, FFieldNew, CFieldNew);

	viewer.selected_data_index = currentSegment + 1;  //Select the last segment off the streamline
	viewer.data().clear();
	viewer.data().set_mesh(VFieldNew, FFieldNew);
	viewer.data().set_colors(CFieldNew);

	//Shade the tail of the streamline
	for (int i = 1; i < streamLengths; i++) {
		//select the next segment off the streamline
		int nextSegment = currentSegment + 1 + i;
		if (nextSegment > streamLengths)
			nextSegment -= (streamLengths);
		viewer.selected_data_index = nextSegment;

		//darken it the futher it is away from the head
		MatrixXd segmentC = CFieldNew * ((1.0 / streamLengths)*(i+1));
		segmentC.conservativeResize(segmentC.rows(), 4);
		segmentC.col(3).setConstant(1.0/i);

		viewer.data().set_colors(segmentC);
	}

	visualization::updateItt(currentSegment, streamLengths, sl_state, sl_state0, currentLifespan, MaxLifespan);

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

	directional::streamlines_init(VMesh, FMesh, raw, sl_data, sl_state);

	//get the color matrix for the mesh
	visualization::create_mask(raw, degree, CMesh);
	
	//triangle mesh
	viewer.data().set_mesh(VMesh, FMesh);
	viewer.data().set_colors(CMesh*0.8);
	viewer.data().show_lines = false;

	// Viewer Settings
	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.core.is_animating = false;
	viewer.core.animation_max_fps = 30.;

	// Draw initial seeds on sample points
	sl_state0 = sl_state;
	directional::streamlines_next(VMesh, FMesh, sl_data, sl_state0);
	Eigen::MatrixXd v = sl_state0.end_point - sl_state0.start_point;
	v.rowwise().normalize();

	//streamline meshes, create 1 for each segment we trace
	for (int i = 0; i < streamLengths; i++) {
		directional::streamlines_next(VMesh, FMesh, sl_data, sl_state);

		Eigen::RowVector3d color(1.0, 1.0, 1.0);

		Eigen::MatrixXd VFieldNew, CFieldNew;
		Eigen::MatrixXi FFieldNew;
		viewer.append_mesh();
		directional::line_cylinders(sl_state.start_point, sl_state.end_point, 0.0005, color.replicate(sl_state.start_point.rows(), 1), 4, VFieldNew, FFieldNew, CFieldNew);
		viewer.data().set_mesh(VFieldNew, FFieldNew);
		viewer.data().set_colors(CFieldNew);
		viewer.data().show_lines = false;
	}

	//Give all streams a random starting age
	currentLifespan.resize(sl_state.start_point.rows());
	for (int i = 0; i < sl_state.start_point.rows(); i++)
		currentLifespan(i) = rand() % (MaxLifespan / 2);

	cout <<
		"Press [space] to toggle animation" << endl;
	viewer.launch();
}
