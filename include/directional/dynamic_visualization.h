#ifndef DIRECTIONAL_DYNAMIC_VISUALIZATION_H
#define DIRECTIONAL_DYNAMIC_VISUALIZATION_H

#include <igl/igl_inline.h>
#include <igl/parula.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/line_cylinders.h>
#include <directional/streamlines.h>
#include <directional/principal_matching.h>

namespace dynamic_visualization
{ 
	struct noodleData {
		directional::StreamlineData sl_data;		//The data of the vectorfield through which the
		directional::StreamlineState sl_state;		//The actual state of the noodle
		directional::StreamlineState sl_state0;		//The state of the noodles when they
		int streamLengths;						//The number of segments a noodle consists off
		int currentSegment;						//The last segment of the noodle which is currently up to be replaced by the front runner
		int MaxLifespan;						//The lifespan of a noodle before it respawns
		int degree;
		Eigen::VectorXi currentLifespan;			//The number off itterations each noodle have been alive
	};
  
  IGL_INLINE bool effort_based_coloring(const Eigen::MatrixXd& VMesh, const Eigen::MatrixXi& FMesh, const dynamic_visualization::noodleData &n_data, const Eigen::MatrixXd& rawField, Eigen::MatrixXd& colors)
  {
    //Color the mesh based on the effort it takes
    Eigen::VectorXd effort;
    //directional::principal_matching(VMesh, FMesh, n_data.sl_data.EV, n_data.sl_data.EF, n_data.sl_data.FE, rawField, n_data.sl_data.matching, effort);

    Eigen::VectorXd scalarColor(FMesh.rows());
    for (int i = 0; i < FMesh.rows(); i++)
    {
      Eigen::Vector3i indices = n_data.sl_data.FE.row(i);
      double value = (n_data.sl_data.effort(indices(0))*n_data.sl_data.effort(indices(0))) +
      (n_data.sl_data.effort(indices(1))*n_data.sl_data.effort(indices(1))) + (n_data.sl_data.effort(indices(2))*n_data.sl_data.effort(indices(2)));
      scalarColor(i) = sqrt(value);
    }
    
    igl::jet(scalarColor, true, colors);
    return true;
  }

	/*IGL_INLINE void create_mask(
		Eigen::MatrixXd VMesh,
		Eigen::MatrixXi FMesh,
		dynamic_visualization::noodleData &n_data,
		Eigen::MatrixXd rawField,
		Eigen::MatrixXd &C,
		const std::function<bool(Eigen::VectorXd&)> userFunc
	) {
		Eigen::VectorXd color;

		//When a function is profided, run this
		if (userFunc) {
			userFunc(color);
			igl::jet(color, true, C);
			return;
		}

		
	}*/

	IGL_INLINE void initialize(
                             igl::opengl::glfw::Viewer &viewer,						//viewer
                             dynamic_visualization::noodleData &n_data,				//noodle data
                             const Eigen::MatrixXd &VMesh,									//verts of the mesh
                             const Eigen::MatrixXi &FMesh,									//Faces of the mesh
                             const int streamLengths,								//How many segments does a noodle consists of
                             const int degree,										//degree of the vectorfield
                             const int MaxLifespan,									//The number of frames a noodle is allowed to life
                             const double percentage,								//Changes the amount of noodles on the mesh
                             const std::function<bool(const Eigen::MatrixXd&, const Eigen::MatrixXi&,  const dynamic_visualization::noodleData &n_data, const Eigen::MatrixXd&, Eigen::MatrixXd&)> userFunc = effort_based_coloring) {
		n_data.streamLengths = streamLengths;
		n_data.degree = degree;
		n_data.MaxLifespan = MaxLifespan;
		n_data.currentSegment = 0;

		// Create a Vector Field
		Eigen::MatrixXcd powerField;
		Eigen::MatrixXd raw;
		Eigen::VectorXi b;
		Eigen::MatrixXd bc;

		b.resize(1);
		b << 0;
		bc.resize(1, 3);
		bc << 1, 1, 1;

		directional::power_field(VMesh, FMesh, b, bc, degree, powerField);

		// Convert it to raw field
		directional::power_to_raw(VMesh, FMesh, powerField, degree, raw, true);
		directional::streamlines_init(VMesh, FMesh, raw, n_data.sl_data, n_data.sl_state, percentage / double(degree));

		//Create a color mask for the imported mesh
		Eigen::MatrixXd CMesh;
		//dynamic_visualization::create_mask(VMesh, FMesh, n_data, raw, CMesh, userFunc);
    userFunc(VMesh, FMesh, n_data, raw,CMesh);

		//Create the imported mesh in the viewer
		viewer.data().set_mesh(VMesh, FMesh);
		viewer.data().set_colors(CMesh);
		viewer.data().show_lines = false;


		// Save the spawn points off the seeds
		n_data.sl_state0 = n_data.sl_state;


		//Create a mesh for each part of the noodle
		for (int i = 0; i < streamLengths; i++) {
			directional::streamlines_next(VMesh, FMesh, n_data.sl_data, n_data.sl_state);

			Eigen::RowVector3d color(1.0, 1.0, 1.0);

			Eigen::MatrixXd VFieldNew, CFieldNew;
			Eigen::MatrixXi FFieldNew;
			viewer.append_mesh();
			directional::line_cylinders(n_data.sl_state.start_point, n_data.sl_state.end_point, 0.0005, color.replicate(n_data.sl_state.start_point.rows(), 1), 4, VFieldNew, FFieldNew, CFieldNew);
			viewer.data().set_mesh(VFieldNew, FFieldNew);
			viewer.data().set_colors(CFieldNew);
			viewer.data().show_lines = false;
		}

		//Give all noodles a random starting age
		n_data.currentLifespan.resize(n_data.sl_state.start_point.rows());
		for (int i = 0; i < n_data.sl_state.start_point.rows(); i++)
			n_data.currentLifespan(i) = rand() % (MaxLifespan / 2);
	}

	IGL_INLINE void update_itteration_values(
		dynamic_visualization::noodleData &n_data
	) {

		//update current segment number
		if (n_data.currentSegment < n_data.streamLengths - 1)
			n_data.currentSegment++;
		else
			n_data.currentSegment = 0;

		//Lifespan update
		for (int i = 0; i < n_data.sl_state.start_point.rows(); i++) {

			if (n_data.currentLifespan(i) < n_data.MaxLifespan)
				n_data.currentLifespan(i)++;
			else {
				n_data.sl_state.start_point.row(i) = n_data.sl_state0.start_point.row(i);
				n_data.sl_state.end_point.row(i) = n_data.sl_state0.end_point.row(i);
				n_data.sl_state.current_direction(i) = n_data.sl_state0.current_direction(i);
				n_data.sl_state.current_face(i) = n_data.sl_state0.current_face(i);
				n_data.currentLifespan(i) = 0;
			}
		}
	}

	IGL_INLINE void shade_noodles(
		dynamic_visualization::noodleData &n_data,
		Eigen::MatrixXd CFieldNew, 
		igl::opengl::glfw::Viewer &viewer
	) {

		//Shade the tail of the streamline
		for (int i = 1; i < n_data.streamLengths; i++) {
			//select the next segment off the streamline
			int nextSegment = n_data.currentSegment + 1 + i;
			if (nextSegment > n_data.streamLengths)
				nextSegment -= (n_data.streamLengths);
			viewer.selected_data_index = nextSegment;

			//darken it the futher it is away from the head
			Eigen::MatrixXd segmentC = CFieldNew * ((1.0 / n_data.streamLengths)*(i + 1));
			segmentC.conservativeResize(segmentC.rows(), 4);
			segmentC.col(3).setConstant(1.0 / i);

			viewer.data().set_colors(segmentC);
		}
	}

	IGL_INLINE void update_noodles(
		dynamic_visualization::noodleData &n_data,
		Eigen::MatrixXd &VMesh, 
		Eigen::MatrixXi &FMesh, 
		igl::opengl::glfw::Viewer &viewer
	) {

		directional::streamlines_next(VMesh, FMesh, n_data.sl_data, n_data.sl_state);

		Eigen::RowVector3d color(1.0, 1.0, 1.0);

		Eigen::MatrixXd VFieldNew, CFieldNew;
		Eigen::MatrixXi FFieldNew;
		directional::line_cylinders(n_data.sl_state.start_point, n_data.sl_state.end_point, 0.0005, color.replicate(n_data.sl_state.start_point.rows(), 1), 4, VFieldNew, FFieldNew, CFieldNew);

		viewer.selected_data_index = n_data.currentSegment + 1;  //Select the last segment off the streamline
		viewer.data().clear();
		viewer.data().set_mesh(VFieldNew, FFieldNew);
		viewer.data().set_colors(CFieldNew);

		shade_noodles(n_data, CFieldNew, viewer);
	}

	IGL_INLINE void update(
		igl::opengl::glfw::Viewer &viewer, 
		dynamic_visualization::noodleData &n_data,
		Eigen::MatrixXd &VMesh, 
		Eigen::MatrixXi &FMesh
	) {

		update_noodles(n_data, VMesh, FMesh, viewer);
		update_itteration_values(n_data);
	}

}

#endif
