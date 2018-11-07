#ifndef DIRECTIONAL_DYNAMIC_VISUALIZATION_H
#define DIRECTIONAL_DYNAMIC_VISUALIZATION_H

#include <igl/igl_inline.h>
#include <igl/parula.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/line_cylinders.h>

namespace dynamic_visualization
{ 
	IGL_INLINE void create_mask(
		Eigen::MatrixXd VMesh,
		Eigen::MatrixXi FMesh,
		directional::StreamlineData &sl_data,
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

		//Color the mesh based on the effort it takes 
		Eigen::VectorXd effort;
		directional::principal_matching(VMesh, FMesh, sl_data.EV, sl_data.EF, sl_data.FE, rawField, sl_data.matching, effort);		
		effort = effort.cwiseAbs();

		color.resize(FMesh.rows());
		for (int i = 0; i < FMesh.rows(); i++) 
		{
			Eigen::Vector3i indices = sl_data.FE.row(i);
			double value = (effort(indices(0))*effort(indices(0))) + (effort(indices(1))*effort(indices(1))) + (effort(indices(2))*effort(indices(2)));
			color(i) = sqrt(value);
		}

		igl::jet(color, true, C);
	}

	IGL_INLINE void initialize(
		igl::opengl::glfw::Viewer &viewer, 
		const int streamLengths, 
		const int degree,
		directional::StreamlineData &sl_data, 
		directional::StreamlineState &sl_state, 
		directional::StreamlineState &sl_state0,
		Eigen::MatrixXd &VMesh, 
		Eigen::MatrixXi &FMesh, 
		Eigen::VectorXi &currentLifespan, 
		const int MaxLifespan
	) {
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
		directional::streamlines_init(VMesh, FMesh, raw, sl_data, sl_state, 1.5 / double(degree));

		//Create a color mask for the imported mesh
		Eigen::MatrixXd CMesh;
		dynamic_visualization::create_mask(VMesh, FMesh, sl_data, raw, CMesh, false);

		//Create the imported mesh in the viewer
		viewer.data().set_mesh(VMesh, FMesh);
		viewer.data().set_colors(CMesh);
		viewer.data().show_lines = false;


		// Save the spawn points off the seeds
		sl_state0 = sl_state;


		//Create a mesh for each part of the noodle
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

		//Give all noodles a random starting age
		currentLifespan.resize(sl_state.start_point.rows());
		for (int i = 0; i < sl_state.start_point.rows(); i++)
			currentLifespan(i) = rand() % (MaxLifespan / 2);
	}

	IGL_INLINE void update_itteration_values(
		int &currentSegment,
		const int streamLengths,
		directional::StreamlineState &sl_state,
		directional::StreamlineState &sl_state0,
		Eigen::VectorXi &currentLifespan,
		const int MaxLifespan
	) {

		//update current segment number
		if (currentSegment < streamLengths - 1)
			currentSegment++;
		else
			currentSegment = 0;

		//Lifespan update
		for (int i = 0; i < sl_state.start_point.rows(); i++) {

			if (currentLifespan(i) < MaxLifespan)
				currentLifespan(i)++;
			else {
				sl_state.start_point.row(i) = sl_state0.start_point.row(i);
				sl_state.end_point.row(i) = sl_state0.end_point.row(i);
				sl_state.current_direction(i) = sl_state0.current_direction(i);
				sl_state.current_face(i) = sl_state0.current_face(i);
				currentLifespan(i) = 0;
			}
		}
	}

	IGL_INLINE void shade_noodles(
		const int streamLengths, 
		int &currentSegment, 
		Eigen::MatrixXd CFieldNew, 
		igl::opengl::glfw::Viewer &viewer
	) {

		//Shade the tail of the streamline
		for (int i = 1; i < streamLengths; i++) {
			//select the next segment off the streamline
			int nextSegment = currentSegment + 1 + i;
			if (nextSegment > streamLengths)
				nextSegment -= (streamLengths);
			viewer.selected_data_index = nextSegment;

			//darken it the futher it is away from the head
			Eigen::MatrixXd segmentC = CFieldNew * ((1.0 / streamLengths)*(i + 1));
			segmentC.conservativeResize(segmentC.rows(), 4);
			segmentC.col(3).setConstant(1.0 / i);

			viewer.data().set_colors(segmentC);
		}
	}

	IGL_INLINE void update_noodles(
		directional::StreamlineData &sl_data, 
		directional::StreamlineState &sl_state, 
		Eigen::MatrixXd &VMesh, 
		Eigen::MatrixXi &FMesh, 
		igl::opengl::glfw::Viewer &viewer, 
		const int streamLengths, 
		int &currentSegment
	) {

		directional::streamlines_next(VMesh, FMesh, sl_data, sl_state);

		Eigen::RowVector3d color(1.0, 1.0, 1.0);

		Eigen::MatrixXd VFieldNew, CFieldNew;
		Eigen::MatrixXi FFieldNew;
		directional::line_cylinders(sl_state.start_point, sl_state.end_point, 0.0005, color.replicate(sl_state.start_point.rows(), 1), 4, VFieldNew, FFieldNew, CFieldNew);

		viewer.selected_data_index = currentSegment + 1;  //Select the last segment off the streamline
		viewer.data().clear();
		viewer.data().set_mesh(VFieldNew, FFieldNew);
		viewer.data().set_colors(CFieldNew);

		shade_noodles(streamLengths, currentSegment, CFieldNew, viewer);
	}

	IGL_INLINE void update(
		igl::opengl::glfw::Viewer &viewer, 
		directional::StreamlineData &sl_data, 
		directional::StreamlineState &sl_state, 
		directional::StreamlineState &sl_state0, 
		Eigen::MatrixXd &VMesh, 
		Eigen::MatrixXi &FMesh, 
		const int streamLengths, 
		int &currentSegment, 
		Eigen::VectorXi &currentLifespan,
		const int MaxLifespan
	) {

		update_noodles(sl_data, sl_state, VMesh, FMesh, viewer, streamLengths, currentSegment);
		update_itteration_values(currentSegment, streamLengths, sl_state, sl_state0, currentLifespan, MaxLifespan);
	}

}

#endif
