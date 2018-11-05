#ifndef DIRECTIONAL_DYNAMIC_VISUALIZATION_H
#define DIRECTIONAL_DYNAMIC_VISUALIZATION_H

#include <igl/igl_inline.h>
#include <igl/parula.h>

namespace dynamic_visualization
{ 
	IGL_INLINE void create_mask(Eigen::MatrixXd field, const int N, Eigen::MatrixXd &C) {
		
		Eigen::VectorXd scalars;
		scalars.resize(field.rows());
		scalars.setZero();
		scalars = field.rowwise().norm();
		igl::parula(-scalars, true, C);
	}

	IGL_INLINE void initialize(igl::opengl::glfw::Viewer &viewer, const int streamLengths, directional::StreamlineData &sl_data, directional::StreamlineState &sl_state, Eigen::MatrixXd &VMesh, Eigen::MatrixXi &FMesh, Eigen::VectorXi &currentLifespan, const int MaxLifespan) {
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
	}

	IGL_INLINE void update_itteration_values(
		int &currentSegment,
		const int streamLengths,
		directional::StreamlineState &sl_state,
		directional::StreamlineState &sl_state0,
		Eigen::VectorXi &currentLifespan,
		const int MaxLifespan
	) {

		//Update the itteration values
		//segment number
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

	IGL_INLINE void shade_noodles(const int streamLengths, int &currentSegment, Eigen::MatrixXd CFieldNew, igl::opengl::glfw::Viewer &viewer) {

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

	IGL_INLINE void update_noodles(directional::StreamlineData &sl_data, directional::StreamlineState &sl_state, Eigen::MatrixXd &VMesh, Eigen::MatrixXi &FMesh, igl::opengl::glfw::Viewer &viewer, const int streamLengths, int &currentSegment) {

		directional::streamlines_next(VMesh, FMesh, sl_data, sl_state);

		Eigen::RowVector3d color(1.0, 1.0, 1.0);

		Eigen::MatrixXd VFieldNew, CFieldNew;
		Eigen::MatrixXi FFieldNew;
		directional::line_cylinders(sl_state.start_point, sl_state.end_point, 0.0005, color.replicate(sl_state.start_point.rows(), 1) /*Eigen::MatrixXd::Constant(sl_state.start_point.rows(),3,1.0)*/, 4, VFieldNew, FFieldNew, CFieldNew);

		viewer.selected_data_index = currentSegment + 1;  //Select the last segment off the streamline
		viewer.data().clear();
		viewer.data().set_mesh(VFieldNew, FFieldNew);
		viewer.data().set_colors(CFieldNew);

		shade_noodles(streamLengths, currentSegment, CFieldNew, viewer);
	}

	IGL_INLINE void update(igl::opengl::glfw::Viewer &viewer, directional::StreamlineData &sl_data, directional::StreamlineState &sl_state, directional::StreamlineState &sl_state0, Eigen::MatrixXd &VMesh, Eigen::MatrixXi &FMesh, const int streamLengths, int &currentSegment, Eigen::VectorXi &currentLifespan,
		const int MaxLifespan) {

		update_noodles(sl_data, sl_state, VMesh, FMesh, viewer, streamLengths, currentSegment);
		update_itteration_values(currentSegment, streamLengths, sl_state, sl_state0, currentLifespan, MaxLifespan);
	}

}

#endif
