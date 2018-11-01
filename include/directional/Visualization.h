#ifndef DIRECTIONAL_VISUALIZATION_H
#define DIRECTIONAL_VISUALIZATION_H

#include <igl/igl_inline.h>
#include <igl/parula.h>

namespace visualization
{
	IGL_INLINE void updateItt(
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

	IGL_INLINE void create_mask(Eigen::MatrixXd field, int N, Eigen::MatrixXd &C) {
		Eigen::VectorXd scalars;
		scalars.resize(field.rows());
		scalars.setZero();
		/*
		std::cout << field.row(3) << std::endl;
		for (int i = 0; i < N; i++) {
			for (int x = 0; x < field.rows(); x++) {
				scalars(x) += field(field.cols()*x + 3 * N) + field(field.cols()*x + 3 * N+1) + field(field.cols()*x + 3 * N+2);
			}

			std::cout << scalars(0) << std::endl;
			std::cout << scalars(1) << std::endl;
			std::cout << scalars(2) << std::endl;
			std::cout << scalars(3) << std::endl << std::endl;
		}
		std::cout << scalars(3) << std::endl << std::endl;
		//scalars /= N;
		*/
		for (int i = 0; i < N; i++) {
			scalars += ((field.col(i * 3)*field.col(i * 3)) + (field.col(i * 3 + 1)*field.col(i * 3 + 1)) + (field.col(i * 3 + 2)*field.col(i * 3 + 2)));
		}

		igl::parula(-scalars, true, C);
	}
}

#endif
