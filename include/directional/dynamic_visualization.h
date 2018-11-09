// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>, 2018 Amir Vaxman <avaxman@gmail.com>, 2018 Lennert Sietsma <lennertsietsma@gmail.com>


#ifndef DIRECTIONAL_DYNAMIC_VISUALIZATION_H
#define DIRECTIONAL_DYNAMIC_VISUALIZATION_H

#include <igl/igl_inline.h>
#include <igl/parula.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/line_cylinders.h>
#include <directional/streamlines.h>

namespace directional
{ 
	struct noodleData {
		directional::StreamlineData sl_data;		//The data of the vectorfield through which the
		directional::StreamlineState sl_state;		//The actual state of the noodle
		directional::StreamlineState sl_state0;		//The state of the noodles when they
		int streamLengths;							//The number of segments a noodle consists off
		int currentSegment;							//The last segment of the noodle which is currently up to be replaced by the front runner
		int MaxLifespan;							//The lifespan of a noodle before it respawns
		int degree;
		Eigen::VectorXi currentLifespan;			//The number off itterations each noodle have been alive 
		Eigen::MatrixXd VNoodles;					//Vertices belonging to the noodles
		Eigen::MatrixXd CNoodles;					//Colors belonging to the noodles
		Eigen::MatrixXi FNoodles;					//Faces belonging to the noodles
	};

	IGL_INLINE void effort_based_coloring(const Eigen::MatrixXd& VMesh, const Eigen::MatrixXi& FMesh, const directional::noodleData &n_data, const Eigen::MatrixXd& rawField, Eigen::MatrixXd& colors)
	{
		//Color the mesh based on the effort it takes by calculating a scalar color per face

		Eigen::VectorXd scalarColor(FMesh.rows());
		for (int i = 0; i < FMesh.rows(); i++)
		{
			Eigen::Vector3i indices = n_data.sl_data.FE.row(i);
			double value = (n_data.sl_data.effort(indices(0))*n_data.sl_data.effort(indices(0))) +
				(n_data.sl_data.effort(indices(1))*n_data.sl_data.effort(indices(1))) + (n_data.sl_data.effort(indices(2))*n_data.sl_data.effort(indices(2)));
			scalarColor(i) = sqrt(value);
		}

		igl::parula(scalarColor, true, colors);
	}

	IGL_INLINE void initialize_noodles(
		directional::noodleData &n_data,				//noodle data
		Eigen::MatrixXd &VMesh,									//verts of the mesh
		Eigen::MatrixXd &CMesh,									//colors of the mesh
		Eigen::MatrixXi &FMesh,									//Faces of the mesh
		const int streamLengths,								//How many segments does a noodle consists of
		const int degree,										//degree of the vectorfield
		const int MaxLifespan,									//The number of frames a noodle is allowed to life
		const double percentage,								//Changes the amount of noodles on the mesh
		const std::function<void(const Eigen::MatrixXd&, const Eigen::MatrixXi&, const directional::noodleData &n_data, const Eigen::MatrixXd&, Eigen::MatrixXd&)> userFunc = effort_based_coloring
	) {
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
		userFunc(VMesh, FMesh, n_data, raw, CMesh);

		// Save the spawn points off the seeds
		n_data.sl_state0 = n_data.sl_state;

		//Create a mesh for each part of the noodle
		for (int i = 0; i < streamLengths; i++) 
		{
			Eigen::MatrixXd VNoodlesNew, CNoodlesNew;
			Eigen::MatrixXi FNoodlesNew;

			directional::streamlines_next(VMesh, FMesh, n_data.sl_data, n_data.sl_state);
			Eigen::RowVector3d color(1.0, 1.0, 1.0);
			directional::line_cylinders(n_data.sl_state.start_point, n_data.sl_state.end_point, 0.0005, color.replicate(n_data.sl_state.start_point.rows(), 1), 4, VNoodlesNew, FNoodlesNew, CNoodlesNew);
			
			Eigen::MatrixXd VNoodlesTemp(n_data.VNoodles.rows() + VNoodlesNew.rows(), VNoodlesNew.cols());
			VNoodlesTemp << n_data.VNoodles, VNoodlesNew;
			n_data.VNoodles = VNoodlesTemp;
			
			Eigen::MatrixXd CNoodlesTemp(n_data.CNoodles.rows() + CNoodlesNew.rows(), CNoodlesNew.cols());
			CNoodlesTemp << n_data.CNoodles, CNoodlesNew;
			n_data.CNoodles = CNoodlesTemp;

			Eigen::MatrixXi FNoodlesTemp(n_data.FNoodles.rows() + FNoodlesNew.rows(), FNoodlesNew.cols());
			FNoodlesTemp << n_data.FNoodles.array(), (FNoodlesNew.array() + (i*VNoodlesNew.rows()));
			n_data.FNoodles = FNoodlesTemp;
		}

		//Give all noodles a random starting age
		n_data.currentLifespan.resize(n_data.sl_state.start_point.rows());
		for (int i = 0; i < n_data.sl_state.start_point.rows(); i++)
			n_data.currentLifespan(i) = rand() % (MaxLifespan / 2);
	}

	IGL_INLINE void update_itteration_values(
		directional::noodleData &n_data
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

	IGL_INLINE void update_noodles(
		directional::noodleData &n_data,
		Eigen::MatrixXd &VMesh,
		Eigen::MatrixXi &FMesh
	) {
		//move the noodle 1 frame in time
		directional::streamlines_next(VMesh, FMesh, n_data.sl_data, n_data.sl_state);

		//default color of the tip of the noodle is white
		Eigen::RowVector3d color(1.0, 1.0, 1.0);

		//Create the new tip of the noodle
		Eigen::MatrixXd VNoodlesNew, CNoodlesNew;
		Eigen::MatrixXi FNoodlesNew;
		directional::line_cylinders(n_data.sl_state.start_point, n_data.sl_state.end_point, 0.0005, color.replicate(n_data.sl_state.start_point.rows(), 1), 4, VNoodlesNew, FNoodlesNew, CNoodlesNew);
 
		//Update the noodle data with the new first segment (faces remain unchanged)
		int noodleVertsPerSegment = n_data.VNoodles.rows() / n_data.streamLengths;
		int noodleColorPerSegment = n_data.CNoodles.rows() / n_data.streamLengths;
		n_data.VNoodles.block(n_data.currentSegment * noodleVertsPerSegment, 0, noodleVertsPerSegment, 3) = VNoodlesNew;
		n_data.CNoodles.block(n_data.currentSegment * noodleColorPerSegment, 0, noodleColorPerSegment, 3) = CNoodlesNew;

		//Shade the tail of the noodles
		n_data.CNoodles *= (1.0 - (1.0 / n_data.streamLengths));

		//This updates the neccesary itteration values
		update_itteration_values(n_data);
	}

}

#endif
