// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_SINGULARITY_SPHERES_H
#define DIRECTIONAL_SINGULARITY_SPHERES_H

#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/jet.h>
#include <directional/point_spheres.h>

namespace directional
{

    // Returns a list of faces, vertices and color values that can be used to draw singularities for non-zero index values.
    // Input:
    //    sources:        #s X 3 point coordinates of the location of the elements
    //    normals:        #n X 3 normals to the locations.
    //    N:              Degree of the field.
    //    avgScale:       Scale of the object (e.g., average edge length in meshes)
    //    singElements:   Singular elements out of "sources".
    //    singIndices:    Their indices
    //    singularityColors: 2*N x 3 colos per positive index in order [-N,..-1, 1, N]
    // Output:
    //    singV:          The vertices of the singularity spheres.
    //    singF:          The faces of the singularity spheres.
    //    singC:         The colors of the singularity spheres.
    void IGL_INLINE singularity_spheres(const Eigen::MatrixXd& sources,
                                        const Eigen::MatrixXd& normals,
                                        const int N,
                                        const double avgScale,
                                        const Eigen::VectorXi& singElements,
                                        const Eigen::VectorXi& singIndices,
                                        const Eigen::MatrixXd singularityColors,
                                        Eigen::MatrixXd& singV,
                                        Eigen::MatrixXi& singF,
                                        Eigen::MatrixXd& singC,
                                        const double radiusRatio)

    {

        Eigen::MatrixXd points(singElements.size(), 3);
        Eigen::MatrixXd pointNormals(singElements.size(), 3);
        Eigen::MatrixXd colors(singElements.size(), 3);
        Eigen::MatrixXd positiveColors=singularityColors.block(singularityColors.rows()/2,0,singularityColors.rows()/2,3);
        Eigen::MatrixXd negativeColors=singularityColors.block(0,0,singularityColors.rows()/2,3);

        for (int i = 0; i < singIndices.rows(); i++)
        {
            points.row(i) = sources.row(singElements(i));
            pointNormals.row(i) =normals.row(singElements(i));
            if (singIndices(i) > 0)
                colors.row(i) = positiveColors.row((singIndices(i)-1 > positiveColors.rows()-1 ? positiveColors.rows()-1  : singIndices(i)-1) );
            else if (singIndices(i)<0)
                colors.row(i) = negativeColors.row((negativeColors.rows()+singIndices(i) > 0 ? negativeColors.rows()+singIndices(i) : 0));
            else
                colors.row(i).setZero(); //this shouldn't have been input

        }
        double radius = radiusRatio*avgScale/5.0;
        directional::point_spheres(points, pointNormals, radius, colors, 8, singV, singF, singC);

    }
}

#endif
