// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_GLYPH_ARROWS_RAW_H
#define DIRECTIONAL_GLYPH_ARROWS_RAW_H

#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <igl/avg_edge_length.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <directional/line_boxes.h>
#include <Eigen/Core>


namespace directional
{
    void IGL_INLINE glyph_arrow(
        const Eigen::Vector3d& localField, 
        const Eigen::Vector3d& normal, 
        const Eigen::Vector3d& baryCenter,
        double width,
        double lengthScaling,
        double arrowWidth,
        double arrowLength,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        int outputVertexOffset = 0,
        int outputFaceOffset = 0)
    {
        const int usedVerticesCount = 7;
        const int usedFacesCount = 5;
        if(V.rows() < outputVertexOffset + usedVerticesCount)
        {
            V.conservativeResize(outputVertexOffset + usedVerticesCount, V.cols());
        }
        if (F.rows() < outputFaceOffset + usedFacesCount)
        {
            F.conservativeResize(outputFaceOffset + usedFacesCount, F.cols());
        }
        const auto fieldLength = localField.norm();
        const auto normedField = localField / fieldLength;

        // Perpendicular direction to field, relative to face normal
        auto left = normal.cross(normedField).normalized(); // To be sure, normalize it


        const auto endPos = baryCenter + normedField * fieldLength * lengthScaling;
        const auto endPosLine = baryCenter + normedField * (fieldLength * lengthScaling-arrowLength);
        // Square line vertices: bottomleft, bottomright, topleft, topright
        V.row(outputVertexOffset) = baryCenter + width * 0.5 * left;
        V.row(outputVertexOffset + 1 ) = baryCenter - width * 0.5 * left;
        V.row(outputVertexOffset + 2) = endPosLine + width * 0.5 * left;
        V.row(outputVertexOffset + 3) = endPosLine - width * 0.5 * left;
        // Arrow head vertices: top, right and left
        V.row(outputVertexOffset + 4) = endPos;
        V.row(outputVertexOffset + 5) = endPosLine - arrowWidth * 0.5 * left; //Right side arrow vertex
        V.row(outputVertexOffset + 6) = endPosLine + arrowWidth * 0.5 * left;
        F.block(outputFaceOffset, 0, usedFacesCount, 3) <<
            // Line faces
            0, 1, 2,
            2, 1, 3,
            // Arrow faces
            2, 3, 4,
            4, 3, 5,
            2, 4, 6;
        // Might not work. Check this!
        F.block(outputFaceOffset, 0, usedFacesCount, 3).array() += outputVertexOffset;
    }
    void IGL_INLINE glyph_arrow_centered(
        const Eigen::VectorXd& localField,
        const Eigen::VectorXd& normal,
        const Eigen::VectorXd& baryCenter,
        double width,
        double lengthScaling,
        double arrowWidth,
        double arrowLength,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        int outputOffset = 0,
        int outputFaceOffset = 0)
    {
        // Offset the bary center to create a centered glyph
        const auto fieldLength = localField.norm() * lengthScaling;
        auto newBaryCenter = baryCenter - 0.5 * fieldLength * localField.normalized();
        // Create the glyph.
        glyph_arrow(
            localField,
            normal,
            newBaryCenter,
            width,
            lengthScaling,
            arrowWidth,
            arrowLength,
            V,
            F,
            outputOffset,
            outputFaceOffset
        );
    }
    
    // Creates mesh elements that comprise glyph drawing of a directional field.
    // Inputs:
    //  V:          #V X 3 vertex coordinates.
    //  F:          #F by 3 face vertex indices.
    //  rawField:   A directional field in raw xyzxyz form
    //  glyphColor: An array of either 1 by 3 color values for each vector, N by 3 colors for each individual directional or #F*N by 3 colours for each individual vector, ordered by #F times vector 1, followed by #F times vector 2 etc.
    //  width, length, height: of the glyphs depicting the directionals

    // Outputs:
    //  fieldV: The vertices of the field mesh
    //  fieldF: The faces of the field mesh
    //  fieldC: The colors of the field mesh

    void IGL_INLINE glyph_arrows_raw(const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        const Eigen::MatrixXd &rawField,
        const Eigen::MatrixXd &glyphColor,
        double width,
        double lengthScaling,
        double arrowWidth,
        double arrowLength,
        Eigen::MatrixXd &fieldV,
        Eigen::MatrixXi &fieldF,
        Eigen::MatrixXd &fieldC)
    {
        // Compute normals for all faces
        Eigen::MatrixXd normals;
        igl::per_face_normals(V, F, normals);
        const int N = rawField.cols() / 3;
        Eigen::MatrixXd barycenters;
        igl::barycenter(V, F, barycenters);
        barycenters += normals * 0.001;

        // See function glyph_arrow() why.
        const int vertsPerGlyph = 7;
        const int facesPerGlyph = 5;

        fieldV.resize(F.rows() * N * vertsPerGlyph, 3);
        fieldF.resize(F.rows() * N * facesPerGlyph, 3);
        // Color per face
        fieldC.resize(F.rows() * N * facesPerGlyph, 3);

        std::cout << "#F:" << F.rows() << ", #field: " << rawField.rows() << std::endl;

        for (int n = 0; n < N; ++n)
        {
            for (int f = 0; f < F.rows(); ++f)
            {
                // Need to cast to vector, otherwise implicit cast to Vector3d in glyph_arrow messes up the block.
                Eigen::RowVector3d localField = rawField.block(f, 3 * n, 1, 3);
                glyph_arrow(
                    localField,
                    normals.row(f),
                    barycenters.row(f),
                    width,
                    lengthScaling,
                    arrowWidth,
                    arrowLength,
                    fieldV,
                    fieldF,
                    (f + n * F.rows()) * vertsPerGlyph,
                    (f + n * F.rows()) * facesPerGlyph
                    );
            }
        }

        // Determine colors
        if(glyphColor.rows() == 1)
        {
            fieldC = glyphColor.replicate(fieldC.rows(), 1);
        }
        // Note: in the degenerate case where N == #F, this colors the glyphs
        else if(glyphColor.rows() == N)
        {
            for(int n = 0; n < N; ++n)
            {
                fieldC.block(n * F.rows() * facesPerGlyph, 0, F.rows() * facesPerGlyph, 3) = glyphColor.row(n).replicate(F.rows() * facesPerGlyph, 1);
            }
        }
        // Color for each directional
        else
        {
            // Throw if not correct?
            assert(glyphColor.rows() == F.rows());
            for (int n = 0; n < N; ++n)
            {
                const auto rowStart = n * F.rows()* facesPerGlyph;
                for (int f = 0; f < F.rows(); ++f)
                {
                    fieldC.block(f * facesPerGlyph + rowStart, 0, facesPerGlyph, 3) = 
                        glyphColor.block(f,3*n,1,3).replicate(facesPerGlyph, 1);
                }
            }
        }
    }
}

#endif
