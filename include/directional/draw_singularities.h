#ifndef DRAW_SINGULARITIES
#define DRAW_SINGULARITIES
#include <cmath>
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <Eigen/Core>
#include <igl/avg_edge_length.h>


namespace directional
{
    // Returns a list of faces, vertices and color values that can be used to draw singularities for non-zero index values.
    // Input:
    //  V:              #V X 3 vertex coordinates.
    //  indices:        #V x 1 index (/N) per vertex (must be 0<index<N-1)
    //  positiveColors: N x 3 colos per positive index
    // negativeColors:  N x 3 colos per negative index
    //  r:              one singularity sphere radius
    // Output:
    //  singV:          The vertices of the singularity spheres.
    //  singF:          The faces of the singularity spheres.
    //  singC:          The colors of the singularity spheres.
    void IGL_INLINE draw_singularities(const Eigen::MatrixXd& V,
                                       const Eigen::VectorXi& indices,
                                       const Eigen::MatrixXd& positiveColors,
                                       const Eigen::MatrixXd& negativeColors,
                                       const double r,
                                       Eigen::MatrixXd &singV,
                                       Eigen::MatrixXi &singF,
                                       Eigen::MatrixXd &singC)
    {
        std::vector<int> singularities;
        for (int i = 0; i < V.rows(); i++)
        {
            if (indices(i))
                singularities.push_back(i);
        }
        
        Eigen::MatrixXd points(singularities.size(), 3);
        Eigen::MatrixXd colors(singularities.size(), 3);
        for (int i = 0; i < singularities.size(); i++)
        {
            points.row(i) = V.row(singularities[i]);
            if (indices(singularities[i]) > 0)
                colors.row(i) = positiveColors.row(std::min(std::round(indices(singularities[i]))-1, (double)positiveColors.rows()-1));
            else
                colors.row(i) = negativeColors.row(std::min(std::abs(std::round(indices(singularities[i])))-1, (double)negativeColors.rows()-1));
            
        }
        directional::point_spheres(points, r, colors, 16, false, singV, singF, singC);
    }
    
}

#endif
