//
// Created by avaxman on 03/04/2024.
//

#ifndef DIRECTIONAL_DEFINITIONS_H
#define DIRECTIONAL_DEFINITIONS_H

namespace directional{
    constexpr double PI = 3.14159265358979311599796346854;

    inline bool segment_segment_intersection(const Eigen::RowVector2d& p,
                                             const Eigen::RowVector2d& u,
                                             const Eigen::RowVector2d& q,
                                             const Eigen::RowVector2d& v,
                                             double& t,
                                             double& s,
                                             const double tol = 1e-6) {

        double det = u(0) * v(1) - u(1) * v(0);

        if (fabs(det) < tol / 100.0)
            return false; // Parallel lines, no intersection

        t = ((q(0) - p(0)) * v(1) - (q(1) - p(1)) * v(0)) / det;
        s = ((q(0) - p(0)) * u(1) - (q(1) - p(1)) * u(0)) / det;

        return (t >= -tol && t <= 1+tol && s >= -tol && s <= 1+tol);
    }

}




#endif //DIRECTIONAL_DEFINITIONS_H
