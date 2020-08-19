#ifndef DIRECTIONAL_RAWFIELD_TO_COLUMNDIRECTIONAL_H
#define DIRECTIONAL_RAWFIELD_TO_COLUMNDIRECTIONAL_H
#include <Eigen/Eigen>
namespace directional
{
    /**
     * Reshapes a rawfield, given as the directionals in [x_0,y_0,z_0,x_1,y_1,z_1,...x_N,y_N,z_N] for N directionals per face in each row
     * into a single column where all vectors of directional 0 are stacked vertically in rows [0, 3 * F), directional 1 in [3 * F, 6 * F) etc.
     * Input:
     * - columnDirectional N * F * 3 x 1 matrix of directional field in column form
     * - N Number of directional fields
     * Output:
     * - rawField F x N * 3 marix containing raw field representation.
     */
    inline void rawfield_to_columndirectional(const Eigen::MatrixXd& rawField, int N, Eigen::VectorXd& columnDirectional)
    {
        columnDirectional.resize(rawField.rows() * rawField.cols());
        const int fCount = rawField.rows();
        for(int n = 0; n < N; ++n)
        {
            for(int f = 0; f < fCount; ++f)
            {
                columnDirectional(n * 3 * fCount + 3 * f) = rawField(f, n * 3);
                columnDirectional(n * 3 * fCount + 3 * f + 1) = rawField(f, n * 3 + 1);
                columnDirectional(n * 3 * fCount + 3 * f + 2) = rawField(f, n * 3 + 2);
            }
        }
    }
}
#endif