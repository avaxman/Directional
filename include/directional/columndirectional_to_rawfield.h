#ifndef DIRECTIONAL_COLUMNDIRECTIONAL_TO_RAWFIELD_H
#define DIRECTIONAL_COLUMNDIRECTIONAL_TO_RAWFIELD_H
#include <Eigen/Eigen>
namespace directional
{
    /**
     * Reshapes a column directional, a single column where all vectors of directional 0 are stacked vertically in rows 
     * [0, 3 * F), directional 1 in [3 * F, 6 * F) etc., to a rawField that contains in each row all the directionals 
     * associated to a face.
     * Input:
     * - columnDirectional N * F * 3 x 1 matrix of directional field in column form
     * - N Number of directional fields
     * Output:
     * - rawField F x N * 3 marix containing raw field representation.
     */
    inline void columndirectional_to_rawfield(const Eigen::VectorXd& columnDirectional, int N, Eigen::MatrixXd& rawField)
    {
        rawField.resize(columnDirectional.size() / (N * 3), N * 3);
        const int fCount = rawField.rows();
        for(int n = 0; n < N; ++n)
        {
            for(int f = 0; f < fCount; ++f)
            {
                rawField(f, n * 3) = columnDirectional(n * 3 * fCount + 3 * f);
                rawField(f, n * 3 + 1) = columnDirectional(n * 3 * fCount + 3 * f + 1);
                rawField(f, n * 3 + 2) = columnDirectional(n * 3 * fCount + 3 * f + 2);
            }
        }
    }
}
#endif