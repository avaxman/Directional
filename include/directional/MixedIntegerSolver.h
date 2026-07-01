#ifndef DIRECTIONAL_MIXED_INTEGER_SOLVER_HEADER_FILE
#define DIRECTIONAL_MIXED_INTEGER_SOLVER_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>
#include "eliminate_constraints.h"

namespace directional {

enum class RoundingTypeEnum { ITERATIVE_ROUNDING, BRANCH_AND_BOUND };

static bool solveUnconstrainedLeastSquares(const Eigen::SparseMatrix<double>& A,
                                           const Eigen::SparseMatrix<double>& M,
                                           const Eigen::VectorXd& b,
                                           Eigen::VectorXd& x,
                                           bool verbose = false,
                                           double regEps = 1e-10)
{
    using namespace Eigen;
    SparseMatrix<double> H = A.transpose() * M * A;
    {
        SparseMatrix<double> eye(H.rows(), H.cols());
        eye.setIdentity();
        H += regEps * eye;
    }

    SimplicialLDLT<SparseMatrix<double>> ldlt;
    ldlt.analyzePattern(H);
    ldlt.factorize(H);

    if (ldlt.info() != Success) {
        if (verbose) std::cout << "[ULS] SimplicialLDLT factorization failed." << std::endl;
        return false;
    }

    VectorXd r = A.transpose() * M * b;
    x = ldlt.solve(r);
    return true;
}

class MixedIntegerSolver {
public:
    RoundingTypeEnum roundingType;

    Eigen::SparseMatrix<double> A, C, M;
    Eigen::VectorXd b;
    Eigen::VectorXd fixedValues;
    Eigen::VectorXi fixedMask;
    Eigen::VectorXi toRoundMask;

    int    numVars;
    bool   verbose;
    double integerTolerance;

    Eigen::VectorXd x;

    MixedIntegerSolver()
        : roundingType(RoundingTypeEnum::ITERATIVE_ROUNDING),
          numVars(0), verbose(false), integerTolerance(1e-10) {}
    ~MixedIntegerSolver() {}

    bool solve()
    {
        using namespace Eigen;
        using namespace std;

        if (numVars <= 0) return false;

        // ----------------------------------------------------------------
        // Phase 1: Combine C and Fixed Variables into a Unified Constraint Matrix
        // ----------------------------------------------------------------
        int numFixedRows = fixedMask.sum();
        SparseMatrix<double> C_combined;
        
        if (numFixedRows > 0) {
            vector<Triplet<double>> totalTriplets;
            
            // 1. Copy original C elements
            for (int k=0; k<C.outerSize(); ++k) {
                for (SparseMatrix<double>::InnerIterator it(C, k); it; ++it) {
                    totalTriplets.emplace_back(it.row(), it.col(), it.value());
                }
            }
            
            // 2. Append fixed rows as identity constraints: 1.0 * x_i = 0
            int currentRow = C.rows();
            for (int i = 0; i < numVars; i++) {
                if (fixedMask(i)) {
                    totalTriplets.emplace_back(currentRow, i, 1.0);
                    currentRow++;
                }
            }
            
            C_combined.resize(C.rows() + numFixedRows, numVars);
            C_combined.setFromTriplets(totalTriplets.begin(), totalTriplets.end());
        } else {
            C_combined = C;
        }

        // ----------------------------------------------------------------
        // Phase 2: Subspace Compression
        // ----------------------------------------------------------------
        SparseMatrix<double> P;
        VectorXi full2subset, subset2full, I_prime;

        if (C_combined.rows() > 0) {
            eliminate_constraints(C_combined, P, full2subset, subset2full, toRoundMask, I_prime, verbose);
        } else {
            P.resize(numVars, numVars);
            P.setIdentity();
            subset2full = VectorXi::LinSpaced(numVars, 0, numVars - 1);
            full2subset = subset2full;
            I_prime = toRoundMask;
        }

        int numSubVars = P.cols();

        // Subspace tracking states
        VectorXi subFixedMask = VectorXi::Zero(numSubVars);
        VectorXd subFixedValues = VectorXd::Zero(numSubVars);
        VectorXi subToRoundMask = I_prime;

        // ----------------------------------------------------------------
        // Phase 3: Iterative Rounding Optimization Loop
        // ----------------------------------------------------------------
        while (true)
        {
            const int nFreeSub = numSubVars - subFixedMask.sum();

            SparseMatrix<double> subVar2FreeSubMat(numSubVars, nFreeSub);
            {
                vector<Triplet<double>> trp;
                int col = 0;
                for (int s = 0; s < numSubVars; s++) {
                    if (!subFixedMask(s)) trp.emplace_back(s, col++, 1.0);
                }
                subVar2FreeSubMat.setFromTriplets(trp.begin(), trp.end());
            }

            SparseMatrix<double> totalProjection = P * subVar2FreeSubMat;
            VectorXd fullFixedContribution = P * subFixedValues;

            SparseMatrix<double> APart = A * totalProjection;
            VectorXd bPart = b - A * fullFixedContribution;

            VectorXd xFreeSub;
            if (!solveUnconstrainedLeastSquares(APart, M, bPart, xFreeSub, verbose)) {
                return false;
            }

            VectorXd xSub = subVar2FreeSubMat * xFreeSub + subFixedValues;
            
            // Pure mapping from subspace to full space
            x = P * xSub;

            if (subToRoundMask.sum() == 0) break;

            double minIntDiff = numeric_limits<double>::max();
            int minIntDiffIdx = -1;
            bool changedFixed = false;

            for (int s = 0; s < numSubVars; s++) {
                if (!subToRoundMask(s)) continue;
                if (subFixedMask(s)) {
                    subToRoundMask(s) = 0;
                    continue;
                }

                changedFixed = true;
                
                // Correctly measure drift on pure subspace integer coordinates
                const double subVal = xSub(s);
                const double intDiff = std::abs(subVal - std::round(subVal));

                if (intDiff < integerTolerance) {
                    subFixedMask(s) = 1;
                    subToRoundMask(s) = 0;
                    subFixedValues(s) = std::round(subVal);
                } else if (intDiff < minIntDiff) {
                    minIntDiff = intDiff;
                    minIntDiffIdx = s;
                }
            }

            if (minIntDiffIdx != -1) {
                // Correctly snap pure subspace coordinate
                const double subVal = xSub(minIntDiffIdx);
                subFixedMask(minIntDiffIdx) = 1;
                subToRoundMask(minIntDiffIdx) = 0;
                subFixedValues(minIntDiffIdx) = std::round(subVal);
            }

            if (!changedFixed) break;
        }

        if (verbose && C.rows() > 0) {
            VectorXd fullResid = C * x;
            cout << "[MIS] Final constraint residual norm (Cx): " << fullResid.norm() << endl;
        }

        return true;
    }
};

} // namespace directional

#endif
