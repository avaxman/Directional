// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2026 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_MIXED_INTEGER_SOLVER_HEADER_FILE
#define DIRECTIONAL_MIXED_INTEGER_SOLVER_HEADER_FILE

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <limits>
#include <cmath>
#include <vector>

// Mixed integer solver for:
//   x = argmin |Ax - b|^2  (with metric M)
//   subject to Cx = d
//   with a subset of x constrained to integers (toRoundMask)
//   and a subset fixed to known values (fixedMask / fixedValues)
//
// FIX 1: Proper diagonal regularization instead of perturbing only (0,0).
// FIX 2: fixedValues is zero-initialized for unfixed variables before each
//         solve, so that C * fixedValues is always a valid partial offset.
// FIX 3: Constraint RHS uses the original C consistently (not the already-
//         reduced CPart), which is correct because fixedValues carries zeros
//         for free variables.
// FIX 4: Constraint feasibility is checked before and after each solve, and
//         the rank of the (reduced) constraint matrix is reported when
//         verbose is true, so contradictory or redundant constraints surface.

namespace directional {

enum class RoundingTypeEnum { ITERATIVE_ROUNDING, BRANCH_AND_BOUND };


// ----------------------------------------------------------------------------
// Constraint diagnostics
// ----------------------------------------------------------------------------

// Returns true if the system CPart * x_free = rhs appears feasible, i.e. the
// rank of [CPart | rhs] equals the rank of CPart.  When verbose is true it
// also prints the rank and a representative residual so you can see how
// ill-conditioned things are.
static bool checkConstraintFeasibility(const Eigen::MatrixXd& CPart,
                                       const Eigen::VectorXd& rhs,
                                       bool verbose,
                                       const std::string& tag = "")
{
    using namespace Eigen;

    FullPivLU<MatrixXd> lu(CPart);
    const int rankC   = lu.rank();
    const int rankAug = FullPivLU<MatrixXd>((MatrixXd(CPart.rows(), CPart.cols() + 1)
                                              << CPart, rhs).finished()).rank();

    if (verbose) {
        std::cout << "[Constraints" << (tag.empty() ? "" : " " + tag) << "] "
                  << "rows=" << CPart.rows()
                  << "  cols=" << CPart.cols()
                  << "  rank(C)=" << rankC
                  << "  rank([C|d])=" << rankAug;
        if (rankC < CPart.rows())
            std::cout << "  ** redundant constraints **";
        if (rankAug > rankC)
            std::cout << "  *** CONTRADICTORY CONSTRAINTS ***";
        std::cout << std::endl;
    }

    return (rankAug == rankC);  // feasible iff augmented rank did not increase
}


// ----------------------------------------------------------------------------
// Core constrained least-squares solve
// ----------------------------------------------------------------------------

// Solves:  min_{x} (Ax-b)^T M (Ax-b)   s.t.  Cx = d
// via a Schur-complement / block-elimination approach.
//
// Returns false if any factorization fails.
static bool solveConstrainedLeastSquares(const Eigen::SparseMatrix<double>& A,
                                         const Eigen::SparseMatrix<double>& M,
                                         const Eigen::VectorXd& b,
                                         const Eigen::SparseMatrix<double>& C,
                                         const Eigen::VectorXd& d,
                                         Eigen::VectorXd& x,
                                         bool verbose = false,
                                         double regEps = 1e-10)
{
    using namespace Eigen;

    // H = A^T M A
    SparseMatrix<double> H = A.transpose() * M * A;

    // FIX 1: proper diagonal regularization
    {
        SparseMatrix<double> eye(H.rows(), H.cols());
        eye.setIdentity();
        H += regEps * eye;
    }

    SimplicialLDLT<SparseMatrix<double>> ldlt;
    ldlt.analyzePattern(H);
    ldlt.factorize(H);

    if (ldlt.info() != Success) {
        if (verbose)
            std::cout << "[CLS] SimplicialLDLT factorization of H failed." << std::endl;
        return false;
    }

    const int m = C.rows();
    VectorXd r = A.transpose() * M * b;

    if (m == 0) {
        // No constraints — unconstrained solve
        x = ldlt.solve(r);
        return true;
    }

    // Build Schur complement  S = C H^{-1} C^T
    MatrixXd Ct  = MatrixXd(C.transpose());
    MatrixXd Y   = ldlt.solve(Ct);           // H^{-1} C^T  (block solve)
    MatrixXd S   = C * Y;                    // m x m, symmetric PSD

    VectorXd Hr       = ldlt.solve(r);
    VectorXd rhs_lam  = C * Hr - d;          // S lambda = C H^{-1} r - d

    // Diagnose the Schur complement for near-singularity
    if (verbose) {
        JacobiSVD<MatrixXd> svd(S);
        const auto& sv = svd.singularValues();
        std::cout << "[CLS] Schur complement size=" << m
                  << "  sigma_max=" << sv(0)
                  << "  sigma_min=" << sv(sv.size()-1)
                  << "  cond=" << sv(0)/sv(sv.size()-1)
                  << std::endl;
    }

    VectorXd lambda;
    LDLT<MatrixXd> sldlt(S);
    if (sldlt.info() == Success) {
        lambda = sldlt.solve(rhs_lam);
    } else {
        if (verbose)
            std::cout << "[CLS] Dense LDLT on Schur complement failed; "
                         "falling back to ColPivHouseholderQR." << std::endl;
        lambda = S.colPivHouseholderQr().solve(rhs_lam);
    }

    x = ldlt.solve(r - C.transpose() * lambda);

    // Constraint satisfaction check after solve
    if (verbose) {
        VectorXd resid = C * x - d;
        std::cout << "[CLS] Post-solve constraint residual norm: "
                  << resid.norm() << std::endl;
    }

    return true;
}


// ----------------------------------------------------------------------------
// MixedIntegerSolver
// ----------------------------------------------------------------------------

class MixedIntegerSolver {
public:
    RoundingTypeEnum roundingType;

    Eigen::SparseMatrix<double> A, C, M;
    Eigen::VectorXd b, d;          // objective rhs and constraint rhs (Cx = d)
    Eigen::VectorXd fixedValues;   // length numVars; entry i used only when fixedMask(i)==1
    Eigen::VectorXi fixedMask;     // 1 = variable is fixed
    Eigen::VectorXi toRoundMask;   // 1 = variable must be rounded to integer

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

        // ----------------------------------------------------------------
        // Sanity checks on dimensions
        // ----------------------------------------------------------------
        if (numVars <= 0) {
            if (verbose) cout << "[MIS] numVars not set." << endl;
            return false;
        }
        if (fixedValues.size() != numVars ||
            fixedMask.size()   != numVars ||
            toRoundMask.size() != numVars) {
            if (verbose) cout << "[MIS] fixedValues / fixedMask / toRoundMask size mismatch." << endl;
            return false;
        }

        // ----------------------------------------------------------------
        // FIX 2: ensure fixedValues is zero for unfixed variables so that
        //         C * fixedValues is always a valid partial-constraint offset.
        // ----------------------------------------------------------------
        for (int i = 0; i < numVars; i++)
            if (!fixedMask(i))
                fixedValues(i) = 0.0;

        if (verbose) {
            cout << "[MIS] Initial fixed variables : " << fixedMask.sum()  << endl;
            cout << "[MIS] Variables to round      : " << toRoundMask.sum() << endl;
        }

        // Use d if supplied, otherwise assume Cx = 0.
        VectorXd dVec = (d.size() == C.rows()) ? d : VectorXd::Zero(C.rows());

        // ----------------------------------------------------------------
        // Main iterative-rounding loop
        // ----------------------------------------------------------------
        while (true)
        {
            const int nFree = numVars - fixedMask.sum();
            if (verbose)
                cout << "[MIS] Solving with " << nFree << " free variables." << endl;

            // Build var2AllMat  (numVars x nFree injection matrix)
            SparseMatrix<double> var2AllMat(numVars, nFree);
            {
                vector<Triplet<double>> trp;
                trp.reserve(nFree);
                int col = 0;
                for (int i = 0; i < numVars; i++)
                    if (!fixedMask(i))
                        trp.emplace_back(i, col++, 1.0);
                var2AllMat.setFromTriplets(trp.begin(), trp.end());
            }

            // Reduced system matrices
            SparseMatrix<double> APart = A * var2AllMat;
            SparseMatrix<double> CPart = C * var2AllMat;

            // FIX 3: RHS adjustment — subtract contribution of fixed variables.
            // fixedValues(i) == 0 for free variables (FIX 2), so this is safe.
            VectorXd bPart = b - A * fixedValues;
            VectorXd dPart = dVec - C * fixedValues;

            // ---------------------------------------------------------------
            // Feasibility check on reduced constraint system (verbose mode)
            // ---------------------------------------------------------------
            if (verbose && CPart.rows() > 0) {
                checkConstraintFeasibility(MatrixXd(CPart), dPart, verbose, "reduced");
            }

            // Solve
            VectorXd xPart;
            if (!solveConstrainedLeastSquares(APart, M, bPart, CPart, dPart,
                                              xPart, verbose)) {
                if (verbose)
                    cout << "[MIS] solveConstrainedLeastSquares() failed!" << endl;
                return false;
            }

            // Recover full solution
            x = var2AllMat * xPart + fixedValues;

            // FIX 4: verify constraint satisfaction on the full solution
            if (verbose && C.rows() > 0) {
                VectorXd fullResid = C * x - dVec;
                cout << "[MIS] Full-solution constraint residual norm: "
                     << fullResid.norm() << endl;
                if (fullResid.norm() > 1e-6)
                    cout << "[MIS] WARNING: constraint residual is large." << endl;
            }

            if (toRoundMask.sum() == 0)
                break;

            // ---------------------------------------------------------------
            // Choose which variable to round next
            // ---------------------------------------------------------------
            double minIntDiff    = numeric_limits<double>::max();
            int    minIntDiffIdx = -1;
            bool   changedFixed  = false;

            for (int i = 0; i < numVars; i++) {
                if (!toRoundMask(i)) continue;

                if (fixedMask(i)) {
                    // Already fixed in a previous iteration
                    toRoundMask(i) = 0;
                    continue;
                }

                changedFixed = true;
                const double val     = x(i);
                const double intDiff = fabs(val - round(val));

                if (intDiff < integerTolerance) {
                    // Close enough — snap and fix
                    fixedMask(i)   = 1;
                    toRoundMask(i) = 0;
                    fixedValues(i) = round(val);
                    if (verbose)
                        cout << "[MIS] Snap-to-integer: var " << i
                             << "  val=" << val
                             << "  -> " << fixedValues(i) << endl;
                } else if (intDiff < minIntDiff) {
                    minIntDiff    = intDiff;
                    minIntDiffIdx = i;
                }
            }

            if (minIntDiffIdx != -1) {
                const double val    = x(minIntDiffIdx);
                const double intVal = round(val);
                fixedMask(minIntDiffIdx)   = 1;
                toRoundMask(minIntDiffIdx) = 0;
                fixedValues(minIntDiffIdx) = intVal;
                if (verbose)
                    cout << "[MIS] Round: var " << minIntDiffIdx
                         << "  val=" << val << "  -> " << intVal
                         << "  (intDiff=" << minIntDiff << ")" << endl;

                // After rounding, check whether the updated fixed-variable
                // offset is still compatible with the constraint row space.
                if (verbose && C.rows() > 0) {
                    // Temporarily build the new free set for the check
                    const int nFreeNext = numVars - fixedMask.sum();
                    SparseMatrix<double> v2a(numVars, nFreeNext);
                    {
                        vector<Triplet<double>> trp;
                        int col = 0;
                        for (int i = 0; i < numVars; i++)
                            if (!fixedMask(i))
                                trp.emplace_back(i, col++, 1.0);
                        v2a.setFromTriplets(trp.begin(), trp.end());
                    }
                    MatrixXd CNext  = MatrixXd(C * v2a);
                    VectorXd dNext  = dVec - C * fixedValues;
                    bool feasible = checkConstraintFeasibility(CNext, dNext, verbose,
                                                              "post-round");
                    if (!feasible)
                        cout << "[MIS] *** Constraint system became infeasible after "
                                "rounding var " << minIntDiffIdx << "! ***" << endl;
                }
            }

            if (!changedFixed)
                break;   // no free integer variables remain
        }

        return true;
    }
};

} // namespace directional

#endif // DIRECTIONAL_MIXED_INTEGER_SOLVER_HEADER_FILE
