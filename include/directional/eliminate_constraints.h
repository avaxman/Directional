#pragma once
#include <Eigen/Sparse>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

inline void eliminate_constraints(
    const Eigen::SparseMatrix<double>& C,
    Eigen::SparseMatrix<double>&       P,
    Eigen::VectorXi&                   full2subset,
    Eigen::VectorXi&                   subset2full,
    const Eigen::VectorXi&             preferSurvive, // This acts as input 'I' (1 for integer, 0 for continuous)
    Eigen::VectorXi&                   I_prime,       // Output: The modified 'I' for the x_subset vector
    bool                               verbose = false)
{
    const int n = C.cols();
    const int m = C.rows();
    const double eps = 1e-10;
    const bool hasPreference = (preferSurvive.size() == n);

    if (verbose) {
        std::cout << "[EC] Running homogeneous integer-safe elimination: m=" << m << " n=" << n << "\n";
    }

    // expr[i] maps: original_var[i] -> sum(coeff * current_active_var)
    std::vector<std::map<int, double>> expr(n);
    for (int i = 0; i < n; i++)
        expr[i][i] = 1.0;

    std::vector<bool> eliminated(n, false);
    std::vector<bool> active_is_integer(n, false);
    if (hasPreference) {
        for(int i = 0; i < n; ++i) {
            if (preferSurvive(i)) active_is_integer[i] = true;
        }
    }

    // Load constraints into row-major structures
    std::vector<std::map<int, double>> rows(m);
    {
        Eigen::SparseMatrix<double, Eigen::RowMajor> Crow(C);
        for (int r = 0; r < m; r++) {
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Crow, r); it; ++it) {
                if (std::abs(it.value()) > eps)
                    rows[r][it.col()] = it.value();
            }
        }
    }

    // =========================================================================
    // PHASE 1: Eliminate Continuous Variables
    // =========================================================================
    for (int r = 0; r < m; r++)
    {
        std::map<int, double> active_row;
        for (auto& [col, val] : rows[r]) {
            for (auto& [active_col, active_val] : expr[col]) {
                active_row[active_col] += val * active_val;
                if (std::abs(active_row[active_col]) < eps)
                    active_row.erase(active_col);
            }
        }

        if (active_row.empty()) continue; // Redundant row (0=0), skip safely.

        int pivot_col = -1;
        double pivot_val = 0.0;
        for (auto& [col, val] : active_row) {
            if (eliminated[col]) continue;
            if (!active_is_integer[col]) {
                if (std::abs(val) > std::abs(pivot_val)) {
                    pivot_col = col;
                    pivot_val = val;
                }
            }
        }

        if (pivot_col == -1) continue;

        eliminated[pivot_col] = true;

        std::map<int, double> pivot_expr;
        for (auto& [col, val] : active_row) {
            if (col == pivot_col || std::abs(val) < eps) continue;
            pivot_expr[col] = -val / pivot_val;
        }

        for (int i = 0; i < n; i++) {
            auto it = expr[i].find(pivot_col);
            if (it == expr[i].end()) continue;
            double coeff = it->second;
            expr[i].erase(it);
            for (auto& [k, v] : pivot_expr) {
                expr[i][k] += coeff * v;
                if (std::abs(expr[i][k]) < eps) expr[i].erase(k);
            }
        }
    }

    // =========================================================================
    // PHASE 2 & 3: Integer Variable Elimination
    // =========================================================================
    for (int r = 0; r < m; r++)
    {
        while (true) {
            std::map<int, double> active_row;
            for (auto& [col, val] : rows[r]) {
                for (auto& [active_col, active_val] : expr[col]) {
                    active_row[active_col] += val * active_val;
                    if (std::abs(active_row[active_col]) < eps)
                        active_row.erase(active_col);
                }
            }

            if (active_row.empty()) break;

            int gcd_pivot_col = -1;
            double gcd_pivot_val = 0.0;

            for (auto& [col, val] : active_row) {
                assert(active_is_integer[col]);
                bool divides_all = true;
                for (auto& [other_col, other_val] : active_row) {
                    double remainder = std::fmod(other_val, val);
                    if (std::abs(remainder) > 1e-5 && std::abs(std::abs(remainder) - std::abs(val)) > 1e-5) {
                        divides_all = false;
                        break;
                    }
                }
                if (divides_all) {
                    gcd_pivot_col = col;
                    gcd_pivot_val = val;
                    break;
                }
            }

            if (gcd_pivot_col != -1) {
                eliminated[gcd_pivot_col] = true;

                std::map<int, double> pivot_expr;
                for (auto& [col, val] : active_row) {
                    if (col == gcd_pivot_col) continue;
                    pivot_expr[col] = -val / gcd_pivot_val;
                }

                for (int i = 0; i < n; i++) {
                    auto it = expr[i].find(gcd_pivot_col);
                    if (it == expr[i].end()) continue;
                    double coeff = it->second;
                    expr[i].erase(it);
                    for (auto& [k, v] : pivot_expr) {
                        expr[i][k] += coeff * v;
                        if (std::abs(expr[i][k]) < eps) expr[i].erase(k);
                    }
                }
                break;
            }

            int col_A = -1, col_B = -1;
            double val_A = std::numeric_limits<double>::max();
            double val_B = std::numeric_limits<double>::max();
            for (auto& [col, val] : active_row) {
                double abs_val = std::abs(val);
                if (abs_val < val_A) {
                    col_B = col_A; val_B = val_A;
                    col_A = col;   val_A = abs_val;
                } else if (abs_val < val_B) {
                    col_B = col;   val_B = abs_val;
                }
            }

            if (col_A == -1 || col_B == -1) break;

            double coeff_A = active_row[col_A];
            double coeff_B = active_row[col_B];
            int q = static_cast<int>(std::round(coeff_B / coeff_A));
            if (q == 0) q = (coeff_B > 0 ? 1 : -1);

            for (int i = 0; i < n; i++) {
                auto it_B = expr[i].find(col_B);
                if (it_B != expr[i].end()) {
                    double val_B_in_expr = it_B->second;
                    expr[i][col_A] += q * val_B_in_expr;
                    if (std::abs(expr[i][col_A]) < eps) expr[i].erase(col_A);
                }
            }
        }
    }

    // =========================================================================
    // PHASE 4: Extract Output Spaces & Matrix Construction
    // =========================================================================
    full2subset.setConstant(n, -1);
    std::vector<int> s2f;
    s2f.reserve(n);

    for (int i = 0; i < n; i++) {
        if (!eliminated[i]) {
            full2subset[i] = (int)s2f.size();
            s2f.push_back(i);
        }
    }

    const int k = (int)s2f.size();
    subset2full.resize(k);
    I_prime.resize(k);

    for (int s = 0; s < k; s++) {
        int full_idx = s2f[s];
        subset2full[s] = full_idx;
        I_prime(s) = active_is_integer[full_idx] ? 1 : 0;
    }

    if (verbose) {
        std::cout << "[EC] Final sub-space details: " << (n - k) << " eliminated, " << k << " active in subset.\n";
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n + 4 * (n - k));
    for (int i = 0; i < n; i++) {
        for (auto& [full_j, coeff] : expr[i]) {
            int s = full2subset[full_j];
            if (s >= 0) {
                triplets.emplace_back(i, s, coeff);
            }
        }
    }

    P.resize(n, k);
    P.setFromTriplets(triplets.begin(), triplets.end());

    // =========================================================================
    // DIAGNOSTIC CHECKS: Validate Nullspace and Integrality Constraints
    // =========================================================================
    if (m > 0) {
        // --- 1. Nullspace Validity Check ---
        Eigen::SparseMatrix<double> CP = C * P;
        double cp_err = 0.0;
        for (int r = 0; r < CP.outerSize(); ++r) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(CP, r); it; ++it) {
                cp_err += it.value() * it.value();
            }
        }
        cp_err = std::sqrt(cp_err);

        Eigen::VectorXd test_sub = Eigen::VectorXd::Random(k);
        Eigen::VectorXd test_full = P * test_sub;
        Eigen::VectorXd test_resid = C * test_full;
        double probe_err = test_resid.norm();

        // --- 2. Integrality Constraint Check ---
        // Verify that any relationship between integer parameters consists purely of strict integer coefficients.
        double max_fractional_error = 0.0;
        int non_integer_count = 0;

        for (int i = 0; i < n; i++) {
            // Only examine variables that were flagged as integer types originally
            if (hasPreference && preferSurvive(i) == 1) {
                for (auto& [active_col, coeff] : expr[i]) {
                    // Check deviation from the nearest whole number
                    double diff = std::abs(coeff - std::round(coeff));
                    if (diff > max_fractional_error) {
                        max_fractional_error = diff;
                    }
                    if (diff > 1e-4) {
                        non_integer_count++;
                    }
                }
            }
        }

        if (verbose) {
            std::cout << "[EC-DIAGNOSTIC] Matrix product norm ||C * P||: " << cp_err << "\n";
            std::cout << "[EC-DIAGNOSTIC] Random subspace probe validation ||C * P * x_sub||: " << probe_err << "\n";
            std::cout << "[EC-DIAGNOSTIC] Maximum fractional entry error inside integer relationships: " << max_fractional_error << "\n";
            
            if (non_integer_count > 0) {
                std::cout << "[EC-DIAGNOSTIC] !!! WARNING: " << non_integer_count << " fractional mappings found affecting integer parameters !!!\n";
            } else {
                std::cout << "[EC-DIAGNOSTIC] Integrality Check Passed: P maps integer dimensions flawlessly.\n";
            }

            if (probe_err > 1e-5 || cp_err > 1e-5) {
                std::cout << "[EC-DIAGNOSTIC] !!! CRITICAL WARNING: P failed to satisfy C structurally !!!\n";
            }
        }
    }
}
