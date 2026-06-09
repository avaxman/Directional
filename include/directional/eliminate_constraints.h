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
    const Eigen::VectorXi&             preferSurvive = Eigen::VectorXi(),
    bool                               verbose = false)
{
    const int n = C.cols();
    const int m = C.rows();
    const double eps = 1e-14;
    const bool hasPreference = (preferSurvive.size() == n);

    if (verbose)
    {
        std::cout << "[EC] eliminate_constraints: m=" << m << " n=" << n
                  << " hasPreference=" << hasPreference;
        if (hasPreference)
            std::cout << " preferSurvive.sum()=" << preferSurvive.sum();
        std::cout << "\n";
        if (hasPreference)
        {
            int nPref = 0;
            for (int i = 0; i < n; i++)
                if (preferSurvive(i)) nPref++;
            std::cout << "[EC] " << nPref << " preferred-survive variables\n";
        }
    }

    std::vector<std::map<int, double>> expr(n);
    for (int i = 0; i < n; i++)
        expr[i][i] = 1.0;

    std::vector<bool> eliminated(n, false);

    std::vector<std::map<int, double>> rows(m);
    {
        Eigen::SparseMatrix<double, Eigen::RowMajor> Crow(C);
        for (int r = 0; r < m; r++)
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Crow, r); it; ++it)
                if (std::abs(it.value()) > eps)
                    rows[r][it.col()] = it.value();
    }

    // ================================================================
    // Sort rows:
    // 0 = mixed (has both continuous and integer vars) - process first
    // 1 = all-continuous
    // 2 = all-integer - process last
    // 3 = empty
    // Within mixed: fewer integer vars first (rank = nInt, so 1 < 2 < 3...)
    // ================================================================
    if (hasPreference)
    {
        auto rowRank = [&](const std::map<int, double>& row) -> int
        {
            int nCont = 0, nInt = 0;
            for (auto& [col, val] : row)
            {
                if (preferSurvive(col)) nInt++;
                else                    nCont++;
            }
            if (nCont == 0 && nInt == 0) return 1000; // empty
            if (nCont == 0)              return 500 + nInt; // all-integer: last
            if (nInt == 0)               return 200;        // all-continuous
            return nInt; // mixed: fewer integer vars first
        };

        std::vector<int> order(m);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(),
            [&](int a, int b){ return rowRank(rows[a]) < rowRank(rows[b]); });

        std::vector<std::map<int, double>> sortedRows(m);
        for (int i = 0; i < m; i++)
            sortedRows[i] = std::move(rows[order[i]]);
        rows = std::move(sortedRows);
    }

    // Substitute pivot_col out of all other expr[] entries.
    // CRITICAL: skip free preferred (integer) variables - their expr is
    // identity {i->1} and must stay that way. Fractions from continuous
    // eliminations must never contaminate free integer variable expressions.
    auto propagate_into_exprs = [&](int piv, const std::map<int, double>& piv_expr)
    {
        for (int i = 0; i < n; i++)
        {
            if (i == piv) continue;
            if (hasPreference && preferSurvive(i) && !eliminated[i]) continue;

            auto it = expr[i].find(piv);
            if (it == expr[i].end()) continue;
            double coeff = it->second;
            expr[i].erase(it);
            for (auto& [k, v] : piv_expr)
            {
                expr[i][k] += coeff * v;
                if (std::abs(expr[i][k]) < eps)
                    expr[i].erase(k);
            }
        }
    };

    int nPivotsContinuous = 0;
    int nPivotsPreferred  = 0;
    int nRowsSkipped      = 0;

    // Track the boundary where pure-integer rows start (for diagnostics)
    bool inIntegerPhase = false;

    for (int r = 0; r < m; r++)
    {
        int    pivot_col      = -1;
        double pivot_val      = 0.0;
        int    pivot_col_fall = -1;
        double pivot_val_fall = 0.0;

        int nFreeInRow     = 0;
        int nFreeContInRow = 0;
        int nFreeIntInRow  = 0;

        for (auto& [col, val] : rows[r])
        {
            if (eliminated[col]) continue;
            nFreeInRow++;
            bool preferred = hasPreference && (preferSurvive(col) != 0);
            if (preferred) nFreeIntInRow++;
            else           nFreeContInRow++;

            if (preferred)
            {
                if (std::abs(val) > std::abs(pivot_val_fall))
                {
                    pivot_col_fall = col;
                    pivot_val_fall = val;
                }
            }
            else
            {
                if (std::abs(val) > std::abs(pivot_val))
                {
                    pivot_col = col;
                    pivot_val = val;
                }
            }
        }

        // Detect transition to pure-integer phase
        if (hasPreference && !inIntegerPhase && nFreeContInRow == 0 && nFreeIntInRow > 0)
        {
            inIntegerPhase = true;
            if (verbose)
            {
                std::cout << "[EC] === entering pure-integer phase at row " << r << " ===\n";
                // Check all integer rows for non-unit coefficients
                int nBadRows = 0;
                for (int r2 = r; r2 < m; r2++)
                {
                    bool rowHasCont = false;
                    bool rowHasFrac = false;
                    for (auto& [col, val] : rows[r2])
                    {
                        if (eliminated[col]) continue;
                        if (!preferSurvive(col)) { rowHasCont = true; break; }
                        if (std::abs(std::abs(val) - 1.0) > 1e-6) rowHasFrac = true;
                    }
                    if (rowHasCont) continue; // skip non-integer rows
                    if (rowHasFrac)
                    {
                        nBadRows++;
                        if (nBadRows <= 5)
                        {
                            std::cout << "[EC]   row " << r2 << " has NON-UNIT integer coefficients: ";
                            for (auto& [col, val] : rows[r2])
                                if (!eliminated[col])
                                    std::cout << val << "*x[" << col << "] ";
                            std::cout << "\n";
                        }
                    }
                }
                if (nBadRows == 0)
                    std::cout << "[EC]   all integer rows have unit coefficients - OK\n";
                else
                    std::cout << "[EC]   " << nBadRows << " integer rows have non-unit coefficients - CONTAMINATED\n";
            }
        }

        if (verbose && r < 110)
        {
            std::cout << "[EC] row " << r
                      << ": freeTotal=" << nFreeInRow
                      << " freeCont="   << nFreeContInRow
                      << " freeInt="    << nFreeIntInRow;
            if (pivot_col >= 0)
                std::cout << " -> pivot=x[" << pivot_col << "] (continuous)"
                          << " val=" << pivot_val;
            else if (pivot_col_fall >= 0)
                std::cout << " -> pivot=x[" << pivot_col_fall << "] (FORCED integer fallback)"
                          << " val=" << pivot_val_fall;
            else
                std::cout << " -> SKIPPED";
            std::cout << "\n";
        }

        if (pivot_col == -1)
        {
            pivot_col = pivot_col_fall;
            pivot_val = pivot_val_fall;
            if (pivot_col >= 0) nPivotsPreferred++;
        }
        else
        {
            nPivotsContinuous++;
        }

        if (pivot_col == -1)
        {
            nRowsSkipped++;
            continue;
        }

        eliminated[pivot_col] = true;

        // Build pivot expression
        std::map<int, double> pivot_expr;
        for (auto& [col, val] : rows[r])
        {
            if (col == pivot_col || std::abs(val) < eps) continue;
            double coeff = -val / pivot_val;
            for (auto& [k, v] : expr[col])
            {
                pivot_expr[k] += coeff * v;
                if (std::abs(pivot_expr[k]) < eps)
                    pivot_expr.erase(k);
            }
        }

        // Diagnostic: warn if integer pivot produces non-unit coefficients
        if (verbose && hasPreference && preferSurvive(pivot_col))
        {
            for (auto& [k, v] : pivot_expr)
                if (std::abs(std::abs(v) - 1.0) > 1e-6 && std::abs(v) > eps)
                    std::cout << "[EC]   WARNING: integer pivot x[" << pivot_col
                              << "] produced non-unit coeff " << v
                              << " for x[" << k << "]\n";
        }

        expr[pivot_col] = pivot_expr;
        propagate_into_exprs(pivot_col, pivot_expr);

        // Propagate into remaining rows
        for (int r2 = r + 1; r2 < m; r2++)
        {
            auto it = rows[r2].find(pivot_col);
            if (it == rows[r2].end()) continue;
            double factor = it->second;
            rows[r2].erase(it);
            double scale = -factor / pivot_val;
            for (auto& [col, val] : rows[r])
            {
                if (col == pivot_col) continue;
                double& entry = rows[r2][col];
                entry += scale * val;
                if (std::abs(entry) < eps)
                    rows[r2].erase(col);
            }
        }
    }

    if (verbose)
    {
        std::cout << "[EC] pivots: " << nPivotsContinuous << " continuous, "
                  << nPivotsPreferred << " forced-integer, "
                  << nRowsSkipped << " rows skipped\n";

        int nPrefElim = 0;
        for (int i = 0; i < n; i++)
        {
            if (hasPreference && preferSurvive(i) && eliminated[i])
            {
                nPrefElim++;
                if (nPrefElim <= 10)
                    std::cout << "[EC] preferred x[" << i << "] eliminated\n";
            }
        }
        if (nPrefElim > 10)
            std::cout << "[EC] ... and " << (nPrefElim-10) << " more\n";
        std::cout << "[EC] total preferred eliminated: " << nPrefElim << "\n";

        // Check final expressions of eliminated integer vars
        if (hasPreference)
        {
            int nBadExpr = 0;
            for (int i = 0; i < n; i++)
            {
                if (!preferSurvive(i) || !eliminated[i]) continue;
                bool bad = false;
                for (auto& [j, c] : expr[i])
                    if (std::abs(std::abs(c) - 1.0) > 1e-6) bad = true;
                if (bad)
                {
                    nBadExpr++;
                    if (nBadExpr <= 5)
                    {
                        std::cout << "[EC] BAD expr[x[" << i << "]] =";
                        for (auto& [j, c] : expr[i])
                            std::cout << " " << c << "*x[" << j << "]";
                        std::cout << "\n";
                    }
                }
            }
            if (nBadExpr == 0)
                std::cout << "[EC] checking expr of eliminated integer vars: all OK (+-1)\n";
            else
                std::cout << "[EC] checking expr: " << nBadExpr << " eliminated integer vars have non-unit coefficients\n";
        }
    }

    // Build index maps
    full2subset.setConstant(n, -1);
    std::vector<int> s2f;
    s2f.reserve(n);
    for (int i = 0; i < n; i++)
    {
        if (!eliminated[i])
        {
            full2subset[i] = (int)s2f.size();
            s2f.push_back(i);
        }
    }
    const int k = (int)s2f.size();
    subset2full.resize(k);
    for (int s = 0; s < k; s++)
        subset2full[s] = s2f[s];

    if (verbose)
        std::cout << "[EC] result: " << (n-k) << " eliminated, " << k << " in subset\n";

    // Assemble P
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(n + 4 * (n - k));
    for (int i = 0; i < n; i++)
    {
        for (auto& [full_j, coeff] : expr[i])
        {
            assert(!eliminated[full_j]);
            int s = full2subset[full_j];
            assert(s >= 0);
            triplets.emplace_back(i, s, coeff);
        }
    }

    P.resize(n, k);
    P.setFromTriplets(triplets.begin(), triplets.end());
    if (verbose){
        std::cout << "P (sparse) =\n";

            for (int k = 0; k < P.outerSize(); ++k)
            {¯
                for (Eigen::SparseMatrix<double>::InnerIterator it(P, k); it; ++it)
                {
                    std::cout << "  (" << it.row() << ", " << it.col()
                              << ") = " << it.value() << "\n";
                }
            }
    }
}
