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



//This is a "home-made" mixed integer solver for quadratic optimization with linear constraints. This might be slow-ish.
//This problem solved is as following:
// x = argmin |Ax - b|^2 with metric matrix M, so that Cx = 0, and where a subset x_int \in Z, and x_fixed = fixedValues;
//fixed values are known in advance, and get removed from the system immediately. Values that have been rounded off join the fixed mask after having been rounded.
//the toRound values are to be fixed by rounding. After having been rounded they are removed from the list. Even if they were previously rounded, it is not checked here.


namespace directional{

enum class RoundingTypeEnum{ITERATIVE_ROUNDING, BRANCH_AND_BOUND};


bool solveConstrainedLeastSquares(const Eigen::SparseMatrix<double>& A,
                                  const Eigen::SparseMatrix<double>& M,
                                  const Eigen::VectorXd& b,
                                  const Eigen::SparseMatrix<double>& C,
                                  const Eigen::VectorXd& d,
                                  Eigen::VectorXd& x,
                                  double regEps = 1e-10)
{
    using namespace Eigen;

    // H = A^T M A
    SparseMatrix<double> H = A.transpose() * M * A;

    // regularization
    H.coeffRef(0, 0) += regEps; // cheaper than building identity matrix

    SimplicialLDLT<SparseMatrix<double>> ldlt;
    ldlt.analyzePattern(H);
    ldlt.factorize(H);

    if (ldlt.info() != Success)
        return false;

    const int n = H.rows();
    const int m = C.rows();

    // r = A^T M b
    VectorXd r = A.transpose() * M * b;

    // -----------------------------
    // FAST PATH: build Y in one shot
    // -----------------------------
    MatrixXd Ct = C.transpose().eval();   // force contiguous dense RHS
    MatrixXd Y  = ldlt.solve(Ct);         // single block solve

    // S = C H^{-1} C^T
    MatrixXd S = C * Y;

    // rhs = C H^{-1} r - d
    VectorXd Hr = ldlt.solve(r);
    VectorXd rhs_lambda = C * Hr - d;

    // Solve Schur system
    VectorXd lambda;

    // Prefer LDLT if possible (S is symmetric PSD)
    LDLT<MatrixXd> sldlt(S);
    if (sldlt.info() == Success)
    {
        lambda = sldlt.solve(rhs_lambda);
    }
    else
    {
        // fallback
        lambda = S.colPivHouseholderQr().solve(rhs_lambda);
    }

    // primal recovery
    x = ldlt.solve(r - C.transpose() * lambda);

    return true;
}

class MixedIntegerSolver{
public:
    RoundingTypeEnum roundingType;
    
    Eigen::SparseMatrix<double> A, C, M;
    Eigen::VectorXd b, fixedValues;
    Eigen::VectorXi fixedMask;
    Eigen::VectorXi toRoundMask;
    
    int numVars;
    bool verbose;
    double integerTolerance;
    
    Eigen::VectorXd x;
    
    MixedIntegerSolver():verbose(false), numVars(0), integerTolerance(10e-12){}
    ~MixedIntegerSolver(){}
    
    
    
    bool solve(){
        
        using namespace Eigen;
        using namespace std;
        
        //Currently online doing iterative rounding
        Eigen::SparseMatrix<double> CPart = C;
        VectorXd xPart;
        
        
        if (verbose)
            std::cout<<"the initial #fixing variables: "<<fixedMask.sum()<<std::endl;
        
        if (verbose)
            std::cout<<"Need to round "<<toRoundMask.sum()<<" integers"<<std::endl;
        
        //This is done at least once, and until the toRound set is emptied
        while(true)
        {
            if (verbose)
                std::cout<<"Solving system with "<<numVars - fixedMask.sum()<<" variables"<<std::endl;
            SparseMatrix<double> var2AllMat;
           //the non-fixed variables to all variables
            var2AllMat.resize(numVars, numVars - fixedMask.sum());
            int varCounter = 0;
            vector<Triplet<double> > var2AllTriplets;
            for(int i = 0; i < numVars; i++)
                if (!fixedMask(i))
                    var2AllTriplets.emplace_back(i, varCounter++, 1.0);
            
            var2AllMat.setFromTriplets(var2AllTriplets.begin(), var2AllTriplets.end());
            
            SparseMatrix<double> APart = A * var2AllMat;
            VectorXd torhs = -A * fixedValues;
            CPart = C * var2AllMat;   //keep reducing rows from the already-reduced constraints matrix (the number of columns stays the same TODO
            VectorXd dFull = -C * fixedValues;
            if (!solveConstrainedLeastSquares(APart, M, b + torhs, CPart, dFull, xPart)){
                if (verbose)
                    std::cout<<"solveConstrainedLeastSquares() has failed!"<<std::endl;
                return false;
            }
            
            x = var2AllMat * xPart.head(numVars - fixedMask.sum()) + fixedValues;
            
            if (toRoundMask.sum()==0)
                break;  //nothing more to do
            
            //Otherwise, find the next roundable integer
            double minIntDiff = std::numeric_limits<double>::max();
            int minIntDiffIndex = -1;
            bool changedFixed = false;
            for (int i = 0; i < numVars; i++)
            {
                if ((toRoundMask(i))&&(fixedMask(i))){  //this variable has already been rounded before (assuming it's integer!)
                    toRoundMask(i) = 0;
                    //if (verbose)
                    //    std::cout<<"Index "<<i<<" already pre-fixed; ignoring"<<std::endl;
                    continue;
                }
                
                if ((toRoundMask(i))&&(!fixedMask(i)))  //a variable to round that has not been fixed already
                {
                    changedFixed = true;  //there is an event that changed fixed variable, and need resolving
                    double currIntDiff =0;
                    double func = x(i);
                    currIntDiff += std::fabs(func - std::round(func));
                    if (currIntDiff<integerTolerance){  //already been rounded to tolerance
                        fixedMask(i) = 1;
                        toRoundMask(i) = 0;
                        //if (verbose)
                        //    std::cout<<"adding tolerance-rounded "<<i<<" with current value "<<func<<" and integer difference "<<currIntDiff<<" to fixed variables"<<std::endl;
                        fixedValues(i) = std::round(func);
                    }else if (currIntDiff < minIntDiff)
                    {
                        minIntDiff = currIntDiff;
                        minIntDiffIndex = i;
                    }
                }
            }
            
            if (minIntDiffIndex != -1)
            {
                fixedMask(minIntDiffIndex) = 1;
                toRoundMask(minIntDiffIndex) = 0;
                double func = x(minIntDiffIndex) ;
                double funcInteger=std::round(func);
                if (verbose)
                    std::cout<<"rounding index "<<minIntDiffIndex<<" from "<<func<<" to "<<funcInteger<<std::endl;
                fixedValues(minIntDiffIndex) = funcInteger;
            }
            
            if (!changedFixed)  //no resolving needed
                break;
        }
        
        
        return true;
    }
    
    
};

}


#endif
