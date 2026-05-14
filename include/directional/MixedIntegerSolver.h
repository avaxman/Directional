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

namespace directional{

enum class RoundingTypeEnum{ITERATIVE_ROUNDING, BRANCH_AND_BOUND};

class MixedIntegerSolver{
public:
    RoundingTypeEnum roundingType;
    
    Eigen::SparseMatrix<double> A, C, M;
    Eigen::VectorXd b, fixedValues;
    Eigen::VectorXi fixedMask, alreadyFixedMask;  //alreadyfixed are values that are already done, fixedMask is all values (including the already-ones) that need to be fixed;
    int numVars;
    bool verbose;
    
    Eigen::VectorXd x;
    
    MixedIntegerSolver():verbose(false), numVars(0){}
    ~MixedIntegerSolver(){}
    
    
    bool solve(){
        
        using namespace Eigen;
        using namespace std;
        
        
        //Currently online doing iterative rounding
        Eigen::SparseMatrix<double> CFull = C;
        VectorXd xPart;//, xprev;
    
        //reducing constraintMat
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
        if (CFull.rows()!=0){
            qrsolver.compute(CFull.transpose());
            int CRank = qrsolver.rank();
            
            //creating sliced permutation matrix
            VectorXi PIndices = qrsolver.colsPermutation().indices();
            
            vector<Triplet<double> > CTriplets;
            for(int k = 0; k < CFull.outerSize(); ++k)
            {
                for(SparseMatrix<double>::InnerIterator it(CFull, k); it; ++it)
                {
                    for(int j = 0; j < CRank; j++)
                        if(it.row() == PIndices(j))
                            CTriplets.emplace_back(j, it.col(), it.value());
                }
            }
            
            CFull.resize(CRank, CFull.cols());
            CFull.setFromTriplets(CTriplets.begin(), CTriplets.end());
        }
        SparseMatrix<double> var2AllMat;
        for(int intIter = 0; intIter < fixedMask.sum(); intIter++)
        {
            //the non-fixed variables to all variables
            var2AllMat.resize(numVars, numVars - alreadyFixedMask.sum());
            int varCounter = 0;
            vector<Triplet<double> > var2AllTriplets;
            for(int i = 0; i < numVars; i++)
            {
                if (!alreadyFixedMask(i)){
                    //for (int j=0;j<intData.d;j++)
                    var2AllTriplets.emplace_back(i, varCounter++, 1.0);
                }
                
            }
            var2AllMat.setFromTriplets(var2AllTriplets.begin(), var2AllTriplets.end());
            
            SparseMatrix<double> Apart = A * var2AllMat;
            VectorXd torhs = -A * fixedValues;
            SparseMatrix<double> AtMA = Apart.transpose() * M * Apart;
            SparseMatrix<double> Cpart = CFull * var2AllMat;
            
            //reducing rank on Cpart
            int CpartRank=0;
            VectorXi PIndices(0);
            if (Cpart.rows()!=0){
                qrsolver.compute(Cpart.transpose());
                CpartRank = qrsolver.rank();
                
                //creating sliced permutation matrix
                PIndices = qrsolver.colsPermutation().indices();
                
                vector<Triplet<double> > CPartTriplets;
                
                for(int k = 0; k < Cpart.outerSize(); ++k)
                {
                    for (SparseMatrix<double>::InnerIterator it(Cpart, k); it; ++it)
                    {
                        for (int j = 0; j < CpartRank; j++)
                            if (it.row() == PIndices(j))
                                CPartTriplets.emplace_back(j, it.col(), it.value());
                    }
                }
                
                Cpart.resize(CpartRank, Cpart.cols());
                Cpart.setFromTriplets(CPartTriplets.begin(), CPartTriplets.end());
            }
            SparseMatrix<double> E(AtMA.rows()+ Cpart.rows(), AtMA.rows() + Cpart.rows());
            
            vector<Triplet<double>> ETriplets;
            for(int k = 0; k < AtMA.outerSize(); ++k)
            {
                for (SparseMatrix<double>::InnerIterator it(AtMA, k); it; ++it)
                    ETriplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
            }
            
            for(int k = 0; k < Cpart.outerSize(); ++k)
            {
                for(SparseMatrix<double>::InnerIterator it(Cpart, k); it; ++it)
                {
                    ETriplets.emplace_back(it.row() + AtMA.rows(), it.col(), it.value());
                    ETriplets.emplace_back(it.col(), it.row() + AtMA.rows(), it.value());
                }
            }
            
            E.setFromTriplets(ETriplets.begin(), ETriplets.end());
            
            //Building Right-hand side with fixed values
            VectorXd rhs = VectorXd::Zero(AtMA.rows() + Cpart.rows());
            rhs.segment(0, AtMA.rows())= Apart.transpose() * M * (b + torhs);
            VectorXd dfull = -CFull * fixedValues;
            VectorXd dPart(CpartRank);
            for(int k = 0; k < CpartRank; k++)
                dPart(k)=dfull(PIndices(k));
            rhs.segment(AtMA.rows(), Cpart.rows()) = dPart;
            
            SparseLU<SparseMatrix<double> > lusolver;
            lusolver.compute(E);
            if(lusolver.info() != Success){
                if (verbose)
                    cout<<"LU decomposition failed!"<<endl;
                return false;
            }
            xPart = lusolver.solve(rhs);
            
            x = var2AllMat * xPart.head(numVars - alreadyFixedMask.sum()) + fixedValues;
            
            if((alreadyFixedMask - fixedMask).sum() == 0)
                break;
            
            double minIntDiff = std::numeric_limits<double>::max();
            int minIntDiffIndex = -1;
            for (int i = 0; i < numVars; i++)
            {
                if ((fixedMask(i)) && (!alreadyFixedMask(i)))
                {
                    double currIntDiff =0;
                    double func = x(i); //fullx.segment(intData.d*i,intData.d);
                    //for (int j=0;j<intData.d;j++)
                    currIntDiff += std::fabs(func - std::round(func));
                    if (currIntDiff < minIntDiff)
                    {
                        minIntDiff = currIntDiff;
                        minIntDiffIndex = i;
                    }
                }
            }
            
            if (minIntDiffIndex != -1)
            {
                alreadyFixedMask(minIntDiffIndex) = 1;
                double func = x(minIntDiffIndex) ;
                double funcInteger=std::round(func);
                if (verbose)
                    std::cout<<"rounding index "<<minIntDiffIndex<<" from "<<func<<" to "<<funcInteger<<std::endl;
                fixedValues(minIntDiffIndex) = funcInteger;
            }
            
            /*xprev.resize(x.rows() - 1);
            varCounter = 0;
            for(int i = 0; i < numVars; i++)
                if (!alreadyFixed(i))
                    xprev(varCounter++) = fullx(i);
            
            xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());*/
        }
        
        return true;
    }
    
    
};

}


#endif
