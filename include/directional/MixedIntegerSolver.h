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
    
        //Initial reduction of the constrain matrix to not contain redundant constraints
        /*SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
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
        }*/
        
        
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
            SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
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
            SparseMatrix<double> AtMA = APart.transpose() * M * APart;
            CPart = C * var2AllMat;   //keep reducing rows from the already-reduced constraints matrix (the number of columns stays the same TODO THAT
            
            //reducing rank on Cpart
            int CPartRank=0;
            VectorXi PIndices(0);
            if (CPart.rows()!=0){
                qrsolver.compute(CPart.transpose());
                CPartRank = qrsolver.rank();
                
                //creating sliced permutation matrix
                PIndices = qrsolver.colsPermutation().indices();
                
                vector<Triplet<double> > CPartTriplets;
                
                for(int k = 0; k < CPart.outerSize(); ++k)
                {
                    for (SparseMatrix<double>::InnerIterator it(CPart, k); it; ++it)
                    {
                        for (int j = 0; j < CPartRank; j++)
                            if (it.row() == PIndices(j))
                                CPartTriplets.emplace_back(j, it.col(), it.value());
                    }
                }
                
                CPart.resize(CPartRank, CPart.cols());
                CPart.setFromTriplets(CPartTriplets.begin(), CPartTriplets.end());
            }
            SparseMatrix<double> E(AtMA.rows()+ CPart.rows(), AtMA.rows() + CPart.rows());
            
            vector<Triplet<double>> ETriplets;
            for(int k = 0; k < AtMA.outerSize(); ++k)
            {
                for (SparseMatrix<double>::InnerIterator it(AtMA, k); it; ++it)
                    ETriplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
            }
            
            for(int k = 0; k < CPart.outerSize(); ++k)
            {
                for(SparseMatrix<double>::InnerIterator it(CPart, k); it; ++it)
                {
                    ETriplets.emplace_back(it.row() + AtMA.rows(), it.col(), it.value());
                    ETriplets.emplace_back(it.col(), it.row() + AtMA.rows(), it.value());
                }
            }
            
            E.setFromTriplets(ETriplets.begin(), ETriplets.end());
            
            //Building Right-hand side with fixed values
            VectorXd rhs = VectorXd::Zero(AtMA.rows() + CPart.rows());
            rhs.segment(0, AtMA.rows())= APart.transpose() * M * (b + torhs);
            VectorXd dfull = -C * fixedValues;
            VectorXd dPart(CPartRank);
            for(int k = 0; k < CPartRank; k++)
                dPart(k)=dfull(PIndices(k));
            rhs.segment(AtMA.rows(), CPart.rows()) = dPart;
            
            SparseLU<SparseMatrix<double> > lusolver;
            lusolver.compute(E);
            if(lusolver.info() != Success){
                if (verbose)
                    cout<<"LU decomposition failed!"<<endl;
                return false;
            }
            xPart = lusolver.solve(rhs);
            
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
                    if (verbose)
                        std::cout<<"Index "<<i<<" already pre-fixed; ignoring"<<std::endl;
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
                        if (verbose)
                            std::cout<<"adding tolerance-rounded "<<i<<" with current value "<<func<<"and integer difference "<<currIntDiff<<" to fixed variables"<<std::endl;
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
        
        //for(int intIter = 0; intIter < fixedMask.sum(); intIter++)
        //{
            
            /*xprev.resize(x.rows() - 1);
            varCounter = 0;
            for(int i = 0; i < numVars; i++)
                if (!alreadyFixed(i))
                    xprev(varCounter++) = fullx(i);
            
            xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());*/
        //}
        
        return true;
    }
    
    
};

}


#endif
