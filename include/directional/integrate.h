// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_INTEGRATE_H
#define DIRECTIONAL_INTEGRATE_H

#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <directional/tree.h>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/principal_matching.h>
#include <directional/setup_integration.h>
//#include <directional/branched_gradient.h>
#include <directional/gradient_matrices.h>
#include <directional/mass_matrices.h>
#include <directional/MixedIntegerSolver.h>


namespace directional
{

// Integrates an N-directional fields into an N-function by solving the seamless Poisson equation. Respects *valid* linear reductions where the field is reducible to an n-field for n<=M, and consequently the function is reducible to an n-function.
// This function only works with face-based fields on triangle meshes.
// Input:
//  field:              The face-based field to be integrated, on the original mesh
//  intData:            Integration data, which must be obtained from directional::setup_integration(). This is altered by the function.
//  meshCut:            Cut mesh (obtained from setup_integration())
// Output:
//  NFunction:          #cV x N parameterization functions per cut vertex (full version with all symmetries unpacked)
//  NCornerFunctions   (3*N) x #F parameterization functions per corner of whole mesh
inline bool integrate(const directional::CartesianField& field,
                      IntegrationData& intData,
                      const directional::TriMesh& meshCut,
                      Eigen::MatrixXd& NFunction,
                      Eigen::MatrixXd& NCornerFunctions)


{
    using namespace Eigen;
    using namespace std;
    
    assert(field.tb->discTangType()==discTangTypeEnum::FACE_SPACES && "Integrate() only works with face-based fields");
    const directional::TriMesh& meshWhole = *((PCFaceTangentBundle*)field.tb)->mesh;
    
    VectorXd edgeWeights = VectorXd::Constant(meshWhole.FE.maxCoeff() + 1, 1.0);
    double paramLength = (meshWhole.V.colwise().maxCoeff()-meshWhole.V.colwise().minCoeff()).norm()*intData.lengthRatio;
    
    //Scaling field uniformly to have average length intData.lengthRatio * bounding_box_diagonal (Default is 0.02)
    VectorXd rawField = field.flatten();
    double avgGradNorm=0;
    for (int i=0;i<meshCut.F.rows();i++)
        for (int j=0;j<intData.N;j++)
            avgGradNorm+=rawField.segment(3*intData.N*i+3*j,3).norm();
    
    avgGradNorm/=(double)(intData.N*meshCut.F.rows());
    
    std::cout<<"paramLength: "<<paramLength<<std::endl;
    std::cout<<"avgGradNorm: "<<avgGradNorm<<std::endl;
    rawField.array()*=1.0/(paramLength*avgGradNorm);
    //paramLength/=avgGradNorm;
    
    int numVars = intData.linRedMat.cols();
    
    //constructing face differentials
    
    SparseMatrix<double> G = directional::conf_gradient_matrix_2D<double>(meshCut, false, intData.N);
    SparseMatrix<double> Mx = directional::face_mass_matrix_2D<double>(meshCut, false, 3*intData.N);
    SparseMatrix<double> GT = G.transpose();
    

    //The variables that should be fixed in the end
    VectorXi fixedMask(numVars);
    fixedMask.setZero();
    
    for (int i=0;i<intData.fixedIndices.size();i++)
        fixedMask(intData.fixedIndices(i)) = 1;
    
    bool roundedSingularities = false;  //if all singularities have been rounded (only relevant to intData.roundSeams=false)
    if(intData.integralSeamless) {
        if (intData.roundSeams) {
            for (int i = 0; i < intData.integerVars.size(); i++)
                for (int j = 0; j < intData.n; j++)
                    fixedMask(intData.n * intData.integerVars(i) + j) = 1;
        }else {
            for (int i = 0; i < intData.singularIndices.size(); i++)
                fixedMask(intData.singularIndices(i)) = 1;
        }
    }
    
    //the variables that were already fixed to begin with
    VectorXi alreadyFixed(numVars);
    alreadyFixed.setZero();
    
    
    for (int i=0;i<intData.fixedIndices.size();i++)
        alreadyFixed(intData.fixedIndices(i)) = 1;
    
    //the values for the fixed variables (size is as all variables)
    VectorXd fixedValues(numVars);
    fixedValues.setZero();  //for everything but the originally fixed values
    for (int i=0;i<intData.fixedValues.size();i++)
        fixedValues(intData.fixedIndices(i))=intData.fixedValues(i);
    
    directional::MixedIntegerSolver mis;
    
    mis.A = G * intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
    mis.M = Mx;
    mis.fixedMask = fixedMask;
    mis.fixedValues = fixedValues;
    mis.alreadyFixedMask = alreadyFixed;
    mis.C = intData.constraintMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
    mis.b = rawField;
    mis.numVars = numVars;
    mis.verbose = intData.verbose;
    bool success = mis.solve();
    if (!success)
        return success;
    VectorXd xFull= mis.x;
    
    //solve again for extra integers for topological unrounded seams
    if ((!intData.roundSeams)&&(!roundedSingularities)&&(intData.integralSeamless)) {
        for (int i = 0; i < intData.integerVars.size(); i++)
            for (int j = 0; j < intData.n; j++)
                mis.fixedMask(intData.n * intData.integerVars(i) + j) = 1;
        
        roundedSingularities = true;
    }
    if (intData.verbose)
        std::cout<<"Rounding extra unrounded topological seams if needed"<<std::endl;
    
    success = mis.solve();
    if (!success)
        return success;
    xFull= mis.x;
    
    //need to integrate this code:
    /****
     //in case all singularities are rounded in the rounding-singularities mode, but there are left unrounded seams (like topological handles).
     if ((alreadyFixed.sum()==fixedMask.sum())&&(!intData.roundSeams)&(!roundedSingularities)&&(intData.integralSeamless)) {
         for (int i = 0; i < intData.integerVars.size(); i++)
             for (int j = 0; j < intData.n; j++)
                 fixedMask(intData.n * intData.integerVars(i) + j) = 1;
         roundedSingularities = true;
     }
     ***//////
    

    /*SparseMatrix<double> Efull = G * intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
    VectorXd x, xprev;
    
    // until then all the N depedencies should be resolved?
    
    //reducing constraintMat
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
    SparseMatrix<double> Cfull = intData.constraintMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
    if (Cfull.rows()!=0){
        qrsolver.compute(Cfull.transpose());
        int CRank = qrsolver.rank();
        
        //creating sliced permutation matrix
        VectorXi PIndices = qrsolver.colsPermutation().indices();
        
        vector<Triplet<double> > CTriplets;
        for(int k = 0; k < Cfull.outerSize(); ++k)
        {
            for(SparseMatrix<double>::InnerIterator it(Cfull, k); it; ++it)
            {
                for(int j = 0; j < CRank; j++)
                    if(it.row() == PIndices(j))
                        CTriplets.emplace_back(j, it.col(), it.value());
            }
        }
        
        Cfull.resize(CRank, Cfull.cols());
        Cfull.setFromTriplets(CTriplets.begin(), CTriplets.end());
    }
    SparseMatrix<double> var2AllMat;
    VectorXd fullx(numVars); fullx.setZero();
    for(int intIter = 0; intIter < fixedMask.sum(); intIter++)
    {
        //the non-fixed variables to all variables
        var2AllMat.resize(numVars, numVars - alreadyFixed.sum());
        int varCounter = 0;
        vector<Triplet<double> > var2AllTriplets;
        for(int i = 0; i < numVars; i++)
        {
            if (!alreadyFixed(i)){
                //for (int j=0;j<intData.d;j++)
                var2AllTriplets.emplace_back(i, varCounter++, 1.0);
            }
            
        }
        var2AllMat.setFromTriplets(var2AllTriplets.begin(), var2AllTriplets.end());
        
        SparseMatrix<double> Epart = Efull * var2AllMat;
        VectorXd torhs = -Efull * fixedValues;
        SparseMatrix<double> EtE = Epart.transpose() * Mx * Epart;
        SparseMatrix<double> Cpart = Cfull * var2AllMat;
        
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
        SparseMatrix<double> A(EtE.rows()+ Cpart.rows(), EtE.rows() + Cpart.rows());
        
        vector<Triplet<double>> ATriplets;
        for(int k = 0; k < EtE.outerSize(); ++k)
        {
            for (SparseMatrix<double>::InnerIterator it(EtE, k); it; ++it)
                ATriplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
        }
        
        for(int k = 0; k < Cpart.outerSize(); ++k)
        {
            for(SparseMatrix<double>::InnerIterator it(Cpart, k); it; ++it)
            {
                ATriplets.emplace_back(it.row() + EtE.rows(), it.col(), it.value());
                ATriplets.emplace_back(it.col(), it.row() + EtE.rows(), it.value());
            }
        }
        
        A.setFromTriplets(ATriplets.begin(), ATriplets.end());
        
        //Right-hand side with fixed values
        VectorXd b = VectorXd::Zero(EtE.rows() + Cpart.rows());
        b.segment(0, EtE.rows())= Epart.transpose() * Mx * (rawField + torhs);
        VectorXd bfull = -Cfull * fixedValues;
        VectorXd bpart(CpartRank);
        for(int k = 0; k < CpartRank; k++)
            bpart(k)=bfull(PIndices(k));
        b.segment(EtE.rows(), Cpart.rows()) = bpart;
        
        SparseLU<SparseMatrix<double> > lusolver;
        lusolver.compute(A);
        if(lusolver.info() != Success){
            if (intData.verbose)
                cout<<"LU decomposition failed!"<<endl;
            return false;
        }
        x = lusolver.solve(b);
        
        fullx = var2AllMat * x.head(numVars - alreadyFixed.sum()) + fixedValues;
        
        
        if((alreadyFixed - fixedMask).sum() == 0)
            break;
        
        double minIntDiff = std::numeric_limits<double>::max();
        int minIntDiffIndex = -1;
        for (int i = 0; i < numVars; i++)
        {
            if ((fixedMask(i)) && (!alreadyFixed(i)))
            {
                double currIntDiff =0;
                double func = fullx(i); //fullx.segment(intData.d*i,intData.d);
                //for (int j=0;j<intData.d;j++)
                currIntDiff += std::fabs(func - std::round(func));
                if (currIntDiff < minIntDiff)
                {
                    minIntDiff = currIntDiff;
                    minIntDiffIndex = i;
                }
            }
        }
        /*
        if (minIntDiffIndex != -1)
        {
            alreadyFixed(minIntDiffIndex) = 1;
            double func = fullx(minIntDiffIndex) ;
            double funcInteger=std::round(func);
            fixedValues(minIntDiffIndex) = /*pinvSymm*projMat**///funcInteger;
    /*}
        //in case all singularities are rounded in the rounding-singularities mode, but there are left unrounded seams (like topological handles).
        if ((alreadyFixed.sum()==fixedMask.sum())&&(!intData.roundSeams)&(!roundedSingularities)&&(intData.integralSeamless)) {
            for (int i = 0; i < intData.integerVars.size(); i++)
                for (int j = 0; j < intData.n; j++)
                    fixedMask(intData.n * intData.integerVars(i) + j) = 1;
            roundedSingularities = true;
        }
        
        xprev.resize(x.rows() - 1);
        varCounter = 0;
        for(int i = 0; i < numVars; i++)
            if (!alreadyFixed(i))
                xprev(varCounter++) = fullx(i);
        
        xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());
    }*/
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd NFunctionVec = intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * xFull;
    NFunction.resize(meshCut.V.rows(), intData.N);
    for(int i = 0; i < NFunction.rows(); i++)
        NFunction.row(i) << NFunctionVec.segment(intData.N * i, intData.N).transpose();
    
    //nFunction = fullx;
    
    //allocating per corner
    NCornerFunctions.resize(meshWhole.F.rows(), intData.N*3);
    for (int i=0;i<meshWhole.F.rows();i++)
        for (int j=0;j<3;j++)
            NCornerFunctions.block(i, intData.N*j, 1, intData.N) = NFunction.row(meshCut.F(i,j));
    
    NFunctionVec = intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * xFull;
    NFunction.resize(meshCut.V.rows(), intData.N);
    for(int i = 0; i < NFunction.rows(); i++)
        NFunction.row(i) << NFunctionVec.segment(intData.N * i, intData.N).transpose();
    
    intData.nVertexFunction = xFull;
    
    
    //allocating per corner
    NCornerFunctions.resize(meshWhole.F.rows(), intData.N*3);
    for (int i=0;i<meshWhole.F.rows();i++)
        for (int j=0;j<3;j++)
            NCornerFunctions.block(i, intData.N*j, 1, intData.N) = NFunction.row(meshCut.F(i,j)).array();
    
    return success;
    
    
}

}

#endif


