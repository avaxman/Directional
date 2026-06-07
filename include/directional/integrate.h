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
    

    //The variables that should be fixed aprior
    VectorXi fixedMask(numVars);
    fixedMask.setZero();
    
    for (int i=0;i<intData.fixedIndices.size();i++)
        fixedMask(intData.fixedIndices(i)) = 1;
    
    //the values for the fixed variables (size is as all variables)
    VectorXd fixedValues(numVars);
    fixedValues.setZero();  //for everything but the originally fixed values
    for (int i=0;i<intData.fixedValues.size();i++)
        fixedValues(intData.fixedIndices(i))=intData.fixedValues(i);
    
    
    //The variables that need to be rounded
    VectorXi toRoundMask(numVars);
    toRoundMask.setZero();
    
    
    bool roundedSingularities = false;  //if all singularities have been rounded (only relevant to intData.roundSeams=false)
    if(intData.integralSeamless) {
        if (intData.roundSeams) {
            for (int i = 0; i < intData.integerVars.size(); i++)
                for (int j = 0; j < intData.n; j++)
                    toRoundMask(intData.n * intData.integerVars(i) + j) = 1;
        }else {
            for (int i = 0; i < intData.singularIndices.size(); i++)
                toRoundMask(intData.singularIndices(i)) = 1;
        }
    }
    
    
    directional::MixedIntegerSolver mis;
    
    mis.A = G * /** intData.featureAlignMat * */intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
    mis.M = Mx;
    mis.fixedMask = fixedMask;
    mis.fixedValues = fixedValues;
    mis.toRoundMask = toRoundMask;
    mis.C = intData.constraintMat;// * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
    if (intData.verbose)
        std::cout<<"Number of constraints: "<<mis.C.rows()<<std::endl;
    mis.b = rawField;
    mis.numVars = numVars;
    mis.verbose = intData.verbose;
    bool success = mis.solve();
    if (!success)
        return success;
    VectorXd xFull= mis.x;
    
    //solve again for extra integers ¯for topological unrounded seams
    //Some of these might have been rounded already, in which case the solver ignores them
    if ((!intData.roundSeams)&&(!roundedSingularities)&&(intData.integralSeamless)) {
        for (int i = 0; i < intData.integerVars.size(); i++)
            for (int j = 0; j < intData.n; j++)
                mis.toRoundMask(intData.n * intData.integerVars(i) + j) = 1;
        
        roundedSingularities = true;
    }
    if (intData.verbose)
        std::cout<<"Rounding extra unrounded topological seams if needed"<<std::endl;
    
    success = mis.solve();
    if (!success)
        return success;
    xFull= mis.x;
    
    std::cout<<"intData.constraintMat * xFull: "<<(intData.constraintMat * xFull).lpNorm<Infinity>()<<std::endl;
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd NFunctionVec = /*intData.featureAlignMat * */intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * xFull;
    NFunction.resize(meshCut.V.rows(), intData.N);
    for(int i = 0; i < NFunction.rows(); i++)
        NFunction.row(i) << NFunctionVec.segment(intData.N * i, intData.N).transpose();
    
    //Checking consistency
    //std::cout<<"intData.constraintMat * NFunctionVec: "<<intData.constraintMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * xFull<<std::endl;
    
    //std::cout<<"NFunctionVec: "<<NFunctionVec<<std::endl;
    VectorXd errorVec = intData.featureAlignConstMat * NFunctionVec;
    std::cout<<"intData.featureAlignConstMat * NFunctionVec: "<<errorVec.lpNorm<Infinity>()<<std::endl;
    
    //nFunction = fullx;
    
    //allocating per corner
    NCornerFunctions.resize(meshWhole.F.rows(), intData.N*3);
    for (int i=0;i<meshWhole.F.rows();i++)
        for (int j=0;j<3;j++)
            NCornerFunctions.block(i, intData.N*j, 1, intData.N) = NFunction.row(meshCut.F(i,j));
    
    NFunctionVec = intData.featureAlignMat * intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * xFull;
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


