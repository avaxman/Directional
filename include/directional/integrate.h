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
#include <igl/igl_inline.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/bounding_box_diagonal.h>
#include <directional/tree.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/principal_matching.h>
#include <directional/setup_integration.h>
#include <directional/branched_gradient.h>
#include <directional/iterative_rounding.h>


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
    IGL_INLINE bool integrate(const directional::CartesianField& field,
                              IntegrationData& intData,
                              const directional::TriMesh& meshCut,
                              Eigen::MatrixXd& NFunction,
                              Eigen::MatrixXd& NCornerFunctions)


    {
        using namespace Eigen;
        using namespace std;

        assert(field.tb->discTangType()==discTangTypeEnum::FACE_SPACES && "Integrate() only works with face-based fields");
        const directional::TriMesh& meshWhole = *((IntrinsicFaceTangentBundle*)field.tb)->mesh;

        VectorXd edgeWeights = VectorXd::Constant(meshWhole.FE.maxCoeff() + 1, 1.0);
        //double length = igl::bounding_box_diagonal(wholeV) * intData.lengthRatio;

        int numVars = intData.linRedMat.cols();
        //constructing face differentials
        vector<Triplet<double> >  d0Triplets;
        vector<Triplet<double> > M1Triplets;
        VectorXd gamma(3 * intData.N * meshWhole.F.rows());
        for(int i = 0; i < meshCut.F.rows(); i++)
        {
            for(int j = 0; j < 3; j++)
            {
                for(int k = 0; k < intData.N; k++)
                {
                    d0Triplets.emplace_back(3 * intData.N * i + intData.N * j + k, intData.N * meshCut.F(i, j) + k, -1.0);
                    d0Triplets.emplace_back(3 * intData.N * i + intData.N * j + k, intData.N * meshCut.F(i, (j + 1) % 3) + k, 1.0);
                    Vector3d edgeVector = (meshCut.V.row(meshCut.F(i, (j + 1) % 3)) - meshCut.V.row(meshCut.F(i, j))).transpose();
                    gamma(3 * intData.N * i + intData.N * j + k) = (field.extField.block(i, 3 * k, 1, 3) * edgeVector)(0, 0);
                    M1Triplets.emplace_back(3 * intData.N * i + intData.N * j + k, 3 * intData.N * i + intData.N * j + k, edgeWeights(meshWhole.FE(i, j)));
                }
            }
        }
        SparseMatrix<double> d0(3 * intData.N * meshWhole.F.rows(), intData.N * meshCut.V.rows());
        d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
        SparseMatrix<double> M1(3 * intData.N * meshWhole.F.rows(), 3 * intData.N *  meshWhole.F.rows());
        M1.setFromTriplets(M1Triplets.begin(), M1Triplets.end());
        SparseMatrix<double> d0T = d0.transpose();

        //creating face vector mass matrix
        std::vector<Triplet<double>> MxTri;
        VectorXd darea;
        igl::doublearea(meshCut.V,meshCut.F,darea);
        for (int i=0;i<meshCut.F.rows();i++)
            for (int j=0;j<intData.N;j++)
                for (int k=0;k<3;k++)
                    MxTri.push_back(Triplet<double>(i*3*intData.N+3*j+k,3*i*intData.N+3*j+k,darea(i)/2.0));

        SparseMatrix<double> Mx(3*intData.N*meshCut.F.rows(), 3*intData.N*meshCut.F.rows());
        Mx.setFromTriplets(MxTri.begin(), MxTri.end());

        //The variables that should be fixed in the end
        VectorXi fixedMask(numVars);
        fixedMask.setZero();


        for (int i=0;i<intData.fixedIndices.size();i++)
            fixedMask(intData.fixedIndices(i)) = 1;

        /*if(false)
          for(int i = 0; i < intData.integerVars.size(); i++)
            for (int j=0;j<intData.n;j++)
              fixedMask(intData.n * intData.integerVars(i)+j) = 1;*/

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

        SparseMatrix<double> Efull = d0 * intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        VectorXd x, xprev;

        // until then all the N depedencies should be resolved?

        //reducing constraintMat
        SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
        SparseMatrix<double> Cfull = intData.constraintMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        VectorXd fullx(numVars);
        fullx.setZero();

        //the results are packets of N functions for each vertex, and need to be allocated for corners
        VectorXd NFunctionVec = intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * fullx;
        NFunction.resize(meshCut.V.rows(), intData.N);
        for(int i = 0; i < NFunction.rows(); i++)
            NFunction.row(i) << NFunctionVec.segment(intData.N * i, intData.N).transpose();

        //nFunction = fullx;

        //allocating per corner
        NCornerFunctions.resize(meshWhole.F.rows(), intData.N*3);
        for (int i=0;i<meshWhole.F.rows();i++)
            for (int j=0;j<3;j++)
                NCornerFunctions.block(i, intData.N*j, 1, intData.N) = NFunction.row(meshCut.F(i,j));

        SparseMatrix<double> G;
        //MatrixXd FN;
        //igl::per_face_normals(cutV, meshCut, FN);
        branched_gradient(meshCut.V,meshCut.F, intData.N, G);
        //cout<<"cutF.rows(): "<<cutF.rows()<<endl;
        SparseMatrix<double> Gd=G*intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        SparseMatrix<double> x2CornerMat=intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        //igl::matlab::MatlabWorkspace mw;
        VectorXi integerIndices(intData.integerVars.size()*intData.n);
        for(int i = 0; i < intData.integerVars.size(); i++)
            for (int j=0;j<intData.n;j++)
                integerIndices(intData.n * i+j) = intData.n * intData.integerVars(i)+j;


        bool success=directional::iterative_rounding(Efull, field.extField, intData.fixedIndices, intData.fixedValues, intData.singularIndices, integerIndices, intData.lengthRatio, gamma, Cfull, Gd, meshCut.faceNormals, intData.N, intData.n, meshCut.V, meshCut.F, x2CornerMat,  intData.integralSeamless, intData.roundSeams, intData.localInjectivity, intData.verbose, fullx);


        if ((!success)&&(intData.verbose))
            cout<<"Rounding has failed!"<<endl;

        //the results are packets of N functions for each vertex, and need to be allocated for corners
        NFunctionVec = intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat * fullx;
        NFunction.resize(meshCut.V.rows(), intData.N);
        for(int i = 0; i < NFunction.rows(); i++)
            NFunction.row(i) << NFunctionVec.segment(intData.N * i, intData.N).transpose();

        intData.nVertexFunction = fullx;

        //nFunction = fullx;

        //cout<<"paramFuncsd: "<<paramFuncsd<<endl;

        //allocating per corner
        NCornerFunctions.resize(meshWhole.F.rows(), intData.N*3);
        for (int i=0;i<meshWhole.F.rows();i++)
            for (int j=0;j<3;j++)
                NCornerFunctions.block(i, intData.N*j, 1, intData.N) = NFunction.row(meshCut.F(i,j)).array();

        return success;


    }

}

#endif


