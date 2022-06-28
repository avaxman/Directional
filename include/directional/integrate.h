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
#include <directional/FaceField.h>
#include <directional/principal_matching.h>
#include <directional/setup_integration.h>
#include <directional/branched_gradient.h>
#include <directional/iterative_rounding.h>


namespace directional
{
  
  // Integrates an N-directional fields into an N-function by solving the seamless Poisson equation. Respects *valid* linear reductions where the field is reducible to an n-field for n<=M, and consequently the function is reducible to an n-function.
  // Input:
  //  wholeV:       #V x 3 vertex coordinates of the original mesh.
  //  wholeF:       #F x 3 face vertex indices of the original mesh.
  //  FE:           #F x 3 faces to edges indices.
  //  rawField:     #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size.
  //  intData:      Integration data, which must be obtained from directional::setup_integration.
  //  cutV:         #cV x 3 vertices of the cut mesh.
  //  cutF:         #F x 3 faces of the cut mesh.
  // Output:
  //  NFunction:    #cV x N parameterization functions per cut vertex (full version with all symmetries unpacked)
  // wholeCornerParamFuncsN   (3*N) x #F parameterization functions per corner of whole mesh
  IGL_INLINE bool integrate(const directional::TriMesh& meshWhole,
                            const directional::FaceField& field,
                            IntegrationData& intData,
                            const directional::TriMesh& meshCut,
                            Eigen::MatrixXd& NFunction,
                            Eigen::MatrixXd& NCornerFunctions)
  
  
  {
    using namespace Eigen;
    using namespace std;
    
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
      SparseMatrix<double> EtE = Epart.transpose() * M1 * Epart;
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
      b.segment(0, EtE.rows())= Epart.transpose() * M1 * (gamma + torhs);
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

      if (minIntDiffIndex != -1)
      {
        alreadyFixed(minIntDiffIndex) = 1;
        double func = fullx(minIntDiffIndex) ;
        double funcInteger=std::round(func);
        fixedValues(minIntDiffIndex) = /*pinvSymm*projMat**/funcInteger;
      }
      
      xprev.resize(x.rows() - 1);
      varCounter = 0;
      for(int i = 0; i < numVars; i++)
        if (!alreadyFixed(i))
          xprev(varCounter++) = fullx(i);
      
      xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());
    }
    
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


