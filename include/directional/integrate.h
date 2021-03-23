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
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/bounding_box_diagonal.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/setup_integration.h>
#include <directional/branched_gradient.h>
#include <directional/iterative_rounding.h>
#include <igl/per_face_normals.h>

namespace directional
{
  
  // Creates an integration of N-functions from a directional field by solving the seamless Poisson equation
  // Input:
  //  wholeV:       #V x 3 vertex coordinates of the original mesh.
  //  wholeF:       #F x 3 face vertex indices of the original mesh.
  //  FE:           #F x 3 faces to edges indices.
  //  rawField:     #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size.
  //  lengthRatio   #edgeLength/bounding_box_diagonal of quad mesh (scaling the gradient).
  //  intData:           Integration data obtained from directional::setup_integration.
  //  cutV:         #cV x 3 vertices of the cut mesh.
  //  cutF:         #F x 3 faces of the cut mesh.
  //  roundIntegers;   which variables (from #V+#T) are rounded iteratively to double integers. for each "x" entry that means that the [4*x,4*x+4] entries of vt will be double integer.
  // Output:
  //  paramFuncsd             #cV x d parameterization functions per cut vertex (the compact version)
  //  paramFuncsN             #cV x N parameterization functions per cut vertex (full version with all symmetries unpacked)
  // wholeCornerParamFuncsN   (3*N) x #F parameterization functions per corner of whole mesh
  IGL_INLINE void integrate(const Eigen::MatrixXd& wholeV,
                            const Eigen::MatrixXi& wholeF,
                            const Eigen::MatrixXi& FE,
                            const Eigen::MatrixXd rawField,
                            const IntegrationData& intData,
                            const Eigen::MatrixXd& cutV,
                            const Eigen::MatrixXi& cutF,
                            Eigen::MatrixXd& paramFuncsd,
                            Eigen::MatrixXd& paramFuncsN,
                            Eigen::MatrixXd& wholeCornerParamFuncsN)
  
  
  {
    using namespace Eigen;
    using namespace std;
    
    VectorXd edgeWeights = VectorXd::Constant(FE.maxCoeff() + 1, 1.0);
    double length = igl::bounding_box_diagonal(wholeV) * intData.lengthRatio;
    
    int N = intData.N; //rawField.cols() / 3;
    int numVars = intData.symmMat.cols();
    //constructing face differentials
    vector<Triplet<double> >  d0Triplets;
    vector<Triplet<double> > M1Triplets;
    VectorXd gamma(3 * N * wholeF.rows());
    for(int i = 0; i < cutF.rows(); i++)
    {
      for(int j = 0; j < 3; j++)
      {
        for(int k = 0; k < N; k++)
        {
          d0Triplets.emplace_back(3 * N * i + N * j + k, N * cutF(i, j) + k, -1.0);
          d0Triplets.emplace_back(3 * N * i + N * j + k, N * cutF(i, (j + 1) % 3) + k, 1.0);
          Vector3d edgeVector = (cutV.row(cutF(i, (j + 1) % 3)) - cutV.row(cutF(i, j))).transpose();
          gamma(3 * N * i + N * j + k) = (rawField.block(i, 3 * k, 1, 3) * edgeVector)(0, 0) / length;
          M1Triplets.emplace_back(3 * N * i + N * j + k, 3 * N * i + N * j + k, edgeWeights(FE(i, j)));
        }
      }
    }
    SparseMatrix<double> d0(3 * N * wholeF.rows(), N * cutV.rows());
    d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
    SparseMatrix<double> M1(3 * N * wholeF.rows(), 3 * N * wholeF.rows());
    M1.setFromTriplets(M1Triplets.begin(), M1Triplets.end());
    SparseMatrix<double> d0T = d0.transpose();
    
    //creating face vector mass matrix
    std::vector<Triplet<double>> MxTri;
    VectorXd darea;
    igl::doublearea(cutV,cutF,darea);
    for (int i=0;i<cutF.rows();i++)
      for (int j=0;j<N;j++)
        for (int k=0;k<3;k++)
          MxTri.push_back(Triplet<double>(i*3*N+3*j+k,3*i*N+3*j+k,darea(i)/2.0));
    
    SparseMatrix<double> Mx(3*N*cutF.rows(), 3*N*cutF.rows());
    Mx.setFromTriplets(MxTri.begin(), MxTri.end());
    
    //The variables that should be fixed in the end
    VectorXi fixedMask(numVars);
    fixedMask.setZero();
    
    
    for (int i=0;i<intData.fixedIndices.size();i++)
      fixedMask(intData.fixedIndices(i)) = 1;
    
    if(false)
      for(int i = 0; i < intData.integerVars.size(); i++)
        for (int j=0;j<intData.n;j++)
          fixedMask(intData.n * intData.integerVars(i)+j) = 1;
    
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
    
    SparseMatrix<double> Efull = d0 * intData.vertexTrans2CutMat * intData.symmMat * intData.singIntSpanMat * intData.intSpanMat;
    VectorXd x, xprev;
    
    // until then all the N depedencies should be resolved?
    
    //reducing constraintMat
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
    SparseMatrix<double> Cfull = intData.constraintMat * intData.symmMat * intData.singIntSpanMat * intData.intSpanMat;
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
      if(lusolver.info() != Success)
        throw std::runtime_error("LU decomposition failed!");
      x = lusolver.solve(b);
      
      fullx = var2AllMat * x.head(numVars - alreadyFixed.sum()) + fixedValues;
      
      //#ifndef NDEBUG
      //cout << "(Cfull * fullx).lpNorm<Infinity>(): "<< (Cfull * fullx).lpNorm<Infinity>() << endl;
      cout << "Initial Poisson error: " << (Efull * fullx - gamma).lpNorm<Infinity>() << endl;
      //#endif
      
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
      
      //#ifndef NDEBUG
      //cout << "variable index: " << minIntDiffIndex << endl;
      //cout << "variable Integer error: " << minIntDiff << endl;
      //#endif
      
      if (minIntDiffIndex != -1)
      {
        alreadyFixed(minIntDiffIndex) = 1;
        double func = fullx(minIntDiffIndex) ;//.segment(intData.d*minIntDiffIndex,intData.d);
        //VectorXd funcInteger(intData.d);
        double funcInteger=std::round(func);
        //cout<<"fullx.segment(intData.d*minIntDiffIndex,intData.d): "<<fullx.segment(intData.d*minIntDiffIndex,intData.d)<<endl;
        //cout<<"func: "<<func<<endl;
        
        //for (int d=0;d<intData.d;d++)
        //  funcInteger(d)=std::round(func(d));
        //cout<<"funcInteger: "<<funcInteger<<endl;
        //fixedValues.segment(intData.d*minIntDiffIndex,intData.d) = /*pinvSymm*projMat**/funcInteger;
        fixedValues(minIntDiffIndex) = /*pinvSymm*projMat**/funcInteger;
        //cout<<" fixedValues.segment(intData.d*minIntDiffIndex,intData.d): "<< fixedValues.segment(intData.d*minIntDiffIndex,intData.d)<<endl;
      }
      
      xprev.resize(x.rows() - 1);
      varCounter = 0;
      for(int i = 0; i < numVars; i++)
        if (!alreadyFixed(i))
          xprev(varCounter++) = fullx(i);
      
      xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());
    }
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd paramFuncsVec = intData.vertexTrans2CutMat * intData.symmMat * intData.singIntSpanMat * intData.intSpanMat * fullx;
    paramFuncsN.conservativeResize(cutV.rows(), intData.N);
    for(int i = 0; i < paramFuncsN.rows(); i++)
      paramFuncsN.row(i) << paramFuncsVec.segment(intData.N * i, N).transpose();
    
    //cout<<"paramFuncsN: "<<paramFuncsN<<endl;
    
    paramFuncsd = fullx;
    
    //allocating per corner
    wholeCornerParamFuncsN.conservativeResize(wholeF.rows(), N*3);
    for (int i=0;i<wholeF.rows();i++)
      for (int j=0;j<3;j++)
        wholeCornerParamFuncsN.block(i, N*j, 1, N) = paramFuncsN.row(cutF(i,j));
    
    //cout<<"fullx: "<<fullx<<endl;
    
    //exporting things to MATLAB
    
    SparseMatrix<double> G;
    MatrixXd FN;
    igl::per_face_normals(cutV, cutF, FN);
    branched_gradient(cutV,cutF, intData.N, G);
    //cout<<"cutF.rows(): "<<cutF.rows()<<endl;
    SparseMatrix<double> Gd=G*intData.vertexTrans2CutMat * intData.symmMat * intData.singIntSpanMat * intData.intSpanMat;
    SparseMatrix<double> x2CornerMat=intData.vertexTrans2CutMat * intData.symmMat * intData.singIntSpanMat * intData.intSpanMat;
    //igl::matlab::MatlabWorkspace mw;
    VectorXi integerIndices(intData.integerVars.size()*intData.n);
    for(int i = 0; i < intData.integerVars.size(); i++)
      for (int j=0;j<intData.n;j++)
        integerIndices(intData.n * i+j) = intData.n * intData.integerVars(i)+j;
    
    //cout<<"integerIndices: "<<integerIndices<<endl;
    /*mw.save(Efull,"A");
     mw.save(fullx,"x0");
     mw.save(rawField,"rawField");
     mw.save_index(intData.fixedIndices,"fixedIndices");
     mw.save(intData.fixedValues,"fixedValues");
     mw.save_index(intData.singularIndices,"singularIndices");
     mw.save(gamma,"b");
     mw.save(Cfull,"C");
     mw.save(Gd,"G");
     //mw.save(Mx,"Mx");
     mw.save(FN,"FN");
     mw.save(intData.N,"N");
     mw.save(lengthRatio, "lengthRatio");
     //mw.save(length,"gradLength");
     mw.save(intData.d,"n");
     //mw.save(intData.intSpanMat,"intSpanMat");
     mw.save(cutV,"V");
     mw.save_index(cutF,"F");
     mw.save(x2CornerMat,"x2CornerMat");
     mw.save_index(integerIndices,"integerIndices");
     mw.write("poisson.mat");*/
    
    //working out MATLAB
    /*MatrixXi fixedIndicesMat(intData.fixedIndices.size(),1);
    fixedIndicesMat.col(0)=intData.fixedIndices;
    MatrixXd fixedValuesMat(intData.fixedValues.size(),1);
    fixedValuesMat.col(0)=intData.fixedValues;
    MatrixXi singularIndicesMat(intData.singularIndices.size(),1);
    singularIndicesMat.col(0)=intData.singularIndices;
    MatrixXi integerIndicesMat(integerIndices.size(),1);
    integerIndicesMat.col(0)=integerIndices;
    MatrixXd lengthRatioMat(1,1); lengthRatioMat(0,0)=lengthRatio;
    MatrixXd bMat(gamma.size(),1);
    bMat.col(0)=gamma;
    MatrixXd x0Mat(fullx.size(),1);
    x0Mat.col(0)=fullx;
    Engine* engine;
    igl::matlab::mlinit(&engine);
    igl::matlab::mlsetmatrix(&engine,"A",Efull);
    igl::matlab::mlsetmatrix(&engine,"rawField",rawField);
    igl::matlab::mlsetmatrix(&engine,"fixedIndices",fixedIndicesMat);
    igl::matlab::mlsetmatrix(&engine,"fixedValues",fixedValuesMat);
    igl::matlab::mlsetmatrix(&engine,"lengthRatio",lengthRatioMat);
    igl::matlab::mlsetmatrix(&engine,"b",bMat);
    igl::matlab::mlsetmatrix(&engine,"C",Cfull);
    igl::matlab::mlsetmatrix(&engine,"G",Gd);
    igl::matlab::mlsetmatrix(&engine,"FN",FN);
    igl::matlab::mlsetscalar(&engine,"N",intData.N);
    igl::matlab::mlsetscalar(&engine,"n",intData.d);
    igl::matlab::mlsetmatrix(&engine,"V",cutV);
    igl::matlab::mlsetmatrix(&engine,"F",cutF);
    igl::matlab::mlsetmatrix(&engine,"x2CornerMat",x2CornerMat);
    igl::matlab::mlsetmatrix(&engine,"x0",x0Mat);
    
    if (roundIntegers){
      igl::matlab::mlsetmatrix(&engine,"singularIndices",singularIndicesMat);
      igl::matlab::mlsetmatrix(&engine,"integerIndices",integerIndicesMat);
    } else{
      igl::matlab::mlsetmatrix(&engine,"singularIndices",Eigen::MatrixXi(0,0));
      igl::matlab::mlsetmatrix(&engine,"integerIndices",Eigen::MatrixXi(0,0));
      
    }
    
    std::string runLine=std::string("run('") + TUTORIAL_SHARED_PATH;
    if (roundSeams)
      igl::matlab::mleval(&engine,runLine + std::string("/../../include/directional/seamless_integration_seams' );"));
    else
      igl::matlab::mleval(&engine,runLine + std::string("/../../include/directional/seamless_integration_singularities' );"));*/
    
    bool success=directional::iterative_rounding(Efull, rawField, intData.fixedIndices, intData.fixedValues, intData.singularIndices, integerIndices, intData.lengthRatio, gamma, Cfull, Gd, FN, intData.N, intData.n, cutV, cutF, x2CornerMat,  intData.integralSeamless, intData.roundSeams, intData.localInjectivity, intData.verbose, fullx);
    
    
    /*MatrixXd fullxMat,successMat;
    igl::matlab::mlgetmatrix(&engine,"xcurr",fullxMat);
    igl::matlab::mlgetmatrix(&engine,"success",successMat);
    fullx=fullxMat.col(0);
    int success=(int)successMat(0,0);*/
    if (!success)
      cout<<"Rounding has failed!"<<endl;
    //igl::matlab::mlclose(&engine);
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    paramFuncsVec = intData.vertexTrans2CutMat * intData.symmMat * intData.singIntSpanMat * intData.intSpanMat * fullx;
    paramFuncsN.conservativeResize(cutV.rows(), intData.N);
    for(int i = 0; i < paramFuncsN.rows(); i++)
      paramFuncsN.row(i) << paramFuncsVec.segment(intData.N * i, N).transpose();
    
    paramFuncsd = fullx;
    
    //cout<<"paramFuncsd: "<<paramFuncsd<<endl;
    
    //allocating per corner
    wholeCornerParamFuncsN.conservativeResize(wholeF.rows(), N*3);
    for (int i=0;i<wholeF.rows();i++)
      for (int j=0;j<3;j++)
        wholeCornerParamFuncsN.block(i, N*j, 1, N) = paramFuncsN.row(cutF(i,j)).array();
    
    
  }
  
}

#endif


