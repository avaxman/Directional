// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PARAMETERIZE_H
#define DIRECTIONAL_PARAMETERIZE_H

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
#include <igl/matlab_format.h>
#include <igl/bounding_box_diagonal.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/setup_parameterization.h>


namespace directional
{

  // Creates a parameterization of (currently supported) (u,v, -u,-v) functions from a directional field by solving the Poisson equation, with custom edge weights
  // Input:
  //  wholeV:       #V x 3 vertex coordinates of the original mesh.
  //  wholeF:       #F x 3 face vertex indices of the original mesh.
  //  FE:           #F x 3 faces to edges indices.
  //  rawField:     #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size.
  //  lengthRatio   #edgeLength/bounding_box_diagonal of quad mesh (scaling the gradient).
  //  pd:           Parameterization data obtained from directional::setup_parameterization.
  //  cutV:         #cV x 3 vertices of the cut mesh.
  //  cutF:         #F x 3 faces of the cut mesh.
  //  roundIntegers;   which variables (from #V+#T) are rounded iteratively to double integers. for each "x" entry that means that the [4*x,4*x+4] entries of vt will be double integer.
  // Output:
  //  paramFuncsd             #cV x d parameterization functions per cut vertex (the compact version)
  //  paramFuncsN             #cV x N parameterization functions per cut vertex (full version with all symmetries unpacked)
  // wholeCornerParamFuncsN   (3*N) x #F parameterization functions per corner of whole mesh
  IGL_INLINE void parameterize(const Eigen::MatrixXd& wholeV,
                               const Eigen::MatrixXi& wholeF,
                               const Eigen::MatrixXi& FE,
                               const Eigen::MatrixXd rawField,
                               const double lengthRatio,
                               const ParameterizationData& pd,
                               const Eigen::MatrixXd& cutV,
                               const Eigen::MatrixXi& cutF,
                               const bool roundIntegers,
                               Eigen::MatrixXd& paramFuncsd,
                               Eigen::MatrixXd& paramFuncsN,
                               Eigen::MatrixXd& wholeCornerParamFuncsN)
  
  
  {
    using namespace Eigen;
    using namespace std;
    
    VectorXd edgeWeights = VectorXd::Constant(FE.maxCoeff() + 1, 1.0);
    double length = igl::bounding_box_diagonal(wholeV) * lengthRatio;
    
    //creating projection operator for N-packet. This is probably extremely redundant process....
    /*MatrixXd pinvSymm = (pd.symmFunc.transpose()*pd.symmFunc).inverse()*pd.symmFunc.transpose();
    MatrixXd symmConst = pd.symmFunc*pinvSymm-MatrixXd::Identity(pd.symmFunc.rows(),pd.symmFunc.rows());
    
    ColPivHouseholderQR<MatrixXd> denseqrsolver(symmConst.transpose());
    denseqrsolver.compute(symmConst.transpose());
    int SRank = denseqrsolver.rank();

    //creating sliced permutation matrix
    VectorXi PIndices = denseqrsolver.colsPermutation().indices();
    //cout<<"PIndices: "<<PIndices<<endl;
    MatrixXd symmConstReduced(SRank, symmConst.cols());
    for (int i=0;i<SRank;i++)
      symmConstReduced.row(i) =symmConst.row(PIndices(i));

    MatrixXd projMat = symmConstReduced.transpose()*(symmConstReduced*symmConstReduced.transpose()).inverse()*symmConstReduced+MatrixXd::Identity(pd.N, pd.N);*/
    
  
    //cout<<"symmConst: "<<symmConst<<endl;
    //cout<<"symmConst: "<<symmConst<<endl;
    //cout<<"symmConstReduced: "<<symmConstReduced<<endl;
    //cout<<"projMat: "<<projMat<<endl;
    
    //TODO: in vertex space, not corner...
    int N = pd.N; //rawField.cols() / 3;
    int numVars = pd.symmMat.cols()/pd.d;
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


    //the variables that should be fixed in the end
    VectorXi fixedMask(numVars);
    fixedMask.setZero();

    fixedMask(0) = 1;  //first vertex is always (0,0)

    if(roundIntegers)
      for(int i = 0; i < pd.integerVars.size(); i++)
        fixedMask(pd.integerVars(i)) = 1;
    
    
    //the variables that were already fixed in the previous iteration
    VectorXi alreadyFixed(numVars);
    alreadyFixed.setZero();

    alreadyFixed(0) = 1;  //first vertex is always (0,0)

    //the values for the fixed variables (size is as all variables)
    VectorXd fixedValues(pd.d*numVars);
    fixedValues.setZero();
    
    SparseMatrix<double> Efull = d0 * pd.vertexTrans2CutMat * pd.symmMat;
    VectorXd x, xprev;

    // until then all the N depedencies should be resolved?

    //reducing constraintMat
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
    SparseMatrix<double> Cfull = pd.constraintMat * pd.symmMat;
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
    SparseMatrix<double> var2AllMat;
    VectorXd fullx(pd.d*numVars); fullx.setZero();
    for(int intIter = 0; intIter < fixedMask.sum(); intIter++)
    {
      //the non-fixed variables to all variables
      var2AllMat.resize(pd.d*numVars, pd.d*numVars - pd.d*alreadyFixed.sum());
      int varCounter = 0;
      vector<Triplet<double> > var2AllTriplets;
      for(int i = 0; i < numVars; i++)
      {
        if (!alreadyFixed(i)){
          for (int j=0;j<pd.d;j++)
            var2AllTriplets.emplace_back(pd.d*i+j, varCounter++, 1.0);
        }
          
      }
      var2AllMat.setFromTriplets(var2AllTriplets.begin(), var2AllTriplets.end());
      
      SparseMatrix<double> Epart = Efull * var2AllMat;
      VectorXd torhs = -Efull * fixedValues;
      SparseMatrix<double> EtE = Epart.transpose() * M1 * Epart;
      SparseMatrix<double> Cpart = Cfull * var2AllMat;
      
      //reducing rank on Cpart
      qrsolver.compute(Cpart.transpose());
      int CpartRank = qrsolver.rank();
  
      //creating sliced permutation matrix
      VectorXi PIndices = qrsolver.colsPermutation().indices();

      vector<Triplet<double> > CPartTriplets;
      VectorXd bpart(CpartRank);
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
      for(int k = 0; k < CpartRank; k++)
        bpart(k)=bfull(PIndices(k));
      b.segment(EtE.rows(), Cpart.rows()) = bpart;

      SparseLU<SparseMatrix<double> > lusolver;
      lusolver.compute(A);
      if(lusolver.info() != Success)
        throw std::runtime_error("LU decomposition failed!");
      x = lusolver.solve(b);

      fullx = var2AllMat * x.head(pd.d*numVars - pd.d*alreadyFixed.sum()) + fixedValues;

//#ifndef NDEBUG
      cout << "(Cfull * fullx).lpNorm<Infinity>(): "<< (Cfull * fullx).lpNorm<Infinity>() << endl;
      cout << "Poisson error: " << (Efull * fullx - gamma).lpNorm<Infinity>() << endl;
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
          VectorXd func = fullx.segment(pd.d*i,pd.d);
          for (int j=0;j<pd.d;j++)
            currIntDiff += std::fabs(0.5 * fullx(i) - std::round(0.5 * fullx(j)));
          if (currIntDiff < minIntDiff)
          {
            minIntDiff = currIntDiff;
            minIntDiffIndex = i;
          }
        }
      }

//#ifndef NDEBUG
      cout << "d-segment index: " << minIntDiffIndex << endl;
      cout << "d-segment Integer error: " << minIntDiff/pd.d << endl;
//#endif
      
      if (minIntDiffIndex != -1)
      {
        alreadyFixed(minIntDiffIndex) = 1;
        VectorXd func = fullx.segment(pd.d*minIntDiffIndex,pd.d);
        VectorXd funcInteger(pd.d);
        //cout<<"fullx.segment(pd.d*minIntDiffIndex,pd.d): "<<fullx.segment(pd.d*minIntDiffIndex,pd.d)<<endl;
        //cout<<"func: "<<func<<endl;
       
        for (int d=0;d<pd.d;d++)
          funcInteger(d)=std::round(0.5*func(d))*2.0;
         //cout<<"funcInteger: "<<funcInteger<<endl;
        fixedValues.segment(pd.d*minIntDiffIndex,pd.d) = /*pinvSymm*projMat**/funcInteger;
        //cout<<" fixedValues.segment(pd.d*minIntDiffIndex,pd.d): "<< fixedValues.segment(pd.d*minIntDiffIndex,pd.d)<<endl;
      }
      
      xprev.resize(x.rows() - pd.d);
      varCounter = 0;
      for(int i = 0; i < numVars; i++)
        if (!alreadyFixed(i))
          xprev.segment(pd.d*varCounter++,pd.d) = fullx.segment(pd.d*i,pd.d);
      
      xprev.tail(Cpart.rows()) = x.tail(Cpart.rows());
    }

    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd paramFuncsVec = pd.vertexTrans2CutMat * pd.symmMat * fullx;
    paramFuncsN.conservativeResize(cutV.rows(), pd.N);
    for(int i = 0; i < paramFuncsN.rows(); i++)
      paramFuncsN.row(i) << paramFuncsVec.segment(pd.N * i, N).transpose();
    
    //cout<<"paramFuncsN: "<<paramFuncsN<<endl;
    
    paramFuncsd = fullx;
    
    //allocating per corner
    wholeCornerParamFuncsN.conservativeResize(wholeF.rows(), N*3);
    for (int i=0;i<wholeF.rows();i++)
      for (int j=0;j<3;j++)
        wholeCornerParamFuncsN.block(i, N*j, 1, N) = paramFuncsN.row(cutF(i,j));
    
    //cout<<"fullx: "<<fullx<<endl;
  }




}
  
#endif

  
