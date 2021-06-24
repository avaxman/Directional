// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef MESH_FUNCTION_HEADER_FILE
#define MESH_FUNCTION_HEADER_FILE


#include <iosfwd>
#include <vector>
#include <set>
#include <math.h>
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include <directional/FunctionMesh.h>

namespace Directional{


  bool mesh_function(const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const Eigen::MatrixXi& EV,
                     const Eigen::MatrixXi& EF,
                     const Eigen::MatrixXi& FE,
                     const Eigen::MatrixXd& nFunction,
                     const Eigen::SparseMatrix<int>& n2NMat,
                     const bool verbose,
                     Eigen::MatrixXd& VOutput,
                     Eigen::VectorXi& DOutput,
                     Eigen::MatrixXi& FOutput,
                     const double exactResolution=10e-9){
    
    
    NFunctionMesher TMesh, FMesh;
    
    Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
    Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
    Eigen::MatrixXd FEsPoly;
    hedra::polygonal_edge_topology(VectorXi::Constant(F.rows(),3), F,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
    hedra::dcel(VectorXi::Constant(F.rows(),3),F,EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);
    
    int capN=(N%2==0  ? N/2 : N);
    MatrixXd reducedCornerFuncs(FMeshWhole.rows(), NFull*3);
    for (int i=0;i<3;i++)
      reducedCornerFuncs.block(0,NFull*i,FMeshWhole.rows(),capN) =wholeCornerParamFuncsN.block(0,N*i,FMeshWhole.rows(),capN);
   
    TMesh.fromHedraDCEL(VectorXi::Constant(F.rows(),3),V, F, EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, VMeshCut, FMeshCut, paramFuncsd, intData.n, intData.N, intData.vertexTrans2CutMatInteger*pd.linRedInteger*intData.singIntSpanMatInteger*intData.intSpanMatInteger,  intData.constraintMatInteger*intData.linRedInteger*intData.singIntSpanMatInteger*intData.intSpanMatInteger, intData.linRed*intData.periodMat, NFunction, pd.integerVars, embNumMat, embDenMat, singVertices);
    
    if (verbose){
      std::cout<<"Generating mesh"<<std::endl;
      TMesh.GenerateMesh(FMesh);
      std::cout<<"Done generating!"<<std::endl;
    }
    
    Eigen::VectorXi genInnerEdges,genTF;
    Eigen::MatrixXi genEV,genEFi, genEF,genFE, genTEdges;
    Eigen::MatrixXd genFEs, genCEdges, genVEdges;
    
    if (verbose){
      std::cout<<"Simplifying Mesh"<<std::endl;
      bool success = FMesh.SimplifyMesh(N);
      if (success)
        std::cout<<"Simplification succeeded!"<<std::endl;
      else std::cout<<"Simplification failee!"<<std::endl;
    }
    
    FMesh.toHedra(V,D,  F, simpFfuncNum);
  }

} //namespace directional






#endif
