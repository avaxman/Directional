// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef MESH_FUNCTION_ISOLINES_HEADER_FILE
#define MESH_FUNCTION_ISOLINES_HEADER_FILE


#include <iosfwd>
#include <vector>
#include <set>
#include <math.h>
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include <directional/polygonal_edge_topology.h>
#include <directional/FunctionMesh.h>

namespace directional{


bool mesh_function_isolines(const Eigen::MatrixXd& wholeV,
                            const Eigen::MatrixXi& wholeF,
                            const Eigen::MatrixXi& EV,
                            const Eigen::MatrixXi& EF,
                            const Eigen::MatrixXi& FE,
                            const Eigen::VectorXd& vertexNFunction,
                            const int N,
                            const Eigen::MatrixXd& cutV,
                            const Eigen::MatrixXi& cutF,
                            const Eigen::SparseMatrix<double>& vertex2CornerMat,
                            const Eigen::SparseMatrix<int>& exactVertex2CornerMat,
                            const Eigen::VectorXi& integerVars,
                            const bool verbose,
                            Eigen::MatrixXd& VOutput,
                            Eigen::VectorXi& DOutput,
                            Eigen::MatrixXi& FOutput,
                            const double exactResolution=10e-9){
  
  
  NFunctionMesher TMesh, FMesh;
  
  Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
  Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
  Eigen::MatrixXd FEsPoly;
  hedra::polygonal_edge_topology(Eigen::VectorXi::Constant(wholeF.rows(),3), wholeF,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
  hedra::dcel(Eigen::VectorXi::Constant(wholeF.rows(),3),wholeF,EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);
  
  TMesh.fromHedraDCEL(Eigen::VectorXi::Constant(wholeF.rows(),3),wholeV, wholeF, EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, cutV, cutF, vertexNFunction,  N, vertex2CornerMat, exactVertex2CornerMat, integerVars);
  
  if (verbose){
    std::cout<<"Generating mesh"<<std::endl;
    TMesh.GenerateMesh(FMesh);
    std::cout<<"Done generating!"<<std::endl;
  }
  
  Eigen::VectorXi genInnerEdges,genTF;
  Eigen::MatrixXi genEV,genEFi, genEF,genFE, genTEdges;
  Eigen::MatrixXd genFEs, genCEdges, genVEdges;
  
  if (verbose)
    std::cout<<"Simplifying Mesh"<<std::endl;
  
  bool success = true ;//FMesh.SimplifyMesh(verbose, N);
  
  if (success){
    if (verbose)
      std::cout<<"Simplification succeeded!"<<std::endl;
    
    FMesh.toHedra(VOutput,DOutput, FOutput);
  } else if (verbose) std::cout<<"Simplification failed!"<<std::endl;
  
  return success;
  
  
}

} //namespace directional






#endif
