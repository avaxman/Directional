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
#include <directional/setup_mesh_function_isolines.h>

namespace directional{


bool mesh_function_isolines(const Eigen::MatrixXd& origV,
                            const Eigen::MatrixXi& origF,
                            const Eigen::MatrixXi& EV,
                            const Eigen::MatrixXi& EF,
                            const Eigen::MatrixXi& FE,
                            const MeshFunctionIsolinesData& mfiData,
                            const bool verbose,
                            Eigen::MatrixXd& VOutput,
                            Eigen::VectorXi& DOutput,
                            Eigen::MatrixXi& FOutput){
  
  
  NFunctionMesher TMesh, FMesh;
  
  Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
  Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
  Eigen::MatrixXd FEsPoly;
  hedra::polygonal_edge_topology(Eigen::VectorXi::Constant(origF.rows(),3), origF,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
  hedra::dcel(Eigen::VectorXi::Constant(origF.rows(),3),origF,EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);
  
  TMesh.fromHedraDCEL(Eigen::VectorXi::Constant(origF.rows(),3),origV, origF, EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, mfiData.cutV, mfiData.cutF, mfiData.vertexNFunction,  mfiData.N, mfiData.orig2CutMat, mfiData.exactOrig2CutMat, mfiData.integerVars);
  
  if (verbose){
    std::cout<<"Generating mesh"<<std::endl;
    TMesh.GenerateMesh(FMesh);
    std::cout<<"Done generating!"<<std::endl;
  }
  
  Eigen::VectorXi genInnerEdges,genTF;
  Eigen::MatrixXi genEV,genEFi, genEF,genFE, genTEdges;
  Eigen::MatrixXd genFEs, genCEdges, genVEdges;
  
  if (verbose)
    std::cout<<"Cleaning Mesh"<<std::endl;
  
  bool success = FMesh.SimplifyMesh(verbose, mfiData.N);
  
  if (success){
    if (verbose)
      std::cout<<"Cleaning succeeded!"<<std::endl;
    
    FMesh.toHedra(VOutput,DOutput, FOutput);
  } else if (verbose) std::cout<<"Cleaning failed!"<<std::endl;
  
  return success;
  
  
}

} //namespace directional






#endif
