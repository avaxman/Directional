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
#include <directional/TriMesh.h>
#include <directional/polygonal_edge_topology.h>
#include <directional/FunctionMesh.h>
#include <directional/setup_mesh_function_isolines.h>

namespace directional{


//Generates a mesh in (V,D,F) format from the integer isolines of a seamless N-function (such as the one computed from the Directional integrator). The mesh is polygonal, not necessarily triangular.
//Inputs:
//  origMesh:     the original whole mesh
//  mfiData:      a MeshFunctionIsolinesData object that is pre-filled with the N-function data (can be generated from the integrator with setup_mesh_function_isolines)
//  verbose:      if to output mesh generation process comments
//  VOutput:      all vertex coordinates of the output polygonal mesh
//  DOutput:     |FOutput| vector of face valences
//  FOutput:      |FOutput| x |max(DOutput)| vertex indices of the face polygons, indexed into VOutput.
bool mesh_function_isolines(const directional::TriMesh& origMesh,
                            const MeshFunctionIsolinesData& mfiData,
                            const bool verbose,
                            Eigen::MatrixXd& VOutput,
                            Eigen::VectorXi& DOutput,
                            Eigen::MatrixXi& FOutput){
  
  
  NFunctionMesher TMesh, FMesh;
  
  //Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
  //Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
  //Eigen::MatrixXd FEsPoly;
  //hedra::polygonal_edge_topology(Eigen::VectorXi::Constant(origMesh.F.rows(),3), origMesh.F,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
  //hedra::dcel(Eigen::VectorXi::Constant(origMesh.F.rows(),3),origMesh.F,EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);
  
  TMesh.fromHedraDCEL(Eigen::VectorXi::Constant(origMesh.F.rows(),3),origMesh.V, origMesh.F, EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, mfiData.cutV, mfiData.cutF, mfiData.vertexNFunction,  mfiData.N, mfiData.orig2CutMat, mfiData.exactOrig2CutMat, mfiData.integerVars);
  
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
