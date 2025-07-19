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
#include <chrono>
#include <math.h>
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/polygonal_edge_topology.h>
#include <directional/setup_mesher.h>
#include <directional/NFunctionMesher.h>
#include <directional/generate_mesh.h>


namespace directional{


//Generates a mesh in (V,D,F) format from the integer isolines of a seamless N-function (such as the one computed from the Directional integrator). The mesh is polygonal, not necessarily triangular.
//Inputs:
//  origMesh:     the original whole mesh
//  mData:      a MesherData object that is pre-filled with the N-function data (should be generated from the integrator with setup_mesher)
//Outputs:
//  VOutput:      all vertex coordinates of the output polygonal mesh
//  DOutput:     |FOutput| vector of face valences
//  FOutput:      |FOutput| x |max(DOutput)| vertex indices of the face polygons, indexed into VOutput.
bool mesher(const directional::TriMesh& origMesh,
            const MesherData& mData,
            Eigen::MatrixXd& VOutput,
            Eigen::VectorXi& DOutput,
            Eigen::MatrixXi& FOutput){
    
    
    NFunctionMesher functionMesher(origMesh, mData); //TMesh, FMesh;
    functionMesher.init();
    
    if (mData.verbose)
        std::cout<<"Generating mesh"<<std::endl;
    
    functionMesher.generate_mesh();
    if (mData.verbose)
        std::cout<<"Done generating!"<<std::endl;
    
    
    Eigen::VectorXi genInnerEdges,genTF;
    Eigen::MatrixXi genEV,genEFi, genEF,genFE, genTEdges;
    Eigen::MatrixXd genFEs, genCEdges, genVEdges;
    
    bool success;
    if (mData.verbose){
        std::cout<<"Cleaning Mesh"<<std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        success = functionMesher.simplify_mesh();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout<<"Mesh simplification time: "<< duration.count()/1e+6 << " seconds" << std::endl;
    }else{
        success = functionMesher.simplify_mesh();
    }
    
    if (success){
        if (mData.verbose)
            std::cout<<"Cleaning succeeded!"<<std::endl;
        
        functionMesher.to_polygonal(VOutput,DOutput, FOutput);
    } else if (mData.verbose) std::cout<<"Cleaning failed!"<<std::endl;
    
    return success;
    
    
}

} //namespace directional






#endif
