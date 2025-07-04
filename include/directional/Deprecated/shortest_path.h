// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SHORTEST_PATH_H
#define DIRECTIONAL_SHORTEST_PATH_H

#include <vector>
#include <set>
#include <queue>
#include <Eigen/Core>
#include <directional/TriMesh.h>

namespace directional
{

/***
 Computing shortest paths from a given vertex
 Input:
 mesh:              On which this gets computed
 initialVertex:     the initial vertex index (into mesh.V)
 targetNodes:   
 ***/

void inline shortest_path(const TriMesh& mesh,
                          const int initialVertex,
                          std::vector<int> targetNodes,
                          std::vector<int> shortestPath){
    
    using namespace std;
    using namespace Eigen;
    double infWeight = mesh.avgEdgeLength*(double)mesh.dcel.halfedges.size()*10000.0;
    Eigen::VectorXd distances = Eigen::VectorXd::Constant(mesh.V.rows(), infWeight);
    int currVertex = 0;
    distances(currVertex)=0.0;
    Eigen::VectorXi isVisited = Eigen::VectorXi::Zero(mesh.V.rows());
    Eigen::VectorXi targetNodeMask = Eigen::VectorXi::Zero(mesh.V.rows());;
    for (int i=0;i<targetNodes.size();i++)  //for easy access of having been visiting
        targetNodeMask(targetNodes[i]) = 1;
    
    if (targetNodeMask[initialVertex]==1)
        return; //initialvertex is already a target vertex
    
    
    int sumVisited = 0;
    set<pair<double, pair<int,int>> > distVertices;
    distVertices.insert(pair<double, pair<int, int>>(0.0,pair<int, int>(currVertex,-1)));
    Eigen::VectorXi predHEList = Eigen::VectorXi::Constant(mesh.V.rows(),-1);
    
    while ((sumVisited<=mesh.V.rows())&&(targetNodeMask[currVertex]==0)){
        //calculating distance to all neighbors
        assert( (!isVisited[currVertex]) && "Trying an already visited vertex!");
        int hebegin = mesh.dcel.vertices[currVertex].halfedge;
        int hecurr = hebegin;
        do{
            int nextVertex =  mesh.dcel.halfedges[mesh.dcel.halfedges[hecurr].next].vertex;
            if (!isVisited(nextVertex)) {
                double length = (mesh.V.row(nextVertex) - mesh.V.row(currVertex)).norm();
                if (distances(nextVertex)>distances(currVertex)+length)
                    distances(nextVertex) = distances(currVertex)+length;
                
                distVertices.insert(pair<double, pair<int, int>>(distances(nextVertex),pair<int,int>(nextVertex, hecurr)));
            }
            if (mesh.dcel.halfedges[hecurr].twin==-1)
                break;
            hecurr = mesh.dcel.halfedges[mesh.dcel.halfedges[hecurr].prev].twin;
        }while (hecurr!=hebegin);
        isVisited[currVertex]=1;
        sumVisited++;
        if (sumVisited>=mesh.V.rows())
            break;  //shouldn't happen! this means we haven't reached any target node!
        //going to next vertex with shortest distance
        int predHE = -1;
        do{
            currVertex = distVertices.begin()->second.first;
            predHE =distVertices.begin()->second.second;
            distVertices.erase(distVertices.begin());
        }while ((isVisited[currVertex])&&(!distVertices.empty()));
        predHEList[currVertex]=predHE;
    }
    
    //Backtracking the path
    assert(targetNodeMask[currVertex]==1 && "Haven't reached any target node!");
    while (currVertex!=initialVertex){
        shortestPath.push_back(predHEList[currVertex]);
        currVertex = mesh.dcel.halfedges[predHEList[currVertex]].vertex;
    }
    
}
}

#endif
