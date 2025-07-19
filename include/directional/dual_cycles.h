//This file is part of Directional, a library for directional field processing.
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_DUAL_CYCLES_H
#define DIRECTIONAL_DUAL_CYCLES_H

#include <Eigen/Core>
#include <vector>
#include <set>
#include <unordered_map>
#include <numbers>
#include <directional/set_diff.h>
#include <directional/matrix_slice.h>
#include "tree.h"


namespace directional
{
// Creates the set of independent dual cycles (closed loops of connected faces that cannot be morphed to each other) on a mesh. Primarily used for index prescription.
// The basis cycle matrix first contains #V-#b cycles for every inner vertex (by order), then #b boundary cycles, and finally 2*g generator cycles around all handles. Total #c cycles.The cycle matrix sums information on the dual edges between the faces, and is indexed into the inner edges alone (excluding boundary)
//input:
//  V: #V by 3 vertices.
//  F: #F by 3 triangles.
//  EV: #E by 2 matrix of edges (vertex indices)
//  EF: #E by 2 matrix of oriented adjacent faces
//output:
//  basisCycles:    #c by #iE basis cycles
//  cycleCurvature:   #c by 1 curvatures of each cycle (for inner-vertex cycles, simply the Gaussian curvature.
//  vertex2cycle:     #v by 1 map between vertex and corresponding cycle (for comfort of input from the user's side; inner vertices map to their cycles, boundary vertices to the bigger boundary cycle.
//  innerEdges:       #iE by 1 the subset of #EV that are inner edges, and with the same ordering as the columns of basisCycles.

inline void dual_cycles(const TriMesh& mesh,
                        Eigen::SparseMatrix<double>& basisCycles,
                        Eigen::VectorXd& cycleCurvature,
                        Eigen::VectorXi& vertex2cycle,
                        Eigen::VectorXi& innerEdges)
{
    using namespace Eigen;
    using namespace std;
    int numV = mesh.F.maxCoeff() + 1;
    int eulerChar = numV - mesh.EV.rows() + mesh.F.rows();
    vertex2cycle.conservativeResize(mesh.V.rows());
    
    
    int numBoundaries=mesh.boundaryLoops.size();
    int numGenerators=2-numBoundaries-eulerChar;
    
    vector<Triplet<double> > basisCycleTriplets(mesh.EV.rows() * 2);
    
    //all 1-ring cycles, including boundaries
    for (int i = 0; i < mesh.EV.rows(); i++) {
        basisCycleTriplets[2*i]=Triplet<double>(mesh.EV(i, 0), i, -1.0);
        basisCycleTriplets[2*i+1]=Triplet<double>(mesh.EV(i, 1), i, 1.0);
    }
    
    //Creating boundary cycles by building a matrix the sums up boundary loops and zeros out boundary vertex cycles - it will be multiplied from the left to basisCyclesMat
    VectorXi isBoundary(mesh.V.rows()); isBoundary.setZero();
    for (int i=0;i<mesh.boundaryLoops.size();i++)
        for (int j=0;j<mesh.boundaryLoops[i].size();j++)
            isBoundary(mesh.boundaryLoops[i][j])=1;
    
    VectorXi pureInnerEdgeMask=VectorXi::Constant(mesh.EV.rows(),1);
    for (int i=0;i<mesh.EV.rows();i++)
        if ((isBoundary(mesh.EV(i,0)))||(isBoundary(mesh.EV(i,1))))
            pureInnerEdgeMask(i)=0;
    
    int currGeneratorCycle=0;
    int currBoundaryCycle=0;
    
    if ((numGenerators!=0)/*||(numBoundaries!=0)*/){
        MatrixXi reducedEV(mesh.EV);
        for (int i = 1; i < reducedEV.rows(); i++)
            if(isBoundary(reducedEV(i,0)) || isBoundary(reducedEV(i, 1)))
                reducedEV(i,0) = -1;
        
        VectorXi primalTreeEdges, primalTreeFathers;
        VectorXi dualTreeEdges, dualTreeFathers;
        tree(reducedEV, primalTreeEdges, primalTreeFathers);
        //creating a set of dual edges that do not cross edges in the primal tree
        VectorXi fullIndices = VectorXi::LinSpaced(mesh.EV.rows(), 0, mesh.EV.rows() - 1);
        VectorXi reducedEFIndices, inFullIndices;
        MatrixXi reducedEF;
        directional::set_diff(fullIndices, primalTreeEdges, reducedEFIndices, inFullIndices);
        VectorXi Two = VectorXi::LinSpaced(2, 0, 1);
        
        //cout<<"mesh.EF: "<<mesh.EF<<endl;
        //cout<<"reducedEFIndices: "<<reducedEFIndices<<endl;
        directional::matrix_slice(mesh.EF, reducedEFIndices, Two, reducedEF);
        //cout<<"reducedEF: "<<reducedEF<<endl;
        tree(reducedEF, dualTreeEdges, dualTreeFathers);
        //converting dualTreeEdges from reducedEF to EF
        for (int i = 0; i < dualTreeEdges.size(); i++)
            dualTreeEdges(i) = inFullIndices(dualTreeEdges(i));
        
        for (int i = 0; i < dualTreeFathers.size(); i++)
            if (dualTreeFathers(i) != -1 && dualTreeFathers(i) != -2)
                dualTreeFathers(i) = inFullIndices(dualTreeFathers(i));
        
        //building tree co-tree based homological cycles
        //finding dual edge which are not in the tree, and following their faces to the end
        VectorXi isinTree = VectorXi::Zero(mesh.EF.rows());
        for (int i = 0; i < dualTreeEdges.size(); i++) {
            isinTree(dualTreeEdges(i)) = 1;
        }
        for (int i = 0; i < primalTreeEdges.size(); i++) {
            isinTree(primalTreeEdges(i)) = 1;
        }
        
        for (int i = 0; i < isinTree.size(); i++) {
            if (isinTree(i))
                continue;
            
            //std::cout<<"New Cycle"<<std::endl;
            //otherwise, follow both end faces to the root and this is the dual cycle
            if (mesh.EF(i, 0) == -1 || mesh.EF(i, 1) == -1)
                continue;
            std::vector<Triplet<double> > candidateTriplets;
            //candidateTriplets.push_back(Triplet<double>(0, i, 1.0));
            Vector2i currLeaves; currLeaves << mesh.EF(i, 0), mesh.EF(i, 1);
            VectorXi visitedOnce = VectorXi::Zero(mesh.EF.rows());  //used to remove the tail from the LCA to the root
            bool isBoundaryCycle=true;
            for (int i = 0; i < 2; i++) { //on leaves
                int currTreeEdge = -1;  //indexing within dualTreeEdges
                int currFace = currLeaves(i);
                currTreeEdge = dualTreeFathers(currFace);
                if (currTreeEdge == -2)
                {
                    break;
                }
                
                while (currTreeEdge != -1) {
                    //std::cout<<"currTreeEdge: "<<currTreeEdge<<"\n"<<std::endl;
                    //determining orientation of current edge vs. face
                    double sign = ((mesh.EF(currTreeEdge, 0) == currFace) != (i == 0) ? 1.0 : -1.0);
                    visitedOnce(currTreeEdge) = 1 - visitedOnce(currTreeEdge);
                    candidateTriplets.push_back(Triplet<double>(0, currTreeEdge, sign));
                    currFace = (mesh.EF(currTreeEdge, 0) == currFace ? mesh.EF(currTreeEdge, 1) : mesh.EF(currTreeEdge, 0));
                    currTreeEdge = dualTreeFathers(currFace);
                }
            }
            
            //only putting in dual edges that are below the LCA
            for (int i=0;i<candidateTriplets.size();i++)
                if ((visitedOnce(candidateTriplets[i].col()))&&(pureInnerEdgeMask(candidateTriplets[i].col())))
                    isBoundaryCycle=false;
            
            if (isBoundaryCycle)
                continue; //ignoring those
            
            int currRow = (isBoundaryCycle ? numV+currBoundaryCycle : numV+currGeneratorCycle);
            (isBoundaryCycle ? currBoundaryCycle++ : currGeneratorCycle++);
            
            basisCycleTriplets.push_back(Triplet<double>(currRow, i, 1.0));
            for (size_t i = 0; i < candidateTriplets.size(); i++)
                if (visitedOnce(candidateTriplets[i].col())){
                    Triplet<double> trueTriplet(currRow,candidateTriplets[i].col(), candidateTriplets[i].value());
                    basisCycleTriplets.push_back(trueTriplet);
                }
        }
        //assert(currBoundaryCycle==numBoundaries && currGeneratorCycle==numGenerators);
    }
    
    numGenerators =currGeneratorCycle;
    
    SparseMatrix<double> sumBoundaryLoops(numV+numBoundaries+numGenerators,numV+numGenerators);
    vector<Triplet<double>> sumBoundaryLoopsTriplets;
    vector<int> innerVerticesList, innerEdgesList;
    VectorXi remainRows, remainColumns;
    
    for (int i=0;i<numV;i++){
        sumBoundaryLoopsTriplets.push_back(Triplet<double>(i, i,1.0-isBoundary[i]));
        if (!isBoundary(i)){
            innerVerticesList.push_back(i);
            vertex2cycle(i)=innerVerticesList.size()-1;
        }
    }
    
    for (int i=0;i<mesh.EV.rows();i++)
        if ((mesh.EF(i,0) != -1)&&(mesh.EF(i,1) != -1))
            innerEdgesList.push_back(i);
    
    //summing up boundary loops
    for (int i=0;i<mesh.boundaryLoops.size();i++)
        for (int j=0;j<mesh.boundaryLoops[i].size();j++){
            sumBoundaryLoopsTriplets.push_back(Triplet<double>(numV+i, mesh.boundaryLoops[i][j],1.0));
            vertex2cycle(mesh.boundaryLoops[i][j])=innerVerticesList.size()+i;
        }
    
    
    //just passing generators through;
    for (int i=numV;i<numV+numGenerators;i++)
        sumBoundaryLoopsTriplets.push_back(Triplet<double>(i, i,1.0));
    
    //Creating a matrix that aggregates basic cycles on boundary vertices
    sumBoundaryLoops.setFromTriplets(sumBoundaryLoopsTriplets.begin(), sumBoundaryLoopsTriplets.end());
    
    basisCycles.resize(numV+numGenerators, mesh.EV.rows());
    basisCycles.setFromTriplets(basisCycleTriplets.begin(), basisCycleTriplets.end());
    basisCycles=sumBoundaryLoops*basisCycles;
    
    //removing rows and columns of boundary vertices
    remainRows.resize(innerVerticesList.size()+numBoundaries+numGenerators);
    remainColumns.resize(innerEdgesList.size());
    for (int i=0;i<innerVerticesList.size();i++)
        remainRows(i)=innerVerticesList[i];
    
    for (int i=0;i<numBoundaries+numGenerators;i++)
        remainRows(innerVerticesList.size()+i)=numV+i;
    
    for (int i=0;i<innerEdgesList.size();i++)
        remainColumns(i)=innerEdgesList[i];
    
    //creating slicing matrices
    std::vector<Triplet<double> > rowSliceTriplets, colSliceTriplets;
    for (int i=0;i<remainRows.size();i++)
        rowSliceTriplets.push_back(Triplet<double>(i, remainRows(i), 1.0));
    for (int i=0;i<remainColumns.size();i++)
        colSliceTriplets.push_back(Triplet<double>(remainColumns(i), i, 1.0));
    
    SparseMatrix<double> rowSliceMat(remainRows.rows(), basisCycles.rows());
    rowSliceMat.setFromTriplets(rowSliceTriplets.begin(), rowSliceTriplets.end());
    
    SparseMatrix<double> colSliceMat(basisCycles.cols(), remainColumns.rows());
    colSliceMat.setFromTriplets(colSliceTriplets.begin(), colSliceTriplets.end());
    
    basisCycles=rowSliceMat*basisCycles*colSliceMat;
    
    innerEdges.conservativeResize(innerEdgesList.size());
    for (int i=0;i<innerEdgesList.size();i++)
        innerEdges(i)=innerEdgesList[i];
    
    //Correct computation of cycle curvature by adding angles
    //getting corner angle sum
    VectorXd allAngles(3*mesh.F.rows());
    for (int i=0;i<mesh.F.rows();i++){
        for (int j=0;j<3;j++){
            RowVector3d edgeVec12=mesh.V.row(mesh.F(i,(j+1)%3))-mesh.V.row(mesh.F(i,j));
            RowVector3d edgeVec13=mesh.V.row(mesh.F(i,(j+2)%3))-mesh.V.row(mesh.F(i,j));
            allAngles(3*i+j)=acos(edgeVec12.normalized().dot(edgeVec13.normalized()));
        }
    }
    
    //for each cycle, summing up all its internal angles negatively  + either 2*pi*|cycle| for internal cycles or pi*|cycle| for boundary cycles.
    cycleCurvature=VectorXd::Zero(basisCycles.rows());
    VectorXi isBigCycle=VectorXi::Ones(basisCycles.rows());  //TODO: retain it rather then reverse-engineer...
    
    for (int i=0;i<mesh.V.rows();i++)  //inner cycles
        if (!isBoundary(i))
            isBigCycle(vertex2cycle(i))=0;
    
    //getting the 4 corners of each edge to allocated later to cycles according to the sign of the edge.
    vector<set<int>> cornerSets(basisCycles.rows());
    vector<set<int>> vertexSets(basisCycles.rows());
    MatrixXi edgeCorners(innerEdges.size(),4);
    for (int i=0;i<innerEdges.rows();i++){
        int inFace1=0;
        while (mesh.F(mesh.EF(innerEdges(i),0),inFace1)!=mesh.EV(innerEdges(i),0))
            inFace1=(inFace1+1)%3;
        int inFace2=0;
        while (mesh.F(mesh.EF(innerEdges(i),1),inFace2)!=mesh.EV(innerEdges(i),1))
            inFace2=(inFace2+1)%3;
        
        edgeCorners(i,0)=mesh.EF(innerEdges(i),0)*3+inFace1;
        edgeCorners(i,1)=mesh.EF(innerEdges(i),1)*3+(inFace2+1)%3;
        edgeCorners(i,2)=mesh.EF(innerEdges(i),0)*3+(inFace1+1)%3;
        edgeCorners(i,3)=mesh.EF(innerEdges(i),1)*3+inFace2;
    }
    
    for (int k=0; k<basisCycles.outerSize(); ++k)
        for (SparseMatrix<double>::InnerIterator it(basisCycles,k); it; ++it){
            cornerSets[it.row()].insert(edgeCorners(it.col(),it.value()<0 ? 0 : 2));
            cornerSets[it.row()].insert(edgeCorners(it.col(),it.value()<0 ? 1 : 3));
            vertexSets[it.row()].insert(mesh.EV(innerEdges(it.col()), it.value()<0 ? 0 : 1));
        }
    
    for (int i=0;i<cornerSets.size();i++){
        if (isBigCycle(i))
            cycleCurvature(i)=std::numbers::pi*(double)(vertexSets[i].size());
        else
            cycleCurvature(i)=2.0*std::numbers::pi;
        for (set<int>::iterator si=cornerSets[i].begin();si!=cornerSets[i].end();si++)
            cycleCurvature(i)-=allAngles(*si);
    }
    
}
}

#endif
