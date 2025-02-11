// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef FUNCTION_MESH_CLASS_HEADER_FILE
#define FUNCTION_MESH_CLASS_HEADER_FILE

#include <set>
#include <math.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <directional/dcel.h>
#include <directional/GMP_definitions.h>
#include <directional/NFunctionMesher.h>
//#include <operators/io_stream.hpp>


namespace directional{

//arranging a line set on a triangle
//triangle is represented by a 3x2 matrix of (CCW) coordinates
//lines are Nx4 matrices of (origin, direction).
//line data is an integer associated with data on the line that gets inherited to the halfedges
//output is the DCEL of the result
//Outer face is deleted in post-process
void NFunctionMesher::arrange_on_triangle(const std::vector<EVector2>& triangle,
                                          const std::vector<std::pair<int, bool>>& triangleData,
                                          const std::vector<LinePencil>& linePencils,
                                          const std::vector<int>& linePencilData,
                                          std::vector<EVector2>& V,
                                          FunctionDCEL& triDcel){
  
  using namespace std;
  using namespace eigen;
  V = triangle;
  
  std::vector<HEData> inData;  //the lines that are inside
  std::vector<Segment2> inSegments;  //parameters of the line segments inside the triangle
  
  //std::cout<<"Triangle coordinates: "<<std::endl;
  /*for (int i=0;i<triangle.size();i++)
   std::cout<<triangle[i][0].get_d()<<" "<<triangle[i][1].get_d()<<std::endl;*/
  for (int i = 0; i < linePencils.size(); i++) {
    //checking for intersections with the triangle
    std::vector<ENumber> inParams, outParams;  //the parameters of intersection
    std::vector<int> inVertices, outVertices, inEdges, outEdges;
    linepencil_triangle_intersection(linePencils[i], triangle, inVertices, inEdges, outVertices, outEdges, inParams, outParams);
    //parameters for generating intersections with all other line pencils
    std::vector<Matrix<ENumber, 2, 2>> I2dts(linePencils.size()-i-1);
    std::Vector<Matrix<Enumber, 2, 1>> t00s(linePencils.size()-i-1);
    ENumber iso1Overlap;  //will not be used since line pencils of parametric lines are not supposed to overlap
    for (int j=i;j<linePencils.size();j++)
      linepencil_intersection(linePencils[i], linePencils[j], t00s[j-i]. I2dt[j-1],iso1Overlap);
    
    //creating the SV
    
    //generating
      
    /*cout<<"inParams: ";
     for (int i=0;i<inParams.size();i++)
     cout<<inParams[i].to_double()<<",";
     cout<<endl;
     cout<<"outParams: ";
     for (int i=0;i<outParams.size();i++)
     cout<<outParams[i].to_double()<<",";
     cout<<endl;
     cout<<"intEdges: ";
     for (int i=0;i<intEdges.size();i++)
     cout<<intEdges[i]<<",";
     cout<<endl;
     cout<<"intFaces: ";
     for (int i=0;i<intFaces.size();i++)
     cout<<intFaces[i]<<",";*/
    //cout<<endl;
    //checking cases of intersection
    for (int j=0;j<linePencils[i].numLines;j++){
      if (!intEdges[j] && !intFaces[j])
        continue;   //no (non-measure-zero) intersection
      
      HEData newData; newData.isFunction=true;
      newData.origNFunctionIndex = i;
      newData.origNFunctionIndex=linePencilData[i];
      inData.push_back(newData);
      EVector2 segSource = linePencils[i].p0 + linePencils[i].pVec*ENumber(j,1) + linePencils[i].direction * inParams[j];
      EVector2 segTarget = linePencils[i].p0 + linePencils[i].pVec*ENumber(j,1) + linePencils[i].direction * outParams[j];
      inSegments.push_back(Segment2(segSource, segTarget));
    }
    
    //std::cout<<"inParam: "<<inParam.get_d()<<std::endl;
    //std::cout<<"outParam: "<<outParam.get_d()<<std::endl;
    //std::cout<<"intersecting segment: ("<<segSource[0].get_d()<<","<<segSource[1].get_d()<<")->("<<segTarget[0].get_d()<<","<<segTarget[1].get_d()<<")"<<std::endl;
  }
  
  //pushing in triangle segments
  for (int i = 0; i < 3; i++) {
    HEData newData;
    newData.isFunction=false;
    newData.origHalfedge=triangleData[i].first;
    newData.origNFunctionIndex = -1;
    inData.push_back(newData);  //no data
    inSegments.push_back(Segment2(triangle[i], triangle[(i + 1) % 3]));
  }
  
  segment_arrangement(inSegments, inData, V, triDcel);
}



void NFunctionMesher::segment_arrangement(const std::vector<Segment2>& segments,
                                          const std::vector<HEData>& data,
                                          std::vector<EVector2>& V,
                                          FunctionDCEL& triDcel) {
  
  //First creating a graph of segment intersection
  
  //std::cout<<"************Computing segment arrangement*************"<<std::endl;
  //std::cout<<"list of segments: "<<std::endl;
  //for (int i=0;i<segments.size();i++)
  //    std::cout<<segments[i]<<" isFunction: "<<data[i].isFunction<<" origHalfedge:  "<<data[i].origHalfedge<<" function index: "<<data[i].origNFunctionIndex<<std::endl;
  //Creating arrangement vertices
  std::vector<EVector2> arrVertices;
  //std::vector<std::vector<int>> VS;  //list of participating segments per vertex
  std::vector<std::set<std::pair<ENumber, int>>> SV(segments.size());  //set of coordinates of intersection per segment
  for (int i=0;i<segments.size();i++) {
    for (int j = i + 1; j < segments.size(); j++) {
      if ((data[i].origNFunctionIndex==data[j].origNFunctionIndex)&&(data[i].isFunction)&&(data[j].isFunction))
        continue;   //don't try to intersect from the same line pencil
      std::vector<std::pair<ENumber, ENumber>> result = segment_segment_intersection(segments[i], segments[j]);
      
      if (result.empty())  //no intersection
        continue;  //that means the segments intersect away from the triangle.
      
      for (int r=0;r<result.size();r++){
        arrVertices.push_back(segments[i].source * (ENumber(1) - result[r].first) + segments[i].target * result[r].first);
        //std::cout<<"Segments ("<<i<<","<<j<<") create arragment vertex at "<<arrVertices[arrVertices.size()-1]<<std::endl;
        /*VS.push_back(std::vector<int>());
         VS[arrVertices.size() - 1].push_back(i);
         VS[arrVertices.size() - 1].push_back(j);*/
        SV[i].insert(std::pair<ENumber, int>(result[r].first, arrVertices.size() - 1));
        SV[j].insert(std::pair<ENumber, int>(result[r].second, arrVertices.size() - 1));
        
      }
      /*if (result.size) {  //pointwise intersection
       // TODO: figure out what happens if more than two lines at the same spot
       arrVertices.push_back(segments[i].first * (1 - t1) + segments[i].second * t1);
       VS[arrVertices.size() - 1].push_back(i);
       VS[arrVertices.size() - 1].push_back(j);
       SV[i].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 1));
       SV[j].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 1));
       }
       
       //Should make this aware of the double and put both into the edge data
       if (result==2) {  //subsegment; now entering two vertices, and letting the edges be entered later
       arrVertices.push_back(segments[i].first * (ENumber(1) - t1) + segments[i].second * t1);
       arrVertices.push_back(segments[j].first * (ENumber(1) - t2) + segments[j].second * t2);
       VS[arrVertices.size() - 2].push_back(i);
       VS[arrVertices.size() - 1].push_back(j);
       VS[arrVertices.size() - 2].push_back(i);
       VS[arrVertices.size() - 1].push_back(j);
       SV[i].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 2));
       SV[j].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 1));
       SV[i].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 2));
       SV[j].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 1));
       }*/
    }
  }
  
  //Creating the arrangement edges
  std::vector<std::pair<int, int>> arrEdges;
  std::vector<std::vector<HEData>> edgeData;
  for (int i=0;i<SV.size();i++){
    for (std::set<std::pair<ENumber, int>>::iterator si = SV[i].begin(); si!=SV[i].end();si++){
      std::set<std::pair<ENumber, int>>::iterator nextsi = si;
      nextsi++;
      if (nextsi!=SV[i].end()) {
        arrEdges.push_back(std::pair<int, int>(si->second, nextsi->second));
        //std::cout<<"Creating an edge ("<<si->second<<", "<<nextsi->second<<")"<<std::endl;
        std::vector<HEData> newEdgeData(1); newEdgeData[0]=data[i];
        edgeData.push_back(newEdgeData);
      }
    }
  }
  
  //unifying vertices with the same coordinates (necessary because some segments may intersect at the same point and segment overlaps
  auto VertexCompare = [](const std::pair<EVector2, int> a, const std::pair<EVector2, int> b) { return a.first < b.first; };
  std::set<std::pair<EVector2, int>, std::function<bool(const std::pair<EVector2, int>, const std::pair<EVector2, int>)>> uniqueVertices(VertexCompare);
  std::vector<int> uniqueVertexMap(arrVertices.size());
  std::vector<EVector2> uniqueArrVertices;
  int uniqueCounter=0;
  for (int i=0;i<arrVertices.size();i++){
    std::pair<EVector2, int> searchElement(arrVertices[i],-1);
    std::set<std::pair<EVector2, int>, decltype(VertexCompare)>::iterator si = uniqueVertices.find(searchElement);
    if (si==uniqueVertices.end()){
      uniqueVertexMap[i]=uniqueCounter;
      std::pair<EVector2, int> newElement = std::pair<EVector2, int>(arrVertices[i], uniqueCounter++);
      uniqueVertices.insert(newElement);
      uniqueArrVertices.push_back(arrVertices[i]);
    } else {
      uniqueVertexMap[i] = si->second;
    }
  }
  
  arrVertices = uniqueArrVertices;
  V = arrVertices;
  for (int i=0;i<arrEdges.size();i++)
    arrEdges[i]=std::pair<int, int>(uniqueVertexMap[arrEdges[i].first], uniqueVertexMap[arrEdges[i].second]);
  
  //std::cout<<"Edges after unifying vertices "<<std::endl;
  
  //unifying edges with the same vertices (aggregating data) or degenerated
  Eigen::VectorXi isDeadEdge=Eigen::VectorXi::Constant(arrEdges.size(),0);
  for (int i=0;i<arrEdges.size();i++) {
    //std::cout<<"("<<arrEdges[i].first<<", "<<arrEdges[i].second<<")"<<std::endl;
    if (arrEdges[i].first==arrEdges[i].second)
      isDeadEdge[i]=1;
    for (int j = i + 1; j < arrEdges.size(); j++) {
      if (((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second)) ||
          ((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second))) {
        isDeadEdge(j) = 1;
        edgeData[i].insert(edgeData[i].end(), edgeData[j].begin(), edgeData[j].end());
      }
    }
  }
  //cleaning dead edges
  std::vector<std::pair<int, int>> newArrEdges;
  std::vector<std::vector<HEData>> newEdgeData;
  for (int i=0;i<arrEdges.size();i++){
    if (isDeadEdge[i])
      continue;
    newArrEdges.push_back(arrEdges[i]);
    newEdgeData.push_back(edgeData[i]);
  }
  arrEdges=newArrEdges;
  edgeData=newEdgeData;
  
  //Generating the DCEL
  triDcel.vertices.resize(arrVertices.size());
  triDcel.edges.resize(arrEdges.size());
  triDcel.halfedges.resize(2*arrEdges.size());
  
  for (int i=0;i<arrVertices.size();i++){
    triDcel.vertices[i].ID = i;
  }
  
  for (int i=0;i<arrEdges.size();i++) {
    triDcel.edges[i].ID = i;
    
    triDcel.halfedges[2*i].ID=2*i;
    triDcel.halfedges[2*i+1].ID=2*i+1;
    
    //Consolidating the edge data
    triDcel.halfedges[2*i].data.isFunction = false;
    triDcel.halfedges[2*i+1].data.isFunction = false;
    for (int j=0;j<edgeData[i].size();j++) {
      if (edgeData[i][j].isFunction) {
        triDcel.halfedges[2*i].data.isFunction =  triDcel.halfedges[2*i+1].data.isFunction = true;
        triDcel.halfedges[2*i].data.origNFunctionIndex =  triDcel.halfedges[2*i+1].data.origNFunctionIndex = edgeData[i][j].origNFunctionIndex;
      }
      if (edgeData[i][j].origHalfedge>=0)
        triDcel.halfedges[2*i].data.origHalfedge =  triDcel.halfedges[2*i+1].data.origHalfedge = edgeData[i][j].origHalfedge;
    }
    
    triDcel.edges[i].halfedge = 2 * i;
    triDcel.halfedges[2*i].vertex=arrEdges[i].first;
    triDcel.halfedges[2*i+1].vertex=arrEdges[i].second;
    triDcel.vertices[arrEdges[i].first].halfedge=2*i;
    triDcel.vertices[arrEdges[i].second].halfedge=2*i+1;
    triDcel.halfedges[2*i].edge=triDcel.halfedges[2*i+1].edge = i;
    triDcel.halfedges[2*i].twin=2*i+1;
    triDcel.halfedges[2*i+1].twin=2*i;
    
  }
  
  //std::cout<<"Clean vertices and edges: "<<std::endl;
  //for (int i=0;i<arrVertices.size();i++)
  //    std::cout<<"arrVertex "<<i<<": "<<arrVertices[i]<<std::endl;
  
  /*for (int i=0;i<arrEdges.size();i++) {
   std::cout << "arrEdge " << i << ": " << arrEdges[i].first << "," << arrEdges[i].second << std::endl;
   std::cout <<" edge data" <<std::endl;
   for (int j=0;j<edgeData[i].size();j++)
   std::cout<<"isFunction: "<<edgeData[i][j].isFunction<<" origHalfedge: "<<edgeData[i][j].origHalfedge<<" function index :"<<edgeData[i][j].origNFunctionIndex<<std::endl;
   }*/
  
  /*for (int i=0;i<triDcel.halfedges.size();i++) {
   std::cout << "halfedge " << i << " of edge " << triDcel.halfedges[i].edge << " isFunction: " <<
   triDcel.halfedges[i].data.isFunction << " origHalfEdge :" <<
   triDcel.halfedges[i].data.origHalfedge << " origNFunctionIndex: " <<
   triDcel.halfedges[i].data.origNFunctionIndex << std::endl;
   }*/
  
  //Orienting segments around each vertex by CCW order
  for (int i=0;i<arrVertices.size();i++) {
    std::vector<std::pair<int,bool>> adjArrEdges;  //second is direction
    for (int j=0;j<arrEdges.size();j++) {
      if (arrEdges[j].first==i)
        adjArrEdges.push_back(std::pair<int,bool>(j,true));
      if (arrEdges[j].second==i)
        adjArrEdges.push_back(std::pair<int,bool>(j,false));
    }//not very efficient but probably not that bad
    
    /*std::cout<<"Orienting vertex "<<i<<std::endl;
     for (int k=0;k<adjArrEdges.size();k++)
     std::cout<<"Adjacent edge "<<adjArrEdges[k].first<<" with vertices "<<arrEdges[adjArrEdges[k].first].first<<","<<arrEdges[adjArrEdges[k].first].second<<std::endl;*/
    
    std::set<std::pair<ENumber, int>> CCWSegments;
    //using this slope function: https://math.stackexchange.com/questions/1450498/rational-ordering-of-vectors
    for (int j=0;j<adjArrEdges.size();j++) {
      EVector2 edgeVec = arrVertices[arrEdges[adjArrEdges[j].first].second] - arrVertices[arrEdges[adjArrEdges[j].first].first];
      edgeVec = (adjArrEdges[j].second ? edgeVec : -edgeVec);
      ENumber slopeFunc = slope_function(edgeVec);
      //std::cout<<"slope: "<<slopeFunc.get_d()<<" for edgeVec "<<edgeVec<<std::endl;
      CCWSegments.insert(std::pair<ENumber, int>(slopeFunc, j));
    }
    //std::cout<<"Ordering of edges"<<std::endl;
    /*for (std::set<std::pair<ENumber, int>>::iterator si = CCWSegments.begin(); si!=CCWSegments.end();si++)
     std::cout<<si->second<<",";
     std::cout<<std::endl;*/
    
    int currHE = -1;
    for (std::set<std::pair<ENumber, int>>::iterator si = CCWSegments.begin(); si!=CCWSegments.end();si++) {
      bool outgoing = adjArrEdges[si->second].second;
      int outCurrHE = (outgoing ? triDcel.edges[adjArrEdges[si->second].first].halfedge  : triDcel.halfedges[triDcel.edges[adjArrEdges[si->second].first].halfedge].twin);
      std::set<std::pair<ENumber, int>>::iterator nextsi = si; nextsi++;
      if (nextsi==CCWSegments.end())
        nextsi = CCWSegments.begin();
      
      outgoing = adjArrEdges[nextsi->second].second;
      int outNextHE = (outgoing ? triDcel.edges[adjArrEdges[nextsi->second].first].halfedge  : triDcel.halfedges[triDcel.edges[adjArrEdges[nextsi->second].first].halfedge].twin);
      triDcel.halfedges[outCurrHE].prev = triDcel.halfedges[outNextHE].twin;
      triDcel.halfedges[triDcel.halfedges[outNextHE].twin].next = outCurrHE;
      
      //triDcel.halfedges[triDcel.halfedges[outCurrHE].twin].next=outNextHE;
      //triDcel.halfedges[outNextHE].prev = triDcel.halfedges[outCurrHE].twin;
    }
  }
  
  //generating faces (at this stage, there is also an outer face)
  int currFace=0;
  for (int i=0;i<triDcel.halfedges.size();i++) {
    if (triDcel.halfedges[i].face != -1)
      continue;  //already been assigned
    
    
    FunctionDCEL::Face newFace;
    newFace.ID = currFace++;
    //std::cout<<"New face "<<newFace.ID<<std::endl;
    //std::cout<<"edges: ";
    int beginHE = i;
    newFace.halfedge=beginHE;
    int currHE = beginHE;
    int counter=0;
    do {
      //std::cout<<currHE<<",";
      triDcel.halfedges[currHE].face = newFace.ID;
      currHE=triDcel.halfedges[currHE].next;
      counter++;
      assert ("something wrong with the face" && counter<10000);
    }while (currHE!=beginHE);
    //std::cout<<std::endl;
    triDcel.faces.push_back(newFace);
  }
  int numFaces=currFace;
  
  //EXPENSIVE! comment out after being sure
  assert("DCEL not consistent!" && triDcel.check_consistency(false, true, true, false));
  
  //Removing the outer face and deleting all associated halfedges
  //identifying it by the only polygon with negative signed area (expensive?)
  int outerFace=-1;
  for (int f = 0;f<numFaces;f++){
    std::vector<EVector2> faceVectors;
    int beginHE = triDcel.faces[f].halfedge;
    int currHE = beginHE;
    do{
      faceVectors.push_back(V[triDcel.halfedges[triDcel.halfedges[currHE].next].vertex] - V[triDcel.halfedges[currHE].vertex]);
      //std::cout<<"face vector :"<<faceVectors[faceVectors.size()-1]<<std::endl;
      currHE= triDcel.halfedges[currHE].next;
    }while(currHE!=beginHE);
    ENumber sfa = signed_face_area(faceVectors);
    //std::cout<<"Signed area of face "<<f<<": "<<sfa.get_d()<<std::endl;
    if (sfa<ENumber(0)){
      outerFace=f;
      break;
    }
  }
  assert("Didn't find outer face!" && outerFace!=-1);
  
  //invalidating outer face
  triDcel.faces[outerFace].valid=false;
  for (int i=0;i< triDcel.halfedges.size();i++) {
    if (triDcel.halfedges[i].face != outerFace)
      continue;
    
    triDcel.halfedges[i].valid=false;
    triDcel.halfedges[triDcel.halfedges[i].twin].twin = -1;
    triDcel.edges[triDcel.halfedges[i].edge].halfedge=triDcel.halfedges[i].twin;
    triDcel.vertices[triDcel.halfedges[i].vertex].halfedge = triDcel.halfedges[triDcel.halfedges[i].twin].next;
    //dcel.VH( dcel.HV( dcel.nextH(i)))= dcel.twinH(i);
    //dcel.VH( dcel.HV(i))= dcel.nextH( dcel.twinH(i));
  }
  
  //removing dead edges
  triDcel.clean_mesh();
  //EXPENSIVE! comment out after being sure
  assert("DCEL not consistent!" && triDcel.check_consistency(false, true, true, false));
}


void NFunctionMesher::generate_mesh(const unsigned long resolution = 1e7) {
  
  using namespace std;
  using namespace Eigen;
  
  //Looping over all triangles, and arranging parametric lines for each
  auto start = std::chrono::high_resolution_clock::now();
  for (int findex = 0; findex < origMesh.F.rows(); findex++) {
    
    //building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face
    vector<ENumber> triExactNFunction = exactNFunction[findex];
    if (findex % 100 == 0){
      auto end = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
      std::cout<<"Triangle "<<findex<<" exactNFunction:"<<std::endl;
      std::cout << "Execution time: " << duration.count()/1e+6 << " seconds" << std::endl;
      start = std::chrono::high_resolution_clock::now();
      
    }
    
    vector<ENumber> minFuncs(mfiData.N);
    vector<ENumber> maxFuncs(mfiData.N);
    for (int k = 0; k < mfiData.N; k++) {
      minFuncs[k] = ENumber(327600);
      maxFuncs[k] = ENumber(-327600);
    }
    
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < mfiData.N; j++) {
        if (triExactNFunction[mfiData.N * i +  j] > maxFuncs[j])
          maxFuncs[j] = triExactNFunction[mfiData.N * i + j];
        if (triExactNFunction[mfiData.N * i + j] < minFuncs[j])
          minFuncs[j] = triExactNFunction[mfiData.N * i + j];
      }
    }
    
    //building the one-triangle arrangement
    std::vector<EVector2> ETriPoints2D(3);
    std::vector<EVector3> ETriPoints3D(3);
    ETriPoints2D[0][0] = ENumber(0);
    ETriPoints2D[0][1] = ENumber(0);
    ETriPoints2D[1][0] = ENumber(1);
    ETriPoints2D[1][1] = ENumber(0);
    ETriPoints2D[2][0] = ENumber(0);
    ETriPoints2D[2][1] = ENumber(1);
    
    for (int i = 0; i < 3; i++) {
      RowVector3d position = origMesh.V.row(origMesh.F(findex, i));
      double tol = 1.0/(double)resolution;
      ENumber x = ENumber(position(0), tol);
      ENumber y = ENumber(position(1), tol);
      ENumber z = ENumber(position(2), tol);
      
      EVector3 xyz;
      xyz[0] = x;
      xyz[1] = y;
      xyz[2] = z;
      ETriPoints3D[i]=xyz;
    }
    
    
    int DomEdge;
    std::vector<std::pair<int, bool>> triangleData;  triangleData.resize(3); //of the triangles, to be put into arrangement data
    int eiterate = origMesh.dcel.faces[findex].halfedge;
    for (int i=0;i<3;i++){
      if ((origMesh.dcel.halfedges[eiterate].twin < 0) || (origMesh.dcel.halfedges[eiterate].twin > eiterate))
        triangleData[i].first  = eiterate;
      else
        triangleData[i].first = origMesh.dcel.halfedges[eiterate].twin;
      
      triangleData[i].second = (origMesh.dcel.halfedges[eiterate].twin < 0);  //is it a boundary
      
      eiterate = origMesh.dcel.halfedges[eiterate].next;
    }
    
    //Generating the parametric lines for the canonical triangle
    vector<LinePencil> linePencils;
    vector<int> linePencilData;
    vector<EVector2> isoDirections(mfiData.N);
    for (int funcIter = 0; funcIter < mfiData.N; funcIter++) {
      
      vector<EInt> isoValues;
      EInt q, r;
      div_mod(minFuncs[funcIter].num, minFuncs[funcIter].den, q, r);
      EInt minIsoValue = q + (minFuncs[funcIter].num < 0 ? -1 : 0);
      div_mod(maxFuncs[funcIter].num, maxFuncs[funcIter].den, q, r);
      EInt maxIsoValue = q + (maxFuncs[funcIter].num < 0 ? 0 : 1);
      for (EInt isoValue = minIsoValue; isoValue <= maxIsoValue; isoValue=isoValue+EInt(1))
        isoValues.push_back(isoValue);
      
      /*std::cout<<"minFuncs[funcIter]: "<<minFuncs[funcIter].to_double()<<std::endl;
       std::cout<<"maxFuncs[funcIter]: "<<maxFuncs[funcIter].to_double()<<std::endl;
       std::cout<<"minIsoValue: "<<minIsoValue.to_string()<<std::endl;
       std::cout<<"maxIsoValue: "<<maxIsoValue.to_string()<<std::endl;*/
      
      //computing gradient of function in plane
      EVector2 e01 = ETriPoints2D[1] - ETriPoints2D[0];
      EVector2 e12 = ETriPoints2D[2] - ETriPoints2D[1];
      EVector2 e20 = ETriPoints2D[0] - ETriPoints2D[2];
      
      //a and b values of lines
      EVector2 gradVector = triExactNFunction[2 * mfiData.N + funcIter] * EVector2({-e01[1], e01[0]}) +
      triExactNFunction[0 * mfiData.N + funcIter] * EVector2({-e12[1], e12[0]}) +
      triExactNFunction[1 * mfiData.N + funcIter] * EVector2({-e20[1], e20[0]});
      
      isoDirections[funcIter] = gradVector;
      
      //Number avgFuncValue = (funcValues[0](funcIter)+funcValues[1](funcIter)+funcValues[2](funcIter))/3.0;
      //TODO: find c = z1*u+z2 of ax+by+c(u) ad then use it to generate all values between floor and ceil.
      
      //pinv of [a 1;b 1;c 1] is [           2*a - b - c,           2*b - a - c,           2*c - b - a]
      //[ b^2 - a*b + c^2 - a*c, a^2 - b*a + c^2 - b*c, a^2 - c*a + b^2 - c*b]/(2*a^2 - 2*a*b - 2*a*c + 2*b^2 - 2*b*c + 2*c^2)
      
      ENumber a = triExactNFunction[0 * mfiData.N + funcIter];
      ENumber b = triExactNFunction[1 * mfiData.N + funcIter];
      ENumber c = triExactNFunction[2 * mfiData.N + funcIter];
      if ((a == b) && (b == c))
        continue;  //that means a degenerate function on the triangle
      
      //cout<<"a,b,c: "<<a.to_double()<<","<<b.to_double()<<","<<c.to_double()<<endl;
      
      ENumber rhs[3];
      rhs[0] = -gradVector[0] * ETriPoints2D[0][0] - gradVector[1] * ETriPoints2D[0][1];
      rhs[1] = -gradVector[0] * ETriPoints2D[1][0] - gradVector[1] * ETriPoints2D[1][1];
      rhs[2] = -gradVector[0] * ETriPoints2D[2][0] - gradVector[1] * ETriPoints2D[2][1];
      
      ENumber invM[2][3];
      invM[0][0] = ENumber(2) * a - b - c;
      invM[0][1] = ENumber(2) * b - a - c;
      invM[0][2] = ENumber(2) * c - b - a;
      invM[1][0] = b * b - a * b + c * c - a * c;
      invM[1][1] = a * a - b * a + c * c - b * c;
      invM[1][2] = a * a - c * a + b * b - c * b;
      for (int row = 0; row < 2; row++)
        for (int col = 0; col < 3; col++)
          invM[row][col] /= (ENumber(2) * (a * a - a * b - a * c + b * b - b * c + c * c));
      
      ENumber x[2];
      x[0] = invM[0][0] * rhs[0] + invM[0][1] * rhs[1] + invM[0][2] * rhs[2];
      x[1] = invM[1][0] * rhs[0] + invM[1][1] * rhs[1] + invM[1][2] * rhs[2];
      
      //sanity check
      /*ENumber error[3];
       error[0] = x[0] * a + x[1] - rhs[0];
       error[1] = x[0] * b + x[1] - rhs[1];
       error[2] = x[0] * c + x[1] - rhs[2];*/
      
      //generating all lines
      LinePencil lp;
      lp.direction[0] = -gradVector[1];
      lp.direction[1] = gradVector[0];
      //lp.pVec = gradVector;  //that means not orthogonal!
      if (gradVector[1] != ENumber(0)) {
        lp.p0[0] = ENumber(0);
        lp.p0[1] = -(x[0]* isoValues[0] + x[1]) / gradVector[1];
        lp.pVec[0] = ENumber(0);
        lp.pVec[1] = -x[0]/gradVector[1];
      } else {
        lp.p0[1] = ENumber(0);
        lp.p0[0] = -(x[0]* isoValues[0] + x[1]) / gradVector[0];
        lp.pVec[1] = ENumber(0);
        lp.pVec[0] = -x[0]/gradVector[0];
      }
      lp.numLines = isoValues.size();
      //lp.minIsoValue = isoValues[0];
      //lp.maxIsoValue = isoValues[isoValues.size()-1];
      linePencils.push_back(lp);
      linePencilData.push_back(funcIter);
      /*for (int isoIndex = 0; isoIndex < isoValues.size(); isoIndex++) {
       ENumber currc = x[0]* isoValues[isoIndex] + x[1];
       EVector2 lineVector;
       
       EVector2 linePoint;
       if (gradVector[1] != ENumber(0)) {
       linePoint[0] = ENumber(0);
       linePoint[1] = -currc / gradVector[1];
       } else {
       linePoint[1] = ENumber(0);
       linePoint[0] = -currc / gradVector[0];
       }
       
       paramLines.push_back(Line2(linePoint, lineVector));
       lineData.push_back(funcIter);
       }*/
      
    }
    
    FunctionDCEL localArrDcel;
    vector<EVector2> localV;
    arrange_on_triangle(ETriPoints2D, triangleData, linePencils, linePencilData, localV, localArrDcel);
    
    //vector<EVector3> ELocalV3D(localV.size());
    // MatrixXd localV3D(localV.size(), 3);
    //converting the vertices to 3D
    for (int i = 0; i < localV.size(); i++) {
      //checking if this is an original vertex
      //std::cout<<"Converting vertex "<<i<<" to 3D "<<std::endl;
      bool isOrigTriangle = false;
      for (int j = 0; j < 3; j++) {
        if (localV[i] == ETriPoints2D[j]) {
          localArrDcel.vertices[i].data.eCoords = ETriPoints3D[j];
          isOrigTriangle = true;
        }
      }
      
      /*if (isOrigTriangle)
       continue;    //this is probably not needed but covered by barycentric coordinates but w/e*/
      
      //finding out barycentric coordinates
      ENumber baryValues[3];
      ENumber sum(0);
      for (int j = 0; j < 3; j++) {
        //ETriangle2D t(vi->point(), ETriPoints2D[(i + 1) % 3], ETriPoints2D[(i + 2) % 3]);
        vector<EVector2> inTri(3);
        inTri[0] = localV[i];
        inTri[1] = ETriPoints2D[(j + 1) % 3];
        inTri[2] = ETriPoints2D[(j + 2) % 3];
        
        baryValues[j] = triangle_area(inTri);
        sum += baryValues[j];
      }
      for (int j = 0; j < 3; j++)
        baryValues[j] /= sum;
      
      //std::cout<<"baryValues vertex "<< baryValues[0].get_d()<<","<< baryValues[1].get_d()<<","<< baryValues[2].get_d()<<std::endl;
      
      localArrDcel.vertices[i].data.eCoords = EVector3({ENumber(0), ENumber(0), ENumber(0)});
      for (int j = 0; j < 3; j++) {
        //std::cout<<"ETriPoints3D["<<j<<"]: "<<ETriPoints3D[j]<<std::endl;
        localArrDcel.vertices[i].data.eCoords =
        localArrDcel.vertices[i].data.eCoords + ETriPoints3D[j] * baryValues[j];
      }
      //localArrDcel.vertices[i].data.eCoords = EVector3({localV[i][0], localV[i][1], 0});
      
      //std::cout<<"localArrDcel.vertices[i].data.eCoords: "<<localArrDcel.vertices[i].data.eCoords<<std::endl;
      
      localArrDcel.vertices[i].data.coords << localArrDcel.vertices[i].data.eCoords[0].to_double(),
      localArrDcel.vertices[i].data.eCoords[1].to_double(),
      localArrDcel.vertices[i].data.eCoords[2].to_double();
      
    }
    
    //aggregating to the general DCEL
    genDcel.aggregate_dcel(localArrDcel);
    //localArrDcel.check_consistency(true, false, false, false);
    //genDcel.check_consistency(true, false, false, false);
  }
  genDcel.check_consistency(true, false, false, false);
}

} //namespace directional


#endif
