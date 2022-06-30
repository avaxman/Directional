// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_GLYPH_LINES_MESH_H
#define DIRECTIONAL_GLYPH_LINES_MESH_H

#include <igl/igl_inline.h>
#include <igl/colon.h>
#include <igl/speye.h>
#include <directional/angled_arrows.h>
#include <Eigen/Core>


namespace directional
{
  
  
  // Creates mesh elements that comprise glyph drawing of a directional field.
  // Inputs:
  //  V:          #V X 3 vertex coordinates.
  //  F:          #F by 3 face vertex indices.
  //  EF:         #E by 2 edge-face adjacency matrix (only used when sparsity !=0)
  //  rawField:   A directional field in raw xyzxyz form
  //  glyphColor: An array of either 1 by 3 color values for each vector, #F by 3 colors for each individual directional or #F*N by 3 colours for each individual vector, ordered by #F times vector 1, followed by #F times vector 2 etc.
  //  length, width,  height: of the glyphs depicting the directionals
  //  N:        The degree of the field.
  
  // Outputs:
  //  fieldV: The vertices of the field mesh
  //  fieldF: The faces of the field mesh
  //  fieldC: The colors of the field mesh
  
  void IGL_INLINE glyph_lines_mesh(const Eigen::MatrixXd& sources,
                                   const Eigen::MatrixXd& normals,
                                   const Eigen::MatrixXi& adjSpaces,
                                   const Eigen::MatrixXd& extField,
                                   const Eigen::MatrixXd& glyphColor,
                                   const double length,
                                   const double width,
                                   const double height,
                                   const int sparsity,
                                   Eigen::MatrixXd &fieldV,
                                   Eigen::MatrixXi &fieldF,
                                   Eigen::MatrixXd &fieldC)
  {
    using namespace Eigen;
    using namespace std;
    
    int N=extField.cols()/3;
    
    double angle = 2*igl::PI/(double)(N);
    if (N==1) angle=igl::PI;
    Eigen::MatrixXd vectorColors, P1, P2;
    
    VectorXi sampledSpaces;
    if (sparsity!=0){
      //creating adjacency matrix
      vector<Triplet<int>> adjTris;
      for (int i=0;i<adjSpaces.rows();i++)
        if ((adjSpaces(i,0)!=-1)&&(adjSpaces(i,1)!=-1)){
          adjTris.push_back(Triplet<int>(adjSpaces(i,0), adjSpaces(i,1),1));
          adjTris.push_back(Triplet<int>(adjSpaces(i,1), adjSpaces(i,0),1));
        }
      
      SparseMatrix<int> adjMat(extField.rows(),extField.rows());
      adjMat.setFromTriplets(adjTris.begin(), adjTris.end());
      SparseMatrix<int> newAdjMat(extField.rows(),extField.rows()),matMult;
      igl::speye(extField.rows(), extField.rows(), matMult);
      for (int i=0;i<sparsity;i++){
        matMult=matMult*adjMat;
        newAdjMat+=matMult;
      }
      
      //cout<<"newAdjMat: "<<newAdjMat<<endl;
      
      adjMat=newAdjMat;
      
      vector<set<int>> ringAdjacencies(extField.rows());
      for (int k=0; k<adjMat.outerSize(); ++k){
        for (SparseMatrix<int>::InnerIterator it(adjMat,k); it; ++it){
          ringAdjacencies[it.row()].insert(it.col());
          ringAdjacencies[it.col()].insert(it.row());
        }
      }
      
      VectorXi sampleMask=VectorXi::Zero(extField.rows());
      for (int i=0;i<extField.rows();i++){
        if (sampleMask(i)!=0) //occupied face
          continue;
        
        sampleMask(i)=2;
        //clearing out all other faces
        for (set<int>::iterator si=ringAdjacencies[i].begin();si!=ringAdjacencies[i].end();si++){
          if (sampleMask(*si)==0)
            sampleMask(*si)=1;
         
        }
      }
      
      vector<int> samplesList;
      for (int i=0;i<sampleMask.size();i++)
        if (sampleMask(i)==2)
          samplesList.push_back(i);
      
      sampledSpaces = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(samplesList.data(), samplesList.size());
      
    } else igl::colon(0,1,extField.rows()-1,sampledSpaces);
    
    MatrixXd vectNormals(sampledSpaces.rows()*N,3);
    P1.resize(sampledSpaces.rows() * N, 3);
    P2.resize(sampledSpaces.rows() * N, 3);
    vectorColors.resize(sampledSpaces.rows() * N, 3);
    
    //normals.array() *= width;
    for (int i=0;i<sampledSpaces.size();i++)
      for (int j=0;j<N;j++){
        P1.row(j*sampledSpaces.size()+i) = sources.row(sampledSpaces(i));
        P2.row(j*sampledSpaces.size()+i) = extField.block(sampledSpaces(i),j*3,1,3);
        vectNormals.row(j*sampledSpaces.size()+i) = normals.row(sampledSpaces(i)).array()*width;
      }
    
    /*P1 = barycenters.replicate(N, 1);
    
    for (int i = 0; i < N; i++)
      P2.middleRows(F.rows()*i, F.rows()) = rawField.middleCols(3*i, 3);*/
    
    P2.array() *= length;
    P2 += P1;
    
    // Duplicate colors so each glyph gets the proper color
    if (glyphColor.rows() == 1)
      vectorColors = glyphColor.replicate(P1.rows(), 1);
    else if ((glyphColor.rows() == extField.rows())&&(glyphColor.cols()==3)){
      MatrixXd sampledColors(sampledSpaces.size(),3);
      for (int i=0;i<sampledSpaces.size();i++)
        sampledColors.row(i)=glyphColor.row(sampledSpaces(i));
      vectorColors = sampledColors.replicate(N, 1);
    }
    else{
      for (int i=0;i<N;i++)
        vectorColors.block(i*sampledSpaces.rows(),0,sampledSpaces.rows(),3)=glyphColor.block(0,3*i,sampledSpaces.rows(),3);
    }
    
    Eigen::MatrixXd Vc, Cc, Vs, Cs;
    Eigen::MatrixXi Fc, Fs;
    directional::angled_arrows(P1,P2,vectNormals, width/length, height, angle, vectorColors, fieldV, fieldF, fieldC);
    
    
  }

  
  //A version without specification of glyph dimensions
  void IGL_INLINE glyph_lines_mesh(const Eigen::MatrixXd& sources,
                                   const Eigen::MatrixXd& normals,
                                   const Eigen::MatrixXi& adjSpaces,
                                   const Eigen::MatrixXd& extField,
                                   const Eigen::MatrixXd &glyphColors,
                                   const double sizeRatio,
                                   const double avgScale,
                                   Eigen::MatrixXd &fieldV,
                                   Eigen::MatrixXi &fieldF,
                                   Eigen::MatrixXd &fieldC,
                                   
                                   const int sparsity=0,
                                   const double offsetRatio = 0.2)
  {
    //double l = igl::avg_edge_length(V, F);
    glyph_lines_mesh(sources, normals, adjSpaces, extField, glyphColors, sizeRatio*avgScale/3.0, sizeRatio*avgScale/15.0,  avgScale*offsetRatio, sparsity, fieldV, fieldF, fieldC);
  }
  
}

#endif
