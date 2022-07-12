// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_INTRINSIC_FACE_TANGENT_BUNDLE_H
#define DIRECTIONAL_INTRINSIC_FACE_TANGENT_BUNDLE_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/boundary_loop.h>
#include <igl/doublearea.h>
#include <directional/dual_cycles.h>
#include <directional/TriMesh.h>


namespace directional{

//This is the interface class for any directional fields represented in cartesian coordinates, of any order N.
class IntrinsicFaceTangentBundle::public TangentBundle{
public:
  
  TriMesh* mesh;
  Eigen::MatrixXd tangentStartAngles;  //where each edge begins on the intrinsic space
  
  virtual discTangTypeEnum discTangType() const {return discTangTypeEnum::VERTEX_SPACES;}
  
  IntrinsicVertexTangentBundle(){}
  ~IntrinsicVertexTangentBundle(){}
  
  void IGL_INLINE init(const TriMesh& _mesh){
    
    typedef std::complex<double> Complex;
    mesh = &_mesh;
    
    //adjacency relation is by dual edges.
    adjSpaces = mesh->EV;
    oneRing = mesh->VE;
    Eigen::VectorXi valence = mesh->vertexValence;
    sources = mesh->V;
    normals = mesh->vertexNormals;
    cycleSources = mesh->barycenters;
    cycleNormals = mesh->faceNormals;
    
    //creating vertex tangent lookup table
    tangentStartAngles.resize(mesh->V.rows(), mesh->vertexValence.maxCoeff());
    for (int i=0;i<mesh->V.rows();i++){
      //TODO: boundaries
      double angleSum = 2*igl::PI - mesh->GaussianCurvature(i);
      tangentStartAngles.col(0).setZero();  //the first angle
      Eigen::RowVector3d prevEdgeVector = mesh->V.row(mesh->HV(mesh->nextH(mesh->VH(i))))-mesh->V.row(i);
      int hebegin = mesh->VH(i);
      int heiterate = mesh->twinH(mesh->prevH(hebegin));
      int j=1;
      do{
        Eigen::RowVector3d currEdgeVector = mesh->V.row(mesh->HV(mesh->nextH(heiterate)))-mesh->V.row(i);
        double angleDiff = std::acos(currEdgeVector.dot(prevEdgeVector)/(prevEdgeVector.norm()*currEdgeVector.norm()));
        tangentStartAngles(i,j)=tangentStartAngles(i,j-1)+2*igl::PI*angleDiff/angleSum;
        heiterate = mesh->twinH(mesh->prevH(heiterate));
        j++;
        prevEdgeVector=currEdgeVector;
      }while (heiterate!=hebegin);
    }
    
    //connection is the ratio of the complex representation of mutual edges
    connection.resize(mesh->EV.rows(),1);  //the difference in the angle representation of edge i from EV(i,0) to EV(i,1)
    Eigen::MatrixXd edgeVectors(mesh->EV.rows(), 3);
    for (int i = 0; i < mesh->EV.rows(); i++) {
      //edgeVectors.row(i) = (mesh->V.row(mesh->EV(i, 1)) - mesh->V.row(mesh->EV(i, 0))).normalized();
      
      //looking up edge in each tangent space
      Complex ef,eg;
      int hebegin = mesh->VH(mesh->EV(i,0));
      int heiterate = hebegin;
      int j=0;
      do{
        if (mesh->HE(heiterate)==i)
          ef = exp(Complex(0,tangentStartAngles(mesh->EV(i, 0),j)));
        heiterate = mesh->twinH(mesh->prevH(heiterate));
        j++;
      }while (heiterate!=hebegin);
      
      hebegin = mesh->VH(mesh->EV(i,1));
      heiterate = hebegin;
      j=0;
      do{
        if (mesh->HE(heiterate)==i)
          eg = exp(Complex(0,tangentStartAngles(mesh->EV(i, 1),j)));
        heiterate = mesh->twinH(mesh->prevH(heiterate));
        j++;
      }while (heiterate!=hebegin);
      
      connection(i) = -eg / ef;
    }
    
    local2Cycle.resize(mesh->F.rows());
    cycleCurvatures=Eigen::VectorXd::Zero(mesh->F.rows());
    cycles.resize(mesh->F.rows(), mesh->EV.rows());  //TODO: higher genus and boundaries
    innerAdjacencies.resize(mesh->EV.rows());
    std::vector<Eigen::Triplet<double>> cyclesTriplets;
    for (int i=0;i<mesh->F.rows();i++){
      local2Cycle(i)=i;
      std::complex<double> complexHolonomy(1.0,0.0);
      for (int j=0;j<3;j++){
        dualCyclesTriplets.push_back(Eigen::Triplet<double>(i,mesh->FE(i,j),mesh->FEs(i,j)));
        if (mesh->FEs(i,j)>0)
          complexHolonomy*=connection(mesh->FE(i,j));
        else
          complexHolonomy/=connection(mesh->FE(i,j));
      }
      cycleCurvatures(i)=arg(complexHolonomy);
    }
    cycles.setFromTriplets(cyclesTriplets.begin(), cyclesTriplets.end());
    
    for (int i=0;i<mesh->EV.rows();i++) //TODO: boundaries
      innerAdjacencies(i)=i;
    //directional::dual_cycles(mesh->V, mesh->F, mesh->EV, mesh->EF, dualCycles, cycleCurvatures, element2Cycle, innerAdjacencies);
    
    //drawing from mesh geometry
    
    /************Smoothness matrices****************/
    stiffnessWeights=Eigen::VectorXd::Zero(mesh->EV.rows());
    
    //cotangent weights
    Eigen::MatrixXd faceCotWeights=Eigen::MatrixXd::Zero(mesh->F.rows(),3);
    for (int i=0;i<mesh->F.rows();i++){
      for (int j=0;j<3;j++){
        Eigen::RowVector3d vec12 = mesh->V.row(mesh->F(i,(j+1)%3))-mesh->V.row(mesh->F(i,j));
        Eigen::RowVector3d vec13 = mesh->V.row(mesh->F(i,(j+2)%3))-mesh->V.row(mesh->F(i,j));
        double cosAngle = vec12.dot(vec13);
        double sinAngle = (vec12.cross(vec13)).norm();
        if (std::abs(sinAngle)>10e-7) //otherwise defaulting to zero
          faceCotWeights(i,j) = cosAngle/sinAngle;
        
        stiffnessWeights(mesh->FE(i,j))+=0.5*faceCotWeights(i,j);
      }
    }
    
    //masses are vertex voronoi areas
    Eigen::MatrixXd dAreas;
    igl::doublearea(mesh->V,mesh->F,dAreas);
    massWeights = Eigen::VectorXd::Zero(mesh->V.rows());
    for (int i=0;i<mesh->F.rows();i++)
      for (int j=0;j<3;j++)
        massWeights(mesh->F(i,j)) += dAreas(i)/6.0;
  
  }
  
  
  //projecting an arbitrary set of extrinsic vectors (e.g. coming from user-prescribed constraints) into intrinsic vectors.
  Eigen::MatrixXd  virtual IGL_INLINE project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const{
    assert(tangentSpaces.rows()==extDirectionals.rows());
    Eigen::MatrixXd intDirectionals(tangentSpaces.rows(),2);
    
    for (int j=0;j<tangentSpaces.rows();j++){
      intDirectionals.row(j)<<(extDirectionals.row(j).array()*mesh->VBx.row(tangentSpaces(j)).array()).sum(),(extDirectionals.row(j).array()*mesh->VBy.row(tangentSpaces(j)).array()).sum();
    }
    return intDirectionals;
  }

  
  //projecting extrinsic to intrinsic
  Eigen::MatrixXd virtual IGL_INLINE project_to_extrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& intDirectionals) const {
    
    assert(tangentSpaces.rows()==intDirectionals.rows());
    Eigen::MatrixXd extDirectionals(tangentSpaces.rows(),3);
    
    extDirectionals.conservativeResize(intDirectionals.rows(),intDirectionals.cols()*3/2);
    for (int i=0;i<intDirectionals.rows();i++)
      for (int j=0;j<intDirectionals.cols();j+=2)
        extDirectionals.block(i,3*j/2,1,3)=mesh->VBx.row(i)*intDirectionals(i,j)+mesh->VBy.row(i)*intDirectionals(i,j+1);
    }
    return extDirectionals;
  }
  
};

}



#endif
  

