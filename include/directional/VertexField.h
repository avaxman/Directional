// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_VERTEXFIELD_H
#define DIRECTIONAL_VERTEXFIELD_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/boundary_loop.h>
#include <igl/doublearea.h>
#include <directional/dual_cycles.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>

namespace directional{

class VertexField: public CartesianField{
public:
    
  VertexField(){}
  VertexField(const TriMesh& _mesh):CartesianField(_mesh){}
  ~VertexField(){}
  
  void IGL_INLINE init_field(const TriMesh& _mesh, const int _fieldType, const int _N){
    
    typedef std::complex<double> Complex;
    mesh = &_mesh;
    fieldType = _fieldType;
    N = _N;
    
    //adjacency relation is by dual edges.
    adjSpaces = mesh->EV;
    oneRing = mesh->VE;
    valence = mesh->vertexValence;
    intField.conservativeResize(mesh->V.rows(),2*N);
    extField.conservativeResize(mesh->V.rows(),3*N);
    
    //creating vertex tangent lookup table
    tangentStartAngles.resize(oneRing.rows(), oneRing.cols());
    for (int i=0;i<mesh->V.rows();i++){
      //TODO: boundaries
      double angleSum = 2*igl::PI - mesh->GaussianCurvature;
      tangentStartAngles.col(0).setZero();  //the first angle
      RowVector3d prevEdgeVector = mesh.EV(mesh.VE(i,0),1)-mesh.EV(mesh.VE(i,0),0);
      prevEdgeVector = (mesh.EV(mesh.VE(i,0),0) == i ? prevEdgeVector : -prevEdgeVector);
      for (int j=1;j<valence(i);j++){
        RowVector3d currEdgeVector = mesh.EV(mesh.VE(i,j),1)-mesh.EV(mesh.VE(i,j),0);
        currEdgeVector = (mesh.EV(mesh.VE(i,j),0) == i ? currEdgeVector : -currEdgeVector);
        angleDiff = std::acos(currEdgeVector.dot(prevEdgeVector));
        tangentStartAngles(i,j)=tangentStartAngles(i,j-1)+2*igl::PI*angleDiff/angleSum;
      }
    }
    
    //connection is the ratio of the complex representation of edges
    connection.resize(mesh->EV.rows(),1);  //the difference in the angle representation of edge i from EV(i,0) to EV(i,1)
    Eigen::MatrixXd edgeVectors(mesh->EV.rows(), 3);
    for (int i = 0; i < mesh->EV.rows(); i++) {
      edgeVectors.row(i) = (mesh->V.row(mesh->EV(i, 1)) - mesh->V.row(mesh->EV(i, 0))).normalized();
      
      //looking up edge in each tangent space
      Complex ef,eg;
      for (int j=0;j<valence(mesh->EV(i, 0));j++)
        if (oneRing(mesh->EV(i, 0),j)==i)
          ef = exp(Complex(0,tangentStartAngles(i,j)));
      
      for (int j=0;j<valence(mesh->EV(i, 1));j++)
        if (oneRing(mesh->EV(i, 1),j)==i)
          ef = exp(Complex(0,tangentStartAngles(i,j)));
      
      connection(i) = -eg / ef;
    }
    
    //TODO: cycles, cycleCurvature
    //directional::dual_cycles(mesh->V, mesh->F, mesh->EV, mesh->EF, dualCycles, cycleCurvatures, element2Cycle, innerAdjacencies);
    
    //drawing from mesh geometry
    
    /************Smoothness matrices****************/
    stiffnessWeights=Eigen::VectorXd::Zero(mesh->EV.rows());
    
    //cotangent weights
    Eigen::MatrixXd faceCotWeights=Eigen::MatrixXd::Zero(mesh.F.rows(),3);
    for (int i=0;i<mesh.F.rows();i++){
      for (int j=0;j<3;j++){
        RowVector3d vec12 = mesh.V(mesh.F(i,(j+1)%3))-mesh.V(mesh.F(i,j));
        RowVector3d vec13 = mesh.V(mesh.F(i,(j+2)%3))-mesh.V(mesh.F(i,j));
        cosAngle = vec12.dot(vec13);
        sinAngle = sqrt(vec12.squaredNorm()*vec13*squaredNorm()-cosAngle);
        if (std::abs(sinAngle))>10e-7) //otherwise defaulting to zero
          faceCotWeights(i,j) = cosAngle/sinAngle;
        
        stiffnessWeights(FE(i,j))+=0.5*faceCotWeights(i,j);
      }
    }
    
    //masses are vertex voronoi areas
    Eigen::MatrixXd dAreas
    igl::doublearea(mesh->V,mesh->F,dAreas);
    massWeights = Eigen::VertexXd::Zero(mesh.V.rows());
    for (int i=0;i<mesh.F.rows();i++)
      for (int j=0;j<mesh.V.rows();j++)
        massWeights(F(i,j)) = dAreas(i)/6.0;
    
    
  }
  
  void IGL_INLINE set_extrinsic_field(const Eigen::MatrixXd& _extField,
                                      const Eigen::MatrixXd& _source=Eigen::MatrixXd(),
                                      const Eigen::VectorXi& _face=Eigen::VectorXi()){
  
    assert(_extField.cols()==3*N);
    
    extField=_extField;
    source = _source;
    face = _face;
    
   
    if (face.rows()==0){
      intField.resize(extField.rows(),2*N);
      for (int i=0;i<N;i++)
        intField.block(0,2*i,intField.rows(),2)<<(extField.block(0,3*i,extField.rows(),3).array()*mesh->Bx.array()).rowwise().sum(),(extField.block(0,3*i,extField.rows(),3).array()*mesh->By.array()).rowwise().sum();
    } /*else {
      intField.resize(face.rows(),2*N);
      for (int i=0;i<N;i++)
        for (int j=0;j<face.rows();j++)
          intField.block(0,2*i,1,2)<<(extField.block(0,3*i,1,3).array()*mesh->Bx.row(face(j)).array()).sum(),(extField.block(0,3*i,1,3).array()*mesh->By.row(face(j)).array()).sum();*/
    
    //TODO: project extrinsic vector on intrinsic face
    }
  }
  
  void IGL_INLINE set_intrinsic_field(const Eigen::MatrixXd& _intField){
    assert (!(fieldType==POWER_FIELD) || (_intField.cols()==2));
    assert ((_intField.cols()==2*N) || !(fieldType==POLYVECTOR_FIELD || fieldType==RAW_FIELD));
    intField = _intField;
    
    //computing extrinsic field
    extField.conservativeResize(intField.rows(),intField.cols()*3/2);
    for (int i=0;i<intField.rows();i++)
      for (int j=0;j<intField.cols();j+=2)
        extField.block(i,3*j/2,1,3)=mesh->Bx.row(i)*intField(i,j)+mesh->By.row(i)*intField(i,j+1);
  }
  
  void IGL_INLINE set_intrinsic_field(const Eigen::MatrixXcd& _intField){
    assert (!(fieldType==POWER_FIELD) || (_intField.cols()==1));
    assert ((_intField.cols()==N) || !(fieldType==POLYVECTOR_FIELD || fieldType==RAW_FIELD));
    
    intField.resize(_intField.rows(),_intField.cols()*2);
    for (int i=0;i<_intField.cols();i++){
      intField.col(2*i)=_intField.col(i).real();
      intField.col(2*i+1)=_intField.col(i).imag();
    }
    //computing extrinsic field
    extField.conservativeResize(_intField.rows(),_intField.cols()*3);
    for (int i=0;i<_intField.rows();i++)
      for (int j=0;j<_intField.cols();j++)
        extField.block(i,3*j,1,3)=mesh->Bx.row(i)*_intField(i,j).real()+mesh->By.row(i)*_intField(i,j).imag();
  }
  
  Eigen::MatrixXd  virtual IGL_INLINE project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const{
    assert(tangentSpaces.rows()==extDirectionals.rows());
    Eigen::MatrixXd intDirectionals(tangentSpaces.rows(),2);
    
    //std::cout<<"tangentSpaces: "<<tangentSpaces<<std::endl;
    for (int j=0;j<tangentSpaces.rows();j++){
      /*std::cout<<"j: "<<j<<std::endl;
      std::cout<<"(extDirectionals.row(j).array()*mesh->Bx.row(tangentSpaces(j)).array()).sum(),(extDirectionals.row(j).array()*mesh->By.row(tangentSpaces(j)).array()).sum():"<<(extDirectionals.row(j).array()*mesh->Bx.row(tangentSpaces(j)).array()).sum()<<","<<(extDirectionals.row(j).array()*mesh->By.row(tangentSpaces(j)).array()).sum()<<std::endl;*/
      intDirectionals.row(j)<<(extDirectionals.row(j).array()*mesh->Bx.row(tangentSpaces(j)).array()).sum(),(extDirectionals.row(j).array()*mesh->By.row(tangentSpaces(j)).array()).sum();
    }
    return intDirectionals;
  }
    
  void IGL_INLINE set_singularities(const Eigen::VectorXi& _singVertices,
                                    const Eigen::VectorXi& _singIndices){
    
    Eigen::VectorXi vertexIndices=Eigen::VectorXi::Zero(mesh->V.rows());
    int singCounter=_singVertices.size();
    for (int i=0;i<_singVertices.size();i++)
      vertexIndices(_singVertices(i))=_singIndices(i);
    
    //removing boundary indices
    std::vector<std::vector<int> > L;
    igl::boundary_loop(mesh->F, L);
    for (int j=0;j<L.size();j++)
      for (int k=0;k<L[j].size();k++){
        if (vertexIndices(L[j][k])) singCounter--;
        vertexIndices(L[j][k])=0;
      }
    
    singCycles.resize(singCounter);
    singIndices.resize(singCounter);
    singCounter=0;
    for (int i=0;i<mesh->V.rows();i++){
      if (vertexIndices(i)){
        singCycles(singCounter)=i;
        singIndices(singCounter++)=vertexIndices(i);
      }
    }
  }
  

};

}



#endif
