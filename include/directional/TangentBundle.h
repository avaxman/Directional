// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_TANGENT_BUNDLE_H
#define DIRECTIONAL_TANGENT_BUNDLE_H

#include <cassert>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>

/***
 This is the interface class for discrete tangent bundles, represented as a graph G_TB of nodes (tangent spaces) and edges (adjacency relations). This represents a d-dimensional tangent bundle.
 Further structure is cycles around which holonomy is measured, and curvature is defined, and consequently one-rings of all edges around a single node.
 The Tangent Bundle is "intrinsic" in the sense that the information required for designing fields does not require any embedding. The class includes variables that include embedding information like sources and normals (assuming the bundle is at most 1-codimensional).
 They are however only used for input/output to the intrinsic variables.
 ***/

namespace directional{

enum class discTangTypeEnum {BASE_CLASS, FACE_SPACES, VERTEX_SPACES};
enum class boundCondTypeEnum {DIRICHLET, NEUMANN};

class TangentBundle {
public:
    
    //In case some methods are only defined for specific tangent bundles
    virtual discTangTypeEnum discTangType() const { return discTangTypeEnum::BASE_CLASS; }
    
    //indicates if the tangent bundle has a meaningful embedding
    virtual bool hasEmbedding() const { return false; }
    
    int intDimension;                                   //Intrinsic dimension of the tangent spaces
    int numSpaces;                                      //The number of tangent spaces
    
    //combinatorics and topology
    Eigen::Matrix<int, Eigen::Dynamic, 2> adjSpaces;    //Adjacent tangent spaces (edges of the tangent bundle graph). This includes boundaries! one side will be "-1" when this is the case. Also see "innerAdjacencies"
    Eigen::MatrixXi oneRing;                            //(Ordered) tangent spaces adjacent around each tangent space
    Eigen::VectorXi innerAdjacencies;                   //Indices into adjSpaces that are not boundary
    Eigen::SparseMatrix<double> cycles;                 //Adjaceny matrix of cycles
    Eigen::VectorXd cycleCurvatures;                    //Curvature of cycles.
    Eigen::VectorXi local2Cycle;                        //Map between local cycles and general cycles
    
    //Geometry
    //the connection between adjacent tangent space. That is, a field is parallel between adjaspaces(i,0) and adjSpaces(i,1) when complex(intField.row(adjSpaceS(i,0))*connection(i))=complex(intField.row(adjSpaceS(i,1))
    Eigen::VectorXcd connection;                                    //#V, Metric connection between adjacent spaces
    Eigen::SparseMatrix<double> connectionMass;                     //The mass matrix of connections, of size #adjSpaces x #adjSpaces
    Eigen::SparseMatrix<double> tangentSpaceMass;                   //The inner-product mass for vectors in tangent spaces, of size #V (self masses) + #E (adjSpaces masses;  optional, usually for high-order fields)
    Eigen::SparseMatrix<double> invTangentSpaceMass;
    double avgAdjLength;                                            //The average distance between adjacent spaces
    
    //Extrinsic components
    Eigen::MatrixXd sources;  //the source point of the extrinsic vectors
    Eigen::MatrixXd normals;  //the normals to the tangent spaces (assuming 1-codimension)
    Eigen::MatrixXd cycleSources;  //source point of cycles
    Eigen::MatrixXd cycleNormals;  //normals to cycles
    
    TangentBundle() {}
    ~TangentBundle() {}
    
    //projecting an arbitrary set of extrinsic vectors (e.g. coming from user-prescribed constraints) into intrinsic vectors.
    Eigen::MatrixXd virtual inline project_to_intrinsic(const Eigen::VectorXi &tangentSpaces, const Eigen::MatrixXd &extDirectionals) const {
        assert(false && "The base class does not have an embedding");
        return Eigen::MatrixXd();
    }
    
    //projecting extrinsic to intrinsic
    Eigen::MatrixXd virtual inline project_to_extrinsic(const Eigen::VectorXi &tangentSpaces, const Eigen::MatrixXd &extDirectionals) const {
        assert(false && "The base class does not have an embedding");
        return Eigen::MatrixXd();
    }
    
    //interpolating field from nodes in palces specified by barycentric coordinates. The interpolator needs to interpret them.
    void virtual inline interpolate(const Eigen::MatrixXi &elemIndices,
                                    const Eigen::MatrixXd &baryCoords,
                                    const Eigen::MatrixXd &intDirectionals,
                                    Eigen::MatrixXd& interpSources,
                                    Eigen::MatrixXd& interpNormals,
                                    Eigen::MatrixXd& interpField)  const{
        interpSources=Eigen::MatrixXd();
        interpNormals=Eigen::MatrixXd();
        interpField=Eigen::MatrixXd();
    }
    
    Eigen::SparseMatrix<double> virtual inline curl_matrix(const boundCondTypeEnum boundCondType,
                                                           const Eigen::VectorXi& matching,
                                                           const bool intrinsic = false) const{
        return Eigen::SparseMatrix<double>();
    }
    
};

}


#endif
