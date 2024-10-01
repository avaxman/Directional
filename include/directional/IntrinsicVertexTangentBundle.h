// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_INTRINSIC_VERTEX_TANGENT_BUNDLE_H
#define DIRECTIONAL_INTRINSIC_VERTEX_TANGENT_BUNDLE_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/sparse_diagonal.h>
#include <directional/dual_cycles.h>
#include <directional/TriMesh.h>
#include <directional/raw_to_polyvector.h>
#include <directional/polyvector_to_raw.h>

/***
 This class implements intrinsic vertex-based fields (defined intrinsically on the curved corners of a 1-ring), where the connection is on edges between neighboring vertices,
 the local cycles are triangles, and the curvature identifies as the holonomy of the parallel transport around this triangle.
 ***/


namespace directional{

//This is the interface class for any directional fields represented in cartesian coordinates, of any order N.
    class IntrinsicVertexTangentBundle : public TangentBundle{
    public:

        const TriMesh* mesh;
        Eigen::MatrixXd tangentStartAngles;  //where each edge begins on the intrinsic space

        virtual discTangTypeEnum discTangType() const {return discTangTypeEnum::VERTEX_SPACES;}
        virtual bool hasCochainSequence() const { return false; }
        virtual bool hasEmbedding() const { return true; }

        IntrinsicVertexTangentBundle(){}
        ~IntrinsicVertexTangentBundle(){}

        void inline init(const TriMesh& _mesh){

            intDimension = 2;
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
                double totalTangentSum = (mesh->isBoundaryVertex(i) ? directional::PI : 2.0*directional::PI);
                double angleSum =  totalTangentSum - mesh->GaussianCurvature(i);
                tangentStartAngles.col(0).setZero();  //the first angle
                Eigen::RowVector3d prevEdgeVector = mesh->V.row(mesh->HV(mesh->nextH(mesh->VH(i))))-mesh->V.row(i);
                int hebegin = mesh->VH(i);  //this should be the first boundary edge in case of boundary
                int heiterate = mesh->twinH(mesh->prevH(hebegin));
                int j=1;
                do{
                    Eigen::RowVector3d currEdgeVector = mesh->V.row(mesh->HV(mesh->nextH(heiterate)))-mesh->V.row(i);
                    double angleDiff = std::acos(currEdgeVector.dot(prevEdgeVector)/(prevEdgeVector.norm()*currEdgeVector.norm()));
                    tangentStartAngles(i,j)=tangentStartAngles(i,j-1)+totalTangentSum*angleDiff/angleSum;
                    heiterate = mesh->twinH(mesh->prevH(heiterate));
                    j++;
                    prevEdgeVector=currEdgeVector;
                }while ((heiterate!=hebegin)&&(heiterate!=-1));
            }

            //connection is the ratio of the complex representation of mutual edges
            connection.resize(mesh->EV.rows(),1);  //the difference in the angle representation of edge i from EV(i,0) to EV(i,1)
            Eigen::MatrixXd edgeVectors(mesh->EV.rows(), 3);
            for (int i = 0; i < mesh->EV.rows(); i++) {
                //edgeVectors.row(i) = (mesh->V.row(mesh->EV(i, 1)) - mesh->V.row(mesh->EV(i, 0))).normalized();
                //looking up edge in each tangent space
                Complex ef,eg;
                for (int j=0;j<mesh->vertexValence(mesh->EV(i,0));j++){
                    if (mesh->VE(mesh->EV(i,0),j)==i)
                        ef = exp(Complex(0,tangentStartAngles(mesh->EV(i, 0),j)));
                }

                for (int j=0;j<mesh->vertexValence(mesh->EV(i,1));j++) {
                  if (mesh->VE(mesh->EV(i, 1), j) == i)
                      eg = exp(Complex(0, tangentStartAngles(mesh->EV(i, 1), j)));
              }

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
                    cyclesTriplets.push_back(Eigen::Triplet<double>(i,mesh->FE(i,j),mesh->FEs(i,j)));
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

            /************masses****************/
            connectionMass.resize(mesh->EV.rows(), mesh->EV.rows());
            tangentSpaceMass.resize(mesh->V.rows(), mesh->V.rows());

            std::vector<Eigen::Triplet<double>> connectionMassTris, tsMassTris;
            Eigen::VectorXd cotWeights = Eigen::VectorXd::Zero(mesh->EV.rows());

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

                    cotWeights(mesh->FE(i,j))+=0.5*faceCotWeights(i,j);
                }
            }

            connectionMass = directional::sparse_diagonal(cotWeights);


            //masses are vertex voronoi areas
            Eigen::VectorXd voronoiMass = Eigen::VectorXd::Zero(mesh->V.rows());
            for (int i=0;i<mesh->F.rows();i++)
                for (int j=0;j<3;j++)
                    voronoiMass(mesh->F(i,j)) += mesh->faceAreas(i)/3.0;

            tangentSpaceMass = directional::sparse_diagonal(voronoiMass);

        }


        //projecting an arbitrary set of extrinsic vectors (e.g. coming from user-prescribed constraints) into intrinsic vectors.
        Eigen::MatrixXd  virtual inline project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const{
            assert(tangentSpaces.rows()==extDirectionals.rows());// || tangentSpaces.rows()==0);

           /* Eigen::VectorXi actualTangentSpaces;
            if (tangentSpaces.rows()==0)
                actualTangentSpaces = Eigen::VectorXi::LinSpaced(sources.rows(), 0, sources.rows()-1);
            else
                actualTangentSpaces = tangentSpaces;*/

            int N = extDirectionals.cols()/3;
            Eigen::MatrixXd intDirectionals(tangentSpaces.rows(),2*N);

            for (int i=0;i<tangentSpaces.rows();i++){
                for (int j=0;j<N;j++)
                    intDirectionals.block(i,2*j,1,2)<<(extDirectionals.block(i,3*j,1,3).array()*mesh->VBx.row(tangentSpaces(i)).array()).sum(),(extDirectionals.block(i,3*j,1,3).array()*mesh->VBy.row(tangentSpaces(i)).array()).sum();
            }
            return intDirectionals;
        }


        //projecting intrinsic to extrinsicc
        Eigen::MatrixXd virtual inline project_to_extrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& intDirectionals) const {

            assert(tangentSpaces.rows()==intDirectionals.rows() || tangentSpaces.rows()==0);
            Eigen::VectorXi actualTangentSpaces;
            if (tangentSpaces.rows()==0)
                actualTangentSpaces = Eigen::VectorXi::LinSpaced(sources.rows(), 0, sources.rows()-1);
            else
                actualTangentSpaces = tangentSpaces;

            int N = intDirectionals.cols()/2;
            Eigen::MatrixXd extDirectionals(actualTangentSpaces.rows(),3);

            extDirectionals.conservativeResize(intDirectionals.rows(),intDirectionals.cols()*3/2);
            for (int i=0;i<intDirectionals.rows();i++)
                for (int j=0;j<intDirectionals.cols();j+=2)
                    extDirectionals.block(i,3*j/2,1,3)=mesh->VBx.row(actualTangentSpaces(i))*intDirectionals(i,j)+mesh->VBy.row(actualTangentSpaces(i))*intDirectionals(i,j+1);

            return extDirectionals;
        }


        //Primitive version that just interpolates inside the triangle - should be changed to a full linear version!
        void inline interpolate(const Eigen::MatrixXi &faceIndices,
                                    const Eigen::MatrixXd &baryCoords,
                                    const Eigen::MatrixXd &intDirectionals,
                                    Eigen::MatrixXd& interpSources,
                                    Eigen::MatrixXd& interpNormals,
                                    Eigen::MatrixXd& interpField) const {

            using namespace Eigen;
            assert(faceIndices.rows()==baryCoords.rows());
            //assert(baryCoords.rows()==intDirectionals.rows());

            int N = intDirectionals.cols()/2;
            interpSources=Eigen::MatrixXd::Zero(faceIndices.rows(),3);
            interpNormals=Eigen::MatrixXd::Zero(faceIndices.rows(),3);
            interpField=Eigen::MatrixXd::Zero(faceIndices.rows(),3*N);

            //interpolating through the polyvector representation to avoid singular issues

            //projecting extrinsic vectors unto corners (again a bit primitive)
            //idea: find rotation matrix vertex->face and then just use intrinsic vector projection without going through extrinsic
            //then interpolate via PV in the center

            MatrixXd cornerField(3*mesh->F.rows(),2*N);

            for (int i=0;i<mesh->F.rows();i++)
                for (int j=0;j<3;j++){
                    MatrixXd VBasis(2,3), FBasis(2,3);
                    VBasis.row(0)=mesh->VBx.row(mesh->F(i,j));
                    VBasis.row(1)=mesh->VBy.row(mesh->F(i,j));
                    FBasis.row(0) =mesh->FBx.row(i);
                    FBasis.row(1)=mesh->FBy.row(i);
                    Matrix2d basisChange = VBasis*FBasis.transpose();
                    for (int k=0;k<N;k++)
                        cornerField.block(3*i+j,2*k,1,2)=intDirectionals.block(mesh->F(i,j),2*k,1,2)*basisChange;
                }

            //std::cout<<"cornerField: "<<cornerField<<std::endl;
            Eigen::MatrixXcd pvDirectionals;
            raw_to_polyvector(cornerField, N, pvDirectionals);
            //std::cout<<"pvDirectionals: "<<pvDirectionals<<std::endl;
            Eigen::MatrixXcd interpPV=Eigen::MatrixXd::Zero(faceIndices.rows(), pvDirectionals.cols());

            //converting to polyvector to blend
            for (int i=0;i<faceIndices.rows();i++) {

                for (int j = 0; j < 3; j++) {
                    interpSources.row(i).array() += mesh->V.row(mesh->F(faceIndices(i), j)).array() * baryCoords(i, j);
                    interpPV.row(i).array() += pvDirectionals.row(mesh->F(faceIndices(i), j)).array() * baryCoords(i, j);
                }

                interpNormals.row(i) = mesh->faceNormals.row(faceIndices(i));
            }

            //std::cout<<"interpPV: "<<interpPV<<std::endl;

            Eigen::MatrixXcd interpFieldInt;
            Eigen::MatrixXd interpPVReal(interpPV.rows(), 2*interpPV.cols());
            for (int j=0;j<interpPV.cols();j++) {
                interpPVReal.col(2 * j) = interpPV.col(j).real();
                interpPVReal.col(2 * j+1) = interpPV.col(j).imag();
            }

            polyvector_to_raw(interpPVReal, N, interpFieldInt);

            for (int i=0;i<faceIndices.rows();i++)
                //This is on faces, so using face-based base
                for (int j=0;j<N;j++)
                    interpField.block(i,j*3,1,3)=interpFieldInt(faceIndices(i),j).real()*mesh->FBx.row(faceIndices(i))+interpFieldInt(faceIndices(i),j).imag()*mesh->FBy.row(faceIndices(i));

        }
    };
}



#endif
  

