// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_TRIMESH_H
#define DIRECTIONAL_TRIMESH_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/barycenter.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/local_basis.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include <igl/boundary_loop.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/avg_edge_length.h>
#include <directional/gaussian_curvature.h>
#include <igl/doublearea.h>
#include <directional/dcel.h>

/***
 This class stores a general-purpose triangle mesh. The triangle mesh can be used to implement a tangent bundle in several ways
 (for instance, face- or vertex-based bundles).
***/

namespace directional{

    class TriMesh{
    public:

        //Basic quantities
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;

        //combinatorial quantities
        Eigen::MatrixXi EF, FE, EV,TT, EFi, VE, VF;
        Eigen::MatrixXd FEs;
        Eigen::VectorXi innerEdges, boundEdges, vertexValence;  //vertexValence is #(outgoing edges) (if boundary, then #faces+1 = vertexvalence)
        Eigen::VectorXi isBoundaryVertex, isBoundaryEdge;

        //DCEL quantities
        Eigen::VectorXi VH,HV,HE,HF,nextH,prevH,twinH;
        Eigen::MatrixXi EH,FH;

        //Geometric quantities
        Eigen::MatrixXd faceNormals;
        Eigen::MatrixXd faceAreas;
        Eigen::MatrixXd vertexNormals;
        Eigen::MatrixXd FBx,FBy;  //local basis vectors per face
        Eigen::MatrixXd VBx,VBy;  //local basis vectors per vertex
        Eigen::MatrixXd barycenters;
        Eigen::VectorXd GaussianCurvature;

        //Measures of the scale of a mesh
        double avgEdgeLength;
        Eigen::RowVector3d minBox, maxBox;   //bounding box

        int eulerChar;
        int numGenerators;

        std::vector<std::vector<int>> boundaryLoops;

        TriMesh(){}
        ~TriMesh(){}

        void IGL_INLINE set_mesh(const Eigen::MatrixXd& _V,
                                 const Eigen::MatrixXi& _F,
                                 const Eigen::MatrixXi& _EV=Eigen::MatrixXi(),
                                 const Eigen::MatrixXi& _FE=Eigen::MatrixXi(),
                                 const Eigen::MatrixXi& _EF=Eigen::MatrixXi()) {

            V = _V;
            F = _F;
            if (_EV.rows() == 0) {
                igl::edge_topology(V, F, EV, FE, EF);
            } else {
                EV = _EV;
                FE = _FE;
                EF = _EF;
            }
            std::vector<int> innerEdgesList, boundEdgesList;
            isBoundaryVertex=Eigen::VectorXi::Zero(V.size());
            isBoundaryEdge=Eigen::VectorXi::Zero(EV.size());
            for (int i = 0; i < EF.rows(); i++) {
                if ((EF(i, 1) == -1) || (EF(i, 0) == -1)) {
                    boundEdgesList.push_back(i);
                    isBoundaryEdge(i)=1;
                    isBoundaryVertex(EV(i,0))=1;
                    isBoundaryVertex(EV(i,1))=1;
                }else
                    innerEdgesList.push_back(i);

            }

            //computing extra combinatorial information
            //Relative location of edges within faces
            EFi = Eigen::MatrixXi::Constant(EF.rows(), 2, -1); // number of an edge inside the face
            for(int i = 0; i < EF.rows(); i++)
            {
                for (int k = 0; k < 2; k++)
                {
                    if (EF(i, k) == -1)
                        continue;
                    for (int j = 0; j < 3; j++)
                        if (FE(EF(i, k), j) == i)
                            EFi(i, k) = j;
                }
            }

            //sign of edge within face
            FEs = Eigen::MatrixXd::Zero(FE.rows(), FE.cols());

            for(int i = 0; i < EF.rows(); i++)
            {
                if(EFi(i, 0) != -1)
                    FEs(EF(i, 0), EFi(i, 0)) = 1.0;
                if(EFi(i,1) != -1)
                    FEs(EF(i, 1), EFi(i, 1)) = -1.0;
            }

            innerEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(innerEdgesList.data(), innerEdgesList.size());
            boundEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(boundEdgesList.data(), boundEdgesList.size());
            igl::barycenter(V, F, barycenters);
            igl::local_basis(V, F, FBx, FBy, faceNormals);
            igl::triangle_triangle_adjacency(F, TT);
            igl::boundary_loop(F, boundaryLoops);
            directional::gaussian_curvature(V,F,isBoundaryVertex, GaussianCurvature);
            igl::doublearea(V,F,faceAreas);
            faceAreas.array()/=2.0;
            eulerChar = V.rows() - EV.rows() + F.rows();
            numGenerators = (2 - eulerChar)/2 - boundaryLoops.size();
            avgEdgeLength=igl::avg_edge_length(V,F);
            minBox = V.colwise().minCoeff();
            maxBox = V.colwise().maxCoeff();

            hedra::dcel(Eigen::VectorXi::Constant(F.rows(),3),F,EV,EF,EFi, innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH);
            vertexValence=Eigen::VectorXi::Zero(V.rows());
            for (int i=0;i<EV.rows();i++){
                vertexValence(EV(i,0))++;
                vertexValence(EV(i,1))++;
            }

            //TODO: adapt to boundaries
            VE.resize(V.rows(),vertexValence.maxCoeff());
            VF.resize(V.rows(),vertexValence.maxCoeff());
            for (int i=0;i<V.rows();i++){
                int counter=0;
                int hebegin = VH(i);
                if (isBoundaryVertex(i)) //winding up hebegin to the first boundary edge
                    while (twinH(hebegin)!=-1)
                        hebegin = nextH(twinH(hebegin));

                //resetting VH for future reference
                VH(i)=hebegin;
                int heiterate = hebegin;
                do {
                    VE(i, counter) = HE(heiterate);
                    VF(i, counter++) = HF(heiterate);
                    if (twinH(prevH(heiterate))==-1) { //last edge before end, adding the next edge
                        VE(i, counter) = HE(prevH(heiterate));  //note counter is already ahead
                        break;
                    }
                    heiterate = twinH(prevH(heiterate));
                }while(hebegin!=heiterate);
            }

            //computing vertex normals by area-weighted aveage of face normals
            vertexNormals=Eigen::MatrixXd::Zero(V.rows(),3);
            for (int i=0;i<F.rows();i++)
                for (int j=0;j<3;j++)
                    vertexNormals.row(F(i,j)).array()+=faceNormals.row(i).array()*faceAreas(i);

            vertexNormals.rowwise().normalize();

            //computing local basis that aligns with the first projected edge of each triangle
            VBx.resize(V.rows(),3);
            VBy.resize(V.rows(),3);
            for (int i=0;i<V.rows();i++){
                Eigen::RowVector3d firstEdge = V.row(HV(nextH(VH(i))))-V.row(i);
                VBx.row(i)=firstEdge-(firstEdge.dot(vertexNormals.row(i)))*vertexNormals.row(i);
                VBx.row(i).normalize();
                Eigen::RowVector3d currx=VBx.row(i);
                Eigen::RowVector3d currn=vertexNormals.row(i);
                VBy.row(i)=currn.cross(currx);
                VBy.row(i).normalize();
            }

        }

    };

}



#endif /* DIRECTIONAL_TRIMESH_H */
