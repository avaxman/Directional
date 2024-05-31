// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_TRIMESH_H
#define DIRECTIONAL_TRIMESH_H

#include <iostream>
#include <vector>
#include <set>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/gaussian_curvature.h>
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
        /*Eigen::VectorXi VH,HV,HE,HF,nextH,prevH,twinH;
        Eigen::MatrixXi EH,FH, VE, VF;*/
        typedef DCEL<int,int,int,int> TriMeshDCEL;
        TriMeshDCEL dcel;

        //writing functions that shorten code
        inline int VH(const int index) const {return dcel.vertices[index].halfedge;}
        inline int twinH(const int index) const {return dcel.halfedges[index].twin;}
        inline int nextH(const int index) const {return dcel.halfedges[index].next;}
        inline int prevH(const int index) const {return dcel.halfedges[index].prev;}
        inline int HV(const int index) const {return dcel.halfedges[index].vertex;}
        inline int HF(const int index) const {return dcel.halfedges[index].face;}
        inline int HE(const int index) const {return dcel.halfedges[index].edge;}
        inline int EH(const int index, const int side) const {return (side==0 ? dcel.edges[index].halfedge : dcel.halfedges[dcel.edges[index].halfedge].twin);}
        inline int FH(const int index) const {return dcel.faces[index].halfedge;}

        //Geometric quantities
        Eigen::MatrixXd faceNormals;
        Eigen::VectorXd faceAreas;
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

        //computing a full HE structure
        void inline compute_edge_quantities(){

            struct ComparePairs {
                bool operator()(const std::pair<std::pair<int, int>, int>& a, const std::pair<std::pair<int, int>, int>& b) const {
                    if (a.first.first == b.first.first) {
                        return a.first.second < b.first.second;
                    } else {
                        return a.first.first < b.first.first;
                    }
                }
            };

            //This is done in the polyscope compatible fashion
            dcel.init(V, F);

            EV.resize(dcel.edges.size(),2);
            EF = Eigen::MatrixXi::Constant(dcel.edges.size(),2,-1);
            EFi.resize(dcel.edges.size(),2);
            FE.resize(F.rows(),3);
            FEs.resize(F.rows(),3);
            TT.resize(F.rows(),3);
            for (int i=0;i<dcel.edges.size();i++){
                EV.row(i)<<dcel.halfedges[dcel.edges[i].halfedge].vertex,
                        dcel.halfedges[dcel.halfedges[dcel.edges[i].halfedge].next].vertex;
                EF(i,0) = dcel.halfedges[dcel.edges[i].halfedge].face;
                if (dcel.halfedges[dcel.edges[i].halfedge].twin!=-1)
                    EF(i,1) = dcel.halfedges[dcel.halfedges[dcel.edges[i].halfedge].twin].vertex;

                EFi(i,0) = (dcel.edges[i].halfedge + 0) % 3;
                if (dcel.halfedges[dcel.edges[i].halfedge].twin !=-1)
                    EFi(i,1) = (dcel.halfedges[dcel.edges[i].halfedge].twin + 0) % 3;
            }
            for (int i=0;i<dcel.faces.size();i++){
                int hebegin = dcel.faces[i].halfedge;
                int heiterate = hebegin;
                int inFaceCounter=0;
                do{
                    FE(i,inFaceCounter) = dcel.halfedges[heiterate].edge;
                    FEs(i,inFaceCounter) = (dcel.edges[dcel.halfedges[heiterate].edge].halfedge==heiterate ? 1 : -1);
                    if (dcel.halfedges[heiterate].twin!=-1)
                        TT(i, inFaceCounter) = dcel.halfedges[dcel.halfedges[heiterate].twin].face;
                    inFaceCounter++;
                    heiterate = dcel.halfedges[heiterate].next;
                }while (heiterate!=hebegin);
            }

            //gathering all boundary halfedges
            std::vector<int> innerEdgesList, boundEdgesList, boundHalfedgesList;
            isBoundaryVertex=Eigen::VectorXi::Zero(V.size());
            isBoundaryEdge=Eigen::VectorXi::Zero(EV.size());
            for (int i = 0; i < dcel.halfedges.size(); i++) {
                if (dcel.halfedges[i].twin == -1){
                    boundEdgesList.push_back(dcel.halfedges[i].edge);
                    boundHalfedgesList.push_back(i);
                    isBoundaryEdge(dcel.halfedges[i].edge)=1;
                    isBoundaryVertex(EV(dcel.halfedges[i].edge,0))=1;
                    isBoundaryVertex(EV(dcel.halfedges[i].edge,1))=1;
                }else
                    innerEdgesList.push_back(i);
            }

            innerEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(innerEdgesList.data(), innerEdgesList.size());
            boundEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(boundEdgesList.data(), boundEdgesList.size());

            //creating boundary loops
            Eigen::VectorXi isVisited = Eigen::VectorXi::Constant(dcel.halfedges.size(),0);
            while (isVisited.sum()!=boundHalfedgesList.size()){
                //choose a first one
                int beginHE;
                for (beginHE=0;beginHE<boundHalfedgesList.size();beginHE++)
                    if (isVisited(boundHalfedgesList[beginHE])==0)
                        break;

                beginHE = boundHalfedgesList[beginHE];  //in the global indexing
                int currHE = beginHE;
                std::vector<int> currBoundaryLoop;
                do{
                    currBoundaryLoop.push_back(dcel.halfedges[currHE].vertex);
                    isVisited(currHE) = 1;
                    //finding the next boundary halfedge
                    currHE = dcel.halfedges[currHE].next;
                    while (dcel.halfedges[currHE].twin!=-1)
                        currHE = dcel.halfedges[dcel.halfedges[currHE].twin].next;
                }while(currHE!=beginHE);
                boundaryLoops.push_back(currBoundaryLoop);
            }

        }

        void inline compute_geometric_quantities(){

            //barycenters and local bases
            barycenters.resize(F.rows(),3);
            FBx.resize(F.rows(),3);
            FBy.resize(F.rows(),3);
            faceNormals.resize(F.rows(),3);
            faceAreas.resize(F.rows());
            for (int i=0;i<F.rows();i++) {
                barycenters.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3;
                Eigen::RowVector3d localx = V.row(F(i,1))-V.row(F(i,0));
                Eigen::RowVector3d localy = V.row(F(i,2))-V.row(F(i,0));
                Eigen::RowVector3d localz = localx.cross(localy);
                faceAreas(i) = localz.norm()/2.0;
                localy = localz.cross(localx);
                FBx.row(i) = localx.normalized();
                FBy.row(i) = localy.normalized();
                faceNormals.row(i) = localz.normalized();
            }

            //igl::triangle_triangle_adjacency(F, TT);
            gaussian_curvature(V,F,isBoundaryVertex, GaussianCurvature);

            //Average edge length
            double sumEdgeLength=0.0;
            for (int i=0;i<EV.rows();i++)
                sumEdgeLength += (V.row(EV(i,0))-V.row(EV(i,1))).norm();

            avgEdgeLength =sumEdgeLength/(double)EV.rows();
        }

        void inline set_mesh(const Eigen::MatrixXd& _V,
                             const Eigen::MatrixXi& _F,
                             const Eigen::MatrixXi& _EV=Eigen::MatrixXi(),
                             const Eigen::MatrixXi& _FE=Eigen::MatrixXi(),
                             const Eigen::MatrixXi& _EF=Eigen::MatrixXi()) {

            V = _V;
            F = _F;
            if (_EV.rows() == 0) {
                compute_edge_quantities();
            } else {
                EV = _EV;
                FE = _FE;
                EF = _EF;
            }

            compute_geometric_quantities();

            eulerChar = V.rows() - EV.rows() + F.rows();
            numGenerators = (2 - eulerChar)/2 - boundaryLoops.size();
            minBox = V.colwise().minCoeff();
            maxBox = V.colwise().maxCoeff();

            //hedra::dcel(Eigen::VectorXi::Constant(F.rows(),3),F,EV,EF,EFi, innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH,twinH);
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
                int hebegin = dcel.vertices[i].halfedge;
                if (isBoundaryVertex(i)) //winding up hebegin to the first boundary edge
                    while (dcel.halfedges[hebegin].twin!=-1)
                        hebegin = dcel.halfedges[dcel.halfedges[hebegin].twin].next;

                //TODO: this rest inside the init function
                //resetting dcel pointer for future reference
                dcel.vertices[i].halfedge=hebegin;
                int heiterate = hebegin;
                do {
                    VE(i, counter) = dcel.halfedges[heiterate].edge;
                    VF(i, counter++) = dcel.halfedges[heiterate].face;
                    if (dcel.halfedges[dcel.halfedges[heiterate].prev].twin==-1){
                        VE(i, counter) =  dcel.halfedges[dcel.halfedges[heiterate].prev].edge;  //note counter is already ahead
                        break;
                    }
                    heiterate = dcel.halfedges[dcel.halfedges[heiterate].prev].twin ;
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
                Eigen::RowVector3d firstEdge = V.row(dcel.halfedges[dcel.halfedges[dcel.vertices[i].halfedge].next].vertex)-V.row(i);
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
