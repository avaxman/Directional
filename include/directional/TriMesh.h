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
        Eigen::MatrixXi EF, FE, EV,TT, EFi;
        Eigen::MatrixXd FEs;
        Eigen::VectorXi innerEdges, boundEdges, vertexValence;  //vertexValence is #(outgoing edges) (if boundary, then #faces+1 = vertexvalence)
        Eigen::VectorXi isBoundaryVertex, isBoundaryEdge;

        //DCEL quantities
        Eigen::VectorXi VH,HV,HE,HF,nextH,prevH,twinH;
        Eigen::MatrixXi EH,FH, VE, VF;

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
            HE.resize(3*F.rows());
            nextH.resize(3*F.rows());
            prevH.resize(3*F.rows());
            HV.resize(3*F.rows());
            VH.resize(V.rows());
            HF.resize(3*F.rows());
            FH.resize(F.rows(), 3);
            twinH = Eigen::VectorXi::Constant(3*F.rows(),-1);
            Eigen::MatrixXi halfedges(3*F.rows(),2);
            for (int i=0;i<F.rows();i++) {
                halfedges.block(3 * i, 0, 3, 2) << F(i, 0), F(i, 1),
                        F(i, 1), F(i, 2),
                        F(i, 2), F(i, 0);
                nextH.segment(3*i,3)<<3*i+1, 3*i+2, 3*i;
                prevH.segment(3*i,3)<<3*i+2, 3*i, 3*i+1;
                HV.segment(3*i,3)<<F.row(i).transpose();
                HF.segment(3*i,3)<<i,i,i;
                FH.row(i)<<3*i, 3*i+1, 3*i+2;
            }

            //finding twins
            typedef std::pair<std::pair<int, int>, int> pairPlusOne;
            std::set<pairPlusOne, ComparePairs> edgeSet;
            std::vector<int> EHList;
            for (int i=0;i<halfedges.rows();i++){
                VH(HV(i))=i;
                std::pair<int,int> oppEdge(halfedges(i,1), halfedges(i,0));
                pairPlusOne oppEdgePlus(oppEdge, -1);
                std::set<pairPlusOne>::iterator si = edgeSet.find(oppEdgePlus);
                if (si == edgeSet.end()) {
                    edgeSet.insert(pairPlusOne(std::pair<int, int>(halfedges(i, 0), halfedges(i, 1)), i));
                    EHList.push_back(i);
                } else {  //found matching twin
                    twinH[si->second] = i;
                    twinH[i] = si->second;
                }
            }

            //std::cout<<"twinH: "<<twinH<<std::endl;
            //creating the edge quantities from the halfedge quantities
            EH = Eigen::Map<Eigen::VectorXi>(EHList.data(), EHList.size());
            HE.resize(halfedges.rows());
            EV.resize(EHList.size(),2);
            EF = Eigen::MatrixXi::Constant(EHList.size(),2,-1);
            EFi.resize(EHList.size(),2);
            //VE.resize(V.size());
            FE.resize(F.size(),3);
            FEs.resize(F.size(),3);
            TT.resize(F.size(),3);
            Eigen::VectorXi HEs(halfedges.rows());
            for (int i=0;i<EH.size();i++){
                EV.row(i)<<halfedges(EH(i)),halfedges(nextH(EH(i)));
                EF(i,0) = HF(EH(i));
                if (twinH(EH(i))!=-1)
                    EF(i,1) = HF(twinH(EH(i)));

                HE(EH(i))=i;
                HEs(EH(i))=1;
                if (twinH((EH(i)))!=-1) {
                    HE(twinH(EH(i))) = i;
                    HEs(twinH(EH(i))) = -1;
                }

                EFi(i,0) = (EH(i) + 0) % 3;
                if (twinH((EH(i)))!=-1)
                    EFi(i,1) = (twinH(EH(i)) + 0) % 3;
            }
            for (int i=0;i<FH.rows();i++){
                for (int j=0;j<3;j++){
                    FE(i,j) = HE(FH(i,j));
                    //checking the orientation of the edge vs. the halfedge
                    FEs(i,j) = HEs(FH(i,j));

                    if (twinH(FH(i,j)) == -1)
                        TT(i, j) = -1;
                    else
                        TT(i, j) = HF(twinH(FH(i,j)));
                }
            }

            //gathering all boundary halfedges
            std::vector<int> innerEdgesList, boundEdgesList, boundHalfedgesList;
            isBoundaryVertex=Eigen::VectorXi::Zero(V.size());
            isBoundaryEdge=Eigen::VectorXi::Zero(EV.size());
            for (int i = 0; i < halfedges.rows(); i++) {
                if (twinH(i) == -1){
                    boundEdgesList.push_back(HE(i));
                    boundHalfedgesList.push_back(i);
                    isBoundaryEdge(HE(i))=1;
                    isBoundaryVertex(EV(HE(i),0))=1;
                    isBoundaryVertex(EV(HE(i),1))=1;
                }else
                    innerEdgesList.push_back(i);
            }


            innerEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(innerEdgesList.data(), innerEdgesList.size());
            boundEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(boundEdgesList.data(), boundEdgesList.size());

            //creating boundary loops
            Eigen::VectorXi isVisited = Eigen::VectorXi::Constant(halfedges.size(),0);
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
                    currBoundaryLoop.push_back(HV(currHE));
                    isVisited(currHE) = 1;
                    //finding the next boundary halfedge
                    currHE =nextH(currHE);
                    while (twinH(currHE)!=-1)
                        currHE = nextH(twinH(currHE));
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


            //computing extra combinatorial information
            //Relative location of edges within faces
            /*EFi = Eigen::MatrixXi::Constant(EF.rows(), 2, -1); // number of an edge inside the face
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
            }*/


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
