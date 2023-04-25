// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>, 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <iostream>
#include <iomanip>
#include <map>
#include <random>
#include <Eigen/Geometry>
#include <igl/edge_topology.h>
#include <igl/sort_vectors_ccw.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/barycenter.h>
#include <igl/slice.h>
#include <igl/speye.h>
#include <igl/avg_edge_length.h>
#include <directional/TriMesh.h>
#include <directional/principal_matching.h>
#include <directional/streamlines.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/IntrinsicVertexTangentBundle.h>

namespace Directional {
    IGL_INLINE void generate_sample_locations(const directional::TriMesh& mesh,
                                              const double distRatio,
                                              Eigen::VectorXi& sampleTris,
                                              Eigen::MatrixXd& samplePoints)
    {

        double minDist = distRatio*igl::avg_edge_length(mesh.V,mesh.F);
        double minDist2=minDist*minDist;


        //first creating a pool of samples to reject from
        //Eigen::VectorXd sumAreas;
        //igl::cumsum(mesh.faceAreas, 0, sumAreas);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::vector<double> faceAreasVec(mesh.faceAreas.data(), mesh.faceAreas.data() + mesh.faceAreas.size());
        std::discrete_distribution<int> distTriangles(faceAreasVec.begin(), faceAreasVec.end());
        std::uniform_real_distribution<double> distBarycentrics(0.0,1.0);
        std::map<int, double> map;
        int nsamples = (int)(10/distRatio)*mesh.F.rows();
        int indirection = 2*distRatio;  //in how far away from the one ring to look

        std::vector<std::vector<Eigen::Vector3d> > samplePool(mesh.F.size());
        std::vector<std::vector<bool>> samplePoolAlive(mesh.F.size());
        //std::cout<<"SamplePool: "<<std::endl;
        for (int i=0;i<nsamples;i++) {
            //random triangle according to area weighting
            int faceIndex = distTriangles(gen);
            //random barycentric coordinate uniformly
            double B1 = distBarycentrics(gen);
            double B2 = distBarycentrics(gen);
            double B3 = distBarycentrics(gen);
            double sum = B1+B2+B3;
            //assuming the unlikely case where they are all zero
            //TODO: convert to actual locations
            Eigen::Vector3d sampleLocation = mesh.V.row(mesh.F(faceIndex,0))*B1/sum+
                                             mesh.V.row(mesh.F(faceIndex,1))*B2/sum+
                                             mesh.V.row(mesh.F(faceIndex,2))*B3/sum;
            samplePool[faceIndex].push_back(sampleLocation);
            samplePoolAlive[faceIndex].push_back(true);
            //std::cout<<faceIndex<<":"<<sampleLocation.transpose()<<std::endl;
        }

        //Choosing samples and deleting everything too close

        //computing "indirection" level adjacencies
        std::vector<Eigen::Triplet<int>> adjTris;
        for (int i=0;i<mesh.EF.rows();i++)
            if ((mesh.EF(i,0)!=-1)&&(mesh.EF(i,1)!=-1)){
                adjTris.push_back(Eigen::Triplet<int>(mesh.EF(i,0), mesh.EF(i,1),1));
                adjTris.push_back(Eigen::Triplet<int>(mesh.EF(i,1), mesh.EF(i,0),1));
            }

        Eigen::SparseMatrix<int> adjMat(mesh.F.rows(),mesh.F.rows());
        adjMat.setFromTriplets(adjTris.begin(), adjTris.end());
        Eigen::SparseMatrix<int> newAdjMat(mesh.F.rows(),mesh.F.rows()),matMult;
        igl::speye(mesh.F.rows(), mesh.F.rows(), matMult);
        for (int i=0;i<indirection;i++){
            matMult=matMult*adjMat;
            newAdjMat+=matMult;
        }

        adjMat=newAdjMat;

        std::vector<std::set<int>> ringAdjacencies(mesh.F.rows());
        for (int k=0; k<adjMat.outerSize(); ++k){
            for (Eigen::SparseMatrix<int>::InnerIterator it(adjMat,k); it; ++it){
                ringAdjacencies[it.row()].insert(it.col());
                ringAdjacencies[it.col()].insert(it.row());
            }
        }

        std::vector<int> sampleTrisVec;
        std::vector<Eigen::Vector3d> samplePointsVec;
        int nlivesamples = nsamples;
        for (int i=0;i<nsamples;i++) {
            int faceIndex = distTriangles(gen);
            int sampleIndex = -1;
            for (int j = 0; j < samplePoolAlive[faceIndex].size(); j++) {
                if (samplePoolAlive[faceIndex][j]) {
                    sampleIndex = j;
                    break;
                }
            }

            if (sampleIndex == -1) //no samples were found
                continue;

            sampleTrisVec.push_back(faceIndex);
            samplePointsVec.push_back(samplePool[faceIndex][sampleIndex]);

            //deleting samples that are less than minDist away
            //only checking faces that are "indirection" level apart
            for (std::set<int>::iterator si = ringAdjacencies[faceIndex].begin();si != ringAdjacencies[faceIndex].end(); si++) {
                int currFace = *si;
                for (int j = 0; j < samplePool[currFace].size(); j++) {
                    if (samplePoolAlive[currFace][j]) {
                        //TODO: checking if close enough and flagging "deleted" if it is
                        double dEuc2 = (samplePool[faceIndex][sampleIndex] - samplePool[currFace][j]).squaredNorm();
                        if (dEuc2 > minDist2)
                            continue;  //too far Euclideanly to delete

                        //Approximate geodesic distance
                        Eigen::Vector3d n1 = mesh.faceNormals.row(faceIndex);
                        Eigen::Vector3d n2 = mesh.faceNormals.row(currFace);
                        Eigen::Vector3d v = (samplePool[faceIndex][sampleIndex] - samplePool[currFace][j])/sqrt(dEuc2);
                        double c1 = n1.dot(v); double c2 = n2.dot(v);
                        double dGeod2;
                        if (abs(c1-c2)<10e-8)
                            dGeod2 = dEuc2/(1-c1*c1);
                        else {
                            dGeod2 = (asin(c2) - asin(c1)) / (c2 - c1);
                            dGeod2 = dGeod2*dGeod2*dEuc2;
                        }

                        if (dGeod2 > minDist2)
                            continue;  //too far geodesically to delete


                        samplePoolAlive[currFace][j] = false;
                    }

                }
            }




            /*std::queue<int> faceQueue;
            faceQueue.push(faceIndex);
            std::set<int> facesVisited;
            facesVisited.insert(faceIndex);

            do{
                int currFace = faceQueue.front();
                faceQueue.pop();
                //std::cout<<"curr face: "<<currFace<<std::endl;
                //finding out if any samples get deleted. If they are, also checking the neighbors, otherwise terminating.
                bool deleted=false;
                for (int j=0;j<samplePool[currFace].size();j++){
                    //if (samplePoolAlive[currFace][j]){
                        //TODO: checking if close enough and flagging "deleted" if it is
                        double dEuc2 = (samplePool[faceIndex][sampleIndex]-samplePool[currFace][j]).squaredNorm();
                        if (dEuc2>minDist2)
                            continue;  //too far Euclideanly to delete

                    //std::cout<<"deleting sample "<<j<<" with squared distance "<<dEuc2<<std::endl;

                        //otherwise, computing geodesic distance TODO

                        //if too close, delete sample (even if it was deleted before)
                        samplePoolAlive[currFace][j]=false;
                        deleted = true;
                    //}
                }
                if (deleted) //queuing neighboring faces since we are still close enough to the original sample
                    for (int j=0;j<3;j++)
                        if ((mesh.TT(currFace, j)!=-1)&&(facesVisited.find(mesh.TT(currFace, j))==facesVisited.end())){
                            faceQueue.push(mesh.TT(currFace,j));
                            facesVisited.insert(mesh.TT(currFace,j));
                        }


            }while (!faceQueue.empty());
            //std::cout<<"done with face"<<std::endl;
        }*/
        }

        sampleTris = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(sampleTrisVec.data(), sampleTrisVec.size());
        samplePoints.resize(samplePointsVec.size(),3);
        for (int i=0;i<samplePointsVec.size();i++)
            samplePoints.row(i)=samplePointsVec[i];

    }
}


IGL_INLINE void directional::streamlines_init(const directional::CartesianField& field,
                                              const Eigen::VectorXi& seedFaces,
                                              const double distRatio,
                                              StreamlineData &data,
                                              StreamlineState &state){
    using namespace Eigen;
    using namespace std;


    // prepare vector field
    // --------------------------

    Eigen::MatrixXd FN;
    Eigen::VectorXi order;
    Eigen::RowVectorXd sorted;
    assert(field.tb->discTangType()==discTangTypeEnum::FACE_SPACES || field.tb->discTangType()==discTangTypeEnum::VERTEX_SPACES);
    //
    if (field.tb->discTangType()==discTangTypeEnum::FACE_SPACES){
        data.stb = *((IntrinsicFaceTangentBundle*)field.tb);
        data.slMesh=data.stb.mesh;
        data.field = field;
        data.slField=data.field.extField;
    }
    if (field.tb->discTangType()==discTangTypeEnum::VERTEX_SPACES){
        IntrinsicVertexTangentBundle* vtb = (IntrinsicVertexTangentBundle*)field.tb;
        data.slMesh=vtb->mesh;

        //converting to a face-based field since streamlines don't handle fully-blown linear interpolation yet.
        data.stb.init(*(vtb->mesh));
        data.field.init(data.stb, fieldTypeEnum::RAW_FIELD, field.N);

        Eigen::MatrixXd baryCoords=Eigen::MatrixXd::Constant(data.slMesh->F.rows(),3,1.0/3.0);
        Eigen::MatrixXd interpSources, interpNormals, interpField;
        vtb->interpolate(Eigen::VectorXi::LinSpaced(data.slMesh->F.rows(),0,data.slMesh->F.rows()),
                         baryCoords,field.intField,interpSources,interpNormals,data.slField);

        //cout<<"data.slField: "<<data.slField<<endl;

        data.field.set_extrinsic_field(data.slField);
        //cout<<"data.field.intField: "<<data.field.intField<<endl;
    }
    //this is currently quite primitive and only samples in mid-faces and considers the entire triangle to be constant
    //It's fine for face-based fields but distorts vertex-based fields considerably.

    /*data.field.extField.setZero(ftb->mesh->F.rows(), field.N * 3);
    for (unsigned i = 0; i < ftb->mesh->F.rows(); ++i){
      const Eigen::RowVectorXd &n = ftb->mesh->faceNormals.row(i);
      Eigen::RowVectorXd temp(1, field.N * 3);
      temp = field.extField.row(i);
      igl::sort_vectors_ccw(temp, n, order, sorted);

      // project vectors to tangent plane
      for (int j = 0; j < field.N; ++j)
      {
        Eigen::RowVector3d pd = sorted.segment(j * 3, 3);
        pd = (pd - (n.dot(pd)) * n).normalized();
        data.field.extField.block(i, j * 3, 1, 3) = pd;
      }
    }*/

    directional::principal_matching(data.field);

    // create seeds for tracing
    // --------------------------

    int nsamples;

    if (seedFaces.rows()==0){
        assert(distRatio>=0);
        Directional::generate_sample_locations(*(data.slMesh),distRatio,data.sampleFaces,data.samplePoints);
    } else {
        data.sampleFaces=seedFaces;
        Eigen::MatrixXd BC_sample;
        igl::slice(data.slMesh->barycenters, data.sampleFaces, 1, data.samplePoints);

    }
    //nsamples = data.nsample = data.sampleFaces.size();

    // initialize state for tracing vector field
    state.currStartPoints= data.samplePoints.replicate(field.N,1);
    state.currElements = data.sampleFaces.replicate(field.N, 1);
    state.currElementTypes.resize(field.N*data.sampleFaces.size());
    for (int i=0;i<state.currElementTypes.size();i++)
        state.currElementTypes[i]=SL_FACE;

    state.currDirectionIndex.setZero(field.N*data.sampleFaces.size());
    state.segmentAlive.resize(field.N*data.sampleFaces.size());
    for (int j = 0; j < field.N; ++j) {
        for (int i = 0; i < data.sampleFaces.size(); ++i) {
            state.currDirectionIndex(j * data.sampleFaces.size() + i) = j;
            state.segmentAlive(j * data.sampleFaces.size() + i) = true;
        }
    }
    state.currTimes.setZero(field.N*data.sampleFaces.size());
    state.currSegmentIndex.setZero(field.N*data.sampleFaces.size());

    //initializing "next" values
    state.nextElements.resize(state.currElements.size());
    state.nextStartPoints.resize(state.currStartPoints.rows(),3);
    state.nextTimes.resize(state.currTimes.size());
    state.nextDirectionIndex.resize(state.currDirectionIndex.size());
    state.nextElementTypes.resize(state.nextElementTypes.size());
    for (int j = 0; j < field.N; ++j) {
        for (int i = 0; i < data.sampleFaces.size(); ++i) {
            int currIndex = j*data.sampleFaces.size()+i;
            int f0 = state.currElements(currIndex);
            int m0 = state.currDirectionIndex(currIndex);
            RowVector3d p=state.currStartPoints.row(currIndex);
            RowVector3d vec = data.slField.block(f0, 3*m0, 1,3);
            int f1, m1;
            double foundIntersection = false;
            for (int k = 0; k < 3; ++k) {
                f1 = data.slMesh->TT(state.currElements(currIndex), k);

                // edge vertices
                const Eigen::RowVector3d &q = data.slMesh->V.row(data.slMesh->F(f0, k));
                const Eigen::RowVector3d &qs = data.slMesh->V.row(data.slMesh->F(f0, (k + 1) % 3));
                // edge direction
                Eigen::RowVector3d s = qs - q;

                double u;
                double t;
                if (igl::segment_segment_intersect(p, vec, q, s, t, u, -1e-6)) {
                    foundIntersection = true;
                    state.nextElements(currIndex) = f1;
                    state.nextTimes(currIndex) = state.currTime + t;
                    state.nextStartPoints.row(currIndex) = p + t * vec;

                    // matching direction on next face
                    int e1 = data.slMesh->FE(f0, k);
                    if (data.slMesh->EF(e1, 0) == f0)
                        m1 = (data.field.matching(e1) + m0) % data.field.N;
                    else
                        m1 = (-data.field.matching(e1) + m0 + data.field.N) % data.field.N;

                    state.nextDirectionIndex(currIndex) = m1;
                    break;
                }
            }
            if (!foundIntersection) {  //something went bad, we couldn't find the next face
                state.segmentAlive(currIndex) = false;
                break;
            }

            //creating new traced segment for the new face
            state.segStart.push_back(state.currStartPoints.row(currIndex));
            state.segEnd.push_back(state.currStartPoints.row(currIndex));
            state.segNormal.push_back(data.slMesh->faceNormals.row(state.currElements(currIndex)));
            state.segOrigFace.push_back(state.currElements(currIndex));
            state.segOrigVector.push_back(j);
            state.segTimeSignatures.push_back(0.0);
            state.currSegmentIndex(currIndex)=state.segStart.size()-1;

        }

    }



    //Creating initial next values and segments

    /*state.P1 =state.start_point;
    state.P2 =state.end_point;
    state.origFace = state.current_face;
    state.origVector.resize(state.origFace.rows(), state.)*/


}

IGL_INLINE void directional::streamlines_next(const StreamlineData & data,
                                              StreamlineState & state,
                                              const double dTime){


    using namespace Eigen;
    using namespace std;



    //IntrinsicFaceTangentBundle* ftb = (IntrinsicFaceTangentBundle*)data.field.tb;  //maye not efficient since virtual lookup, or negligible?


    //Going through all ongoing streamlines. Those where the currTime+dTime < nextTime only extend their segment. Otherwise tracing forward through triangles until this happens
    for (int i = 0; i < data.field.N; ++i) {
        for (int j = 0; j < data.sampleFaces.size(); ++j) {
            int currIndex = i * data.sampleFaces.size() + j;
            int f0 = state.currElements(currIndex);
            int m0 = state.currDirectionIndex(currIndex);
            if (!state.segmentAlive(currIndex))
                continue;
            bool inCurrFace = false;
            bool keepTracing = true;
            cout<<"currIndex: "<<currIndex<<endl;
            do{
                RowVector3d vec = data.slField.block(f0, 3*m0, 1,3);
                RowVector3d p = state.currStartPoints.row(currIndex);
                if (vec.squaredNorm() < 10e-8) {   //we don't stuck in a local minima;
                    state.segmentAlive(currIndex) = false;
                    break;
                }
                cout<<"state.currTime: "<<state.currTime<<endl;
                cout<<"dTime: "<<dTime<<endl;
                cout<<"state.currTimes(currIndex): "<<state.currTimes(currIndex)<<endl;
                cout<<"state.nextTimes(currIndex): "<<state.nextTimes(currIndex)<<endl;

                if ((state.currTime+dTime>=state.currTimes(currIndex)) && (state.currTime+dTime<state.nextTimes(currIndex))) {  //updating segment within this face

                    double timeDiffFromStart =
                            state.currTime + dTime - state.currTimes(currIndex);

                    state.segEnd[state.currSegmentIndex(currIndex)]=state.segStart[state.currSegmentIndex(currIndex)]+timeDiffFromStart*vec;
                    cout<<"Stopping mid-face"<<endl;
                    break;
                } else {//trace forward
                    //finishing previous segment
                    cout<<"Tracing forward"<<endl;
                    state.segEnd[state.currSegmentIndex(currIndex)] = state.nextStartPoints.row(currIndex);
                    state.currElements(currIndex) = state.nextElements(currIndex);
                    state.currTimes(currIndex) = state.nextTimes(currIndex);
                    state.currStartPoints.row(currIndex) = state.nextStartPoints.row(currIndex);
                    state.currDirectionIndex(currIndex) = state.nextDirectionIndex(currIndex);
                    f0 = state.currElements(currIndex);
                    m0 = state.currDirectionIndex(currIndex);
                    //Updating the next element
                    int f1, m1;
                    double foundIntersection = false;
                    for (int k = 0; k < 3; ++k) {
                        f1 = data.slMesh->TT(state.currElements(currIndex), k);

                        // edge vertices
                        const Eigen::RowVector3d &q = data.slMesh->V.row(data.slMesh->F(f0, k));
                        const Eigen::RowVector3d &qs = data.slMesh->V.row(data.slMesh->F(f0, (k + 1) % 3));
                        // edge direction
                        Eigen::RowVector3d s = qs - q;

                        double u;
                        double t;
                        if (igl::segment_segment_intersect(p, vec, q, s, t, u, -1e-6)) {
                            foundIntersection = true;
                            state.nextElements(currIndex) = f1;
                            cout<<"Found intersection after: "<<t<<endl;
                            state.nextTimes(currIndex) = state.currTimes(currIndex) + t;
                            state.nextStartPoints.row(currIndex) = p + t * vec;

                            // matching direction on next face
                            int e1 = data.slMesh->FE(f0, k);
                            if (data.slMesh->EF(e1, 0) == f0)
                                m1 = (data.field.matching(e1) + m0) % data.field.N;
                            else
                                m1 = (-data.field.matching(e1) + m0 + data.field.N) % data.field.N;

                            state.nextDirectionIndex(currIndex) = m1;
                            break;
                        }
                    }
                    if (!foundIntersection) {  //something went bad, we couldn't find the next face
                         state.segmentAlive(currIndex) = false;
                         break;
                    }

                    //creating new traced segment for the new face
                    state.segStart.push_back(state.currStartPoints.row(currIndex));
                    state.segEnd.push_back(state.currStartPoints.row(currIndex));
                    state.segNormal.push_back(data.slMesh->faceNormals.row(state.currElements(currIndex)));
                    state.segOrigFace.push_back(state.currElements(currIndex));
                    state.segOrigVector.push_back(i);
                    state.segTimeSignatures.push_back(state.currTimes(currIndex));
                    state.currSegmentIndex(currIndex)=state.segStart.size()-1;
                }
            }while(keepTracing);
        }
    }
    state.currTime+=dTime;

}

