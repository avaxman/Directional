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
        std::discrete_distribution<double> distTriangles(faceAreasVec.begin(), faceAreasVec.end());
        std::uniform_real_distribution<double> distBarycentrics(0.0,1.0);
        std::map<int, double> map;
        int nsamples = (int)(10/distRatio)*mesh.F.rows();

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
        std::vector<int> sampleTrisVec;
        std::vector<Eigen::Vector3d> samplePointsVec;
        int nlivesamples = nsamples;
        for (int i=0;i<nsamples;i++) {
            int faceIndex = distTriangles(gen);
            int sampleIndex=-1;
            for (int j = 0; j < samplePoolAlive[faceIndex].size(); j++) {
                if (samplePoolAlive[faceIndex][j]) {
                    sampleIndex = j;
                    break;
                }
            }

            if (sampleIndex==-1) //no samples were found
                continue;

            sampleTrisVec.push_back(faceIndex);
            samplePointsVec.push_back(samplePool[faceIndex][sampleIndex]);
            //deleting samples that are less than minDist away
            std::queue<int> faceQueue;
            faceQueue.push(faceIndex);
            std::set<int> facesVisited;
            facesVisited.insert(faceIndex);

            do{
                int currFace = faceQueue.front();
                faceQueue.pop();
                //std::cout<<"curr face: "<<currFace<<std::endl;
                //finding out if any samples get deleted. If they are, also checking the neighbors, otherwise terminating.
                bool deleted=true;
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

    state.numSteps=0;

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
    nsamples = data.nsample = data.sampleFaces.size();

    // initialize state for tracing vector field

    state.start_point = data.samplePoints.replicate(field.N,1);
    state.end_point = state.start_point;

    state.current_face = data.sampleFaces.replicate(1, field.N);

    /*state.P1 =state.start_point;
    state.P2 =state.end_point;
    state.origFace = state.current_face;
    state.origVector.resize(state.origFace.rows(), state.)*/

    state.current_direction.setZero(nsamples, field.N);
    for (int i = 0; i < nsamples; ++i)
        for (int j = 0; j < field.N; ++j)
            state.current_direction(i, j) = j;

}

IGL_INLINE void directional::streamlines_next(const StreamlineData & data,
                                              StreamlineState & state){


    using namespace Eigen;
    using namespace std;

    int nsample = data.nsample;

    //IntrinsicFaceTangentBundle* ftb = (IntrinsicFaceTangentBundle*)data.field.tb;  //maye not efficient since virtual lookup, or negligible?
    state.start_point = state.end_point;
    state.numSteps++;

    Eigen::VectorXi currFace=Eigen::VectorXi::Zero(state.start_point.rows());
    Eigen::VectorXi currVector=Eigen::VectorXi::Zero(state.start_point.rows());
    Eigen::VectorXi currTime=Eigen::VectorXi::Zero(state.start_point.rows());

    for (int i = 0; i < data.field.N; ++i)
    {
        for (int j = 0; j < nsample; ++j)
        {
            currFace(j + nsample * i)=data.sampleFaces(j);
            currVector(j + nsample * i)=i;
            currTime(j + nsample * i)=state.numSteps;

            int f0 = state.current_face(j,i);
            if (f0 == -1) // reach boundary
                continue;
            int m0 = state.current_direction(j, i);

            // the starting point of the vector
            const Eigen::RowVector3d &p = state.start_point.row(j + nsample * i);

            // the direction where we are trying to go
            const Eigen::RowVector3d &r = data.slField.block(f0, 3 * m0, 1, 3);


            // new state,
            int f1, m1;

            for (int k = 0; k < 3; ++k)
            {
                f1 = data.slMesh->TT(f0, k);

                // edge vertices
                const Eigen::RowVector3d &q = data.slMesh->V.row(data.slMesh->F(f0, k));
                const Eigen::RowVector3d &qs = data.slMesh->V.row(data.slMesh->F(f0, (k + 1) % 3));
                // edge direction
                Eigen::RowVector3d s = qs - q;

                double u;
                double t;
                if (igl::segment_segment_intersect(p, r, q, s, t, u, -1e-6))
                {
                    // point on next face
                    state.end_point.row(j + nsample * i) = p + t * r;
                    state.current_face(j,i) = f1;

                    // matching direction on next face
                    int e1 = data.slMesh->FE(f0, k);
                    if (data.slMesh->EF(e1, 0) == f0)
                        m1 = (data.field.matching(e1)+m0)%data.field.N; //m1 = data.match_ab(e1, m0);
                    else
                        m1 = (-data.field.matching(e1)+m0+data.field.N)%data.field.N;  //data.match_ba(e1, m0);

                    state.current_direction(j, i) = m1;
                    break;
                }

            }
        }
    }

    //aggregating
    state.P1.conservativeResize(state.P1.rows()+state.start_point.rows(),3);
    state.P2.conservativeResize(state.P2.rows()+state.end_point.rows(),3);
    state.P1.block(state.P1.rows()-state.start_point.rows(),0,state.start_point.rows(),3)=state.start_point;
    state.P2.block(state.P2.rows()-state.end_point.rows(),0,state.end_point.rows(),3)=state.end_point;
    state.origFace.conservativeResize(state.origFace.size()+state.start_point.rows());
    state.origVector.conservativeResize(state.origVector.size()+state.start_point.rows());
    state.timeSignature.conservativeResize(state.timeSignature.size()+state.start_point.rows());
    state.origFace.tail(state.start_point.rows())=currFace;
    state.origVector.tail(state.end_point.rows())=currVector;
    state.timeSignature.tail(state.end_point.rows())=currTime;

}
