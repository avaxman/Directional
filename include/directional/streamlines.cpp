// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>, 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Geometry>
#include <igl/edge_topology.h>
#include <igl/sort_vectors_ccw.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/barycenter.h>
#include <igl/slice.h>
#include <igl/speye.h>
#include <directional/principal_matching.h>
#include <directional/streamlines.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/definitions.h>


namespace Directional {
IGL_INLINE void generate_sample_locations(const Eigen::MatrixXi& F,
                                          const Eigen::MatrixXi& EF,
                                          const int ringDistance,
                                          Eigen::VectorXi& samples)
{
  //creating adjacency matrix
  std::vector<Eigen::Triplet<int>> adjTris;
  for (int i=0;i<EF.rows();i++)
    if ((EF(i,0)!=-1)&&(EF(i,1)!=-1)){
      adjTris.push_back(Eigen::Triplet<int>(EF(i,0), EF(i,1),1));
      adjTris.push_back(Eigen::Triplet<int>(EF(i,1), EF(i,0),1));
    }
  
  Eigen::SparseMatrix<int> adjMat(F.rows(),F.rows());
  adjMat.setFromTriplets(adjTris.begin(), adjTris.end());
  Eigen::SparseMatrix<int> newAdjMat(F.rows(),F.rows()),matMult;
  igl::speye(F.rows(), F.rows(), matMult);
  for (int i=0;i<ringDistance;i++){
    matMult=matMult*adjMat;
    newAdjMat+=matMult;
  }
  
  //cout<<"newAdjMat: "<<newAdjMat<<endl;
  
  adjMat=newAdjMat;
  
  std::vector<std::set<int>> ringAdjacencies(F.rows());
  for (int k=0; k<adjMat.outerSize(); ++k){
    for (Eigen::SparseMatrix<int>::InnerIterator it(adjMat,k); it; ++it){
      ringAdjacencies[it.row()].insert(it.col());
      ringAdjacencies[it.col()].insert(it.row());
    }
  }

  Eigen::VectorXi sampleMask=Eigen::VectorXi::Zero(F.rows());
  for (int i=0;i<F.rows();i++){
    if (sampleMask(i)!=0) //occupied face
      continue;
    
    sampleMask(i)=2;
    //clearing out all other faces
    //cout<<"closeby set to face "<<i<<endl;
    for (std::set<int>::iterator si=ringAdjacencies[i].begin();si!=ringAdjacencies[i].end();si++){
      if (sampleMask(*si)==0)
        sampleMask(*si)=1;
      //cout<<*si<<" ";
    }
    //cout<<endl;
    
  }
  
  std::vector<int> samplesList;
  for (int i=0;i<sampleMask.size();i++)
    if (sampleMask(i)==2)
      samplesList.push_back(i);
  
  samples = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(samplesList.data(), samplesList.size());
}
}


IGL_INLINE void directional::streamlines_init(const directional::CartesianField& field,
                                              const Eigen::VectorXi& seedLocations,
                                              const int ringDistance,
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
  assert(field.tb->discTangType()==discTangTypeEnum::FACE_SPACES);
  IntrinsicFaceTangentBundle* ftb = (IntrinsicFaceTangentBundle*)field.tb;
  data.field = field;
  data.field.extField.setZero(ftb->mesh->F.rows(), field.N * 3);
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
  }
  Eigen::VectorXd effort;
  Eigen::VectorXi singIndices, singVertices;
  directional::principal_matching(data.field);
  
  // create seeds for tracing
  // --------------------------
 
  int nsamples;
  
  if (seedLocations.rows()==0){
    assert(ringDistance>=0);
    Directional::generate_sample_locations(ftb->mesh->F,ftb->mesh->EF,ringDistance,data.samples);
    nsamples = data.nsample = data.samples.size();
    nsamples = data.nsample;
    /*Eigen::VectorXd r;
    r.setRandom(nsamples, 1);
    r = (1 + r.array()) / 2.;
    samples = (r.array() * F.rows()).cast<int>();
    data.nsample = nsamples;*/
  } else {
    data.samples=seedLocations;
    nsamples = data.nsample = seedLocations.size();
  }
  
  Eigen::MatrixXd BC_sample;
  igl::slice(ftb->mesh->barycenters, data.samples, 1, BC_sample);
  
  // initialize state for tracing vector field
  
  state.start_point = BC_sample.replicate(field.N,1);
  state.end_point = state.start_point;
  
  state.current_face = data.samples.replicate(1, field.N);
  
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

  IntrinsicFaceTangentBundle* ftb = (IntrinsicFaceTangentBundle*)data.field.tb;  //maye not efficient since virtual lookup, or negligible?
  state.start_point = state.end_point;
  state.numSteps++;
  
  Eigen::VectorXi currFace=Eigen::VectorXi::Zero(state.start_point.rows());
  Eigen::VectorXi currVector=Eigen::VectorXi::Zero(state.start_point.rows());
  Eigen::VectorXi currTime=Eigen::VectorXi::Zero(state.start_point.rows());
  
  for (int i = 0; i < data.field.N; ++i)
  {
    for (int j = 0; j < nsample; ++j)
    {
      currFace(j + nsample * i)=data.samples(j);
      currVector(j + nsample * i)=i;
      currTime(j + nsample * i)=state.numSteps;
      
      int f0 = state.current_face(j,i);
      if (f0 == -1) // reach boundary
        continue;
      int m0 = state.current_direction(j, i);
      
      // the starting point of the vector
      const Eigen::RowVector3d &p = state.start_point.row(j + nsample * i);
      
      // the direction where we are trying to go
      const Eigen::RowVector3d &r = data.field.extField.block(f0, 3 * m0, 1, 3);
      
      
      // new state,
      int f1, m1;
      
      for (int k = 0; k < 3; ++k)
      {
        f1 = ftb->mesh->TT(f0, k);
        
        // edge vertices
        const Eigen::RowVector3d &q = ftb->mesh->V.row(ftb->mesh->F(f0, k));
        const Eigen::RowVector3d &qs = ftb->mesh->V.row(ftb->mesh->F(f0, (k + 1) % 3));
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
          int e1 = ftb->mesh->FE(f0, k);
          if (ftb->mesh->EF(e1, 0) == f0)
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
