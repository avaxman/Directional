// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef MESH_FUNCTION_HEADER_FILE
#define MESH_FUNCTION_HEADER_FILE


#include "Definitions.h"
#include <iosfwd>
#include <vector>
#include <set>
#include <math.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <iostream>
#include <fstream>
#include <Eigen/Sparse>
#include <directional/FunctionMesh.h>

namespace Directional{


class Vertex{
public:
  int ID;
  Point3D Coordinates;
  EPoint3D ECoordinates;
  int AdjHalfedge;
  
  bool isFunction;
  bool Valid;
  
  Vertex():ID(-1), AdjHalfedge(-1), isFunction(false), Valid(true){}
  ~Vertex(){}
};

class Halfedge{
public:
  int ID;
  int Origin;
  int Next;
  int Prev;
  int Twin;
  int AdjFace;
  Eigen::VectorXd paramFuncs;
  vector<ENumber> exactParamFuncs;
  bool isFunction;
  bool Valid;
  
  //Parametric function values
  int OrigHalfedge;
  int OrigParamFunc;  //the original parameteric function assooicated with this edge
  //int prescribedAngleDiff;
  double prescribedAngle;  //the actual prescribed angle
  
  
  Halfedge():ID(-1), Origin(-1), Next(-1), Prev(-1), Twin(-1), AdjFace(-1), isFunction(false), Valid(true), OrigHalfedge(-1), OrigParamFunc(-1),  prescribedAngle(-1.0){}
  ~Halfedge(){}
};


class Face{
public:
  int ID;
  int AdjHalfedge;
  Vector3D Normal;
  Point3D Centroid;
  Vector3D Basis1, Basis2;
  bool Valid;
  
  Face():ID(-1), AdjHalfedge(-1), Valid(true){}
  ~Face(){}
};



class functionMesh{
public:
  vector<Vertex> Vertices;
  vector<Halfedge> Halfedges;
  vector<Face> Faces;
  
  vector<int> TransVertices;
  vector<int> InStrip;
  vector<set<int> > VertexChains;
  
  ofstream DebugLog;
  
  bool JoinFace(int heindex);
  void UnifyEdges(int heindex);
  bool CheckMesh(bool checkHalfedgeRepetition, bool CheckTwinGaps, bool checkPureBoundary);
  void CleanMesh();
  void ComputeTwins();
  void WalkBoundary(int &CurrEdge);
  void RemoveVertex(int vindex, std::deque<int>& removeVertexQueue);
  void RemoveEdge(int heindex);
  void RemoveFace(int findex, int heindex);
  void TestUnmatchedTwins();
  
  
  void GenerateMesh(const Eigen::VectorXd& funcOrientations, Mesh& HexMesh);
  
  bool SimplifyHexMesh(int N);
  
  void RemoveDegree2Faces();
  
  void Allocate(int NumofVertices, int NumofFaces, int NumofHEdges)
  {
    Vertices.resize(NumofVertices);
    Faces.resize(NumofFaces);
    Halfedges.resize(NumofHEdges);
  }
  

  //produces y = M*x
  void exactSparseMult(const Eigen::SparseMatrix<int> M, const std::vector<ENumber>& x,std::vector<ENumber>& y){
    y.resize(M.rows());
    
    for (int i=0;i<y.size();i++)
      y[i]=ENumber(0);
    
    for (int k=0; k<M.outerSize(); ++k)
      for (Eigen::SparseMatrix<int>::InnerIterator it(M,k); it; ++it)
        y[it.row()]+=ENumber((long)it.value())*x[it.col()];
  }
  
  void exactDenseMult(const Eigen::MatrixXi &nM, const Eigen::MatrixXi& dM, const std::vector<ENumber>& x, std::vector<ENumber>& y)
  {
    y.resize(nM.rows());
    for (int i=0;i<y.size();i++)
      y[i]=ENumber(0);
    for (int i=0;i<nM.rows();i++)
      for (int j=0;j<nM.cols();j++)
        y[i]+=x[j]*ENumber(nM(i,j), dM(i,j));
  }
  
  
  
  void fromHedraDCEL(const Eigen::VectorXi& D,
                     const Eigen::MatrixXd& V,
                     const Eigen::MatrixXi& F,
                     const Eigen::MatrixXi& EV,
                     const Eigen::MatrixXi& FE,
                     const Eigen::MatrixXi& EF,
                     const Eigen::MatrixXi& EFi,
                     const Eigen::MatrixXd& FEs,
                     const Eigen::VectorXi& innerEdges,
                     const Eigen::VectorXi& VH,
                     const Eigen::MatrixXi& EH,
                     const Eigen::MatrixXi& FH,
                     const Eigen::VectorXi& HV,
                     const Eigen::VectorXi& HE,
                     const Eigen::VectorXi& HF,
                     const Eigen::VectorXi& nextH,
                     const Eigen::VectorXi& prevH,
                     const Eigen::VectorXi& twinH,
                     const Eigen::MatrixXd& cutV,
                     const Eigen::MatrixXi& cutF,
                     const Eigen::MatrixXd& paramFuncsd,
                     const int d,
                     const int N,
                     const Eigen::SparseMatrix<int>& d2NMat,
                     const Eigen::SparseMatrix<int>& constraintMatInteger,
                     const Eigen::MatrixXi& symmFunc,
                     const Eigen::MatrixXd& cornerParamFuncs,
                     const Eigen::VectorXi& integerVars,
                     const Eigen::MatrixXi& embNumMat,
                     const Eigen::MatrixXi& embDenMat,
                     const Eigen::VectorXi& singVertices){
    
    using namespace std;
    using namespace CGAL;
    Vertices.resize(V.rows());
    Halfedges.resize(HE.rows());
    Faces.resize(F.rows());
    
    int NEmb = embNumMat.rows();
    
    int capN = (N%2==0 ? N/2: N);
    int capEmbN=(NEmb%2==0 ? NEmb/2 : NEmb);
    int NFull=capN+capEmbN;
    
    for (int i=0;i<V.rows();i++){
      Vertices[i].Coordinates=Point3D(V(i,0), V(i,1), V(i,2));
      Vertices[i].AdjHalfedge=VH(i);
      Vertices[i].ID=i;
    }
    
    for (int i=0;i<HE.rows();i++){
      Halfedges[i].ID=i;
      Halfedges[i].Origin=HV(i);
      Halfedges[i].Next=nextH(i);
      Halfedges[i].Prev=prevH(i);
      Halfedges[i].Twin=twinH(i);
      Halfedges[i].AdjFace=HF(i);
    }
    
    for (int i=0;i<FH.rows();i++)
      for (int j=0;j<FH.cols();j++)
        Halfedges[FH(i,j)].paramFuncs = cornerParamFuncs.block(i, NFull*j, 1, NFull).transpose();
    
    for (int i=0;i<F.rows();i++){
      Faces[i].ID=i;
      Faces[i].AdjHalfedge=FH(i);
      /*Faces[i].NumVertices=D(i);
       for (int j=0;j<D(i);j++)
       Faces[i].Vertices[j]=F(i,j);*/
    }
    
    //computing exact rational corner values by quantizing the free variables d and then manually performing the sparse matrix multiplication
    
    //resolution is set to 10e-10 of bounding box of mesh
    vector<Point3D> coordList;
    for (int i=0;i<Vertices.size();i++)
      coordList.push_back(Vertices[i].Coordinates);
    
    Bbox_3 boundBox = CGAL::bbox_3  ( coordList.begin(), coordList.end());
    
    double minRange = 3276700.0;
    for (int i=0;i<2;i++)
      minRange=std::min(minRange, boundBox.max(i)-boundBox.min(i));
    
    long Resolution=1e8; //pow(10,ceil(10/log10(minRange)));
    //cout<<"Resolution: "<<Resolution<<endl;
    
    vector<ENumber> exactParamFuncsd(paramFuncsd.size());
    for (int i=0;i<paramFuncsd.size();i++){
      exactParamFuncsd[i]=ENumber((long)round(paramFuncsd(i)*Resolution),Resolution);
      //cout<<"rounding diff of var "<<i<<" is "<<exactParamFuncsd[i].to_double()-paramFuncsd(i)<<endl;;
    }
    
    for (int i=0;i<integerVars.size();i++){
      for (int j=0;j<d;j++){
        exactParamFuncsd[d*integerVars(i)+j]=ENumber((long)round(paramFuncsd(d*integerVars(i)+j)));
        //cout<<"rounding diff of integer var "<<d*integerVars(i)+j<<" is "<<exactParamFuncsd[d*integerVars(i)+j].to_double()-paramFuncsd(d*integerVars(i)+j)<<endl;
      }
    }
    
    vector<ENumber> exactParamFuncsVec;
    exactSparseMult(d2NMat, exactParamFuncsd,exactParamFuncsVec);
    
    vector<ENumber> constraintError;
    exactSparseMult(constraintMatInteger, exactParamFuncsd,constraintError);
    
    ENumber MaxError(0);
    
    for (int i=0;i<constraintError.size();i++)
      if (abs(constraintError[i])>MaxError)
        MaxError=abs(constraintError[i]);
    
    cout<<"constraintMatInteger*exactParamFuncsd MaxError: "<<MaxError<<endl;
    
    //introducing offset
    Eigen::VectorXi offset= symmFunc.rowwise().sum();
    int offDen=(N%3==0 ? 3 : 2);
    if (N%6==0) {offDen=1; offset.setZero();}
    
    
    //ONLY FOR EXAMPLE!!!
    //offDen=1; offset.setZero();
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    vector<vector<ENumber> > exactParamFuncsN(cutV.rows());
    for(int i = 0; i < cutV.rows(); i++){
      exactParamFuncsN[i].resize(N);
      for (int j=0;j<N;j++)
        exactParamFuncsN[i][j] =exactParamFuncsVec[N * i+j]+ENumber(offset(j),offDen);
    }
    
    //allocating per corner
    vector<vector<ENumber> > exactWholeCornerParamFuncsN(F.rows());
    
    
    //assuming the functions are ordered where the sign symmetry follows in the later half of the functions

    for (int i=0;i<F.rows();i++){
      exactWholeCornerParamFuncsN[i].resize((capN+capEmbN)*3);
      for (int j=0;j<3;j++)
        for (int k=0;k<capN;k++)
          exactWholeCornerParamFuncsN[i][(capN+capEmbN)*j+k] = exactParamFuncsN[cutF(i,j)][k];
    }
    
    //putting in embedded functions
   
    for (int i=0;i<F.rows();i++){
      for (int j=0;j<3;j++){
        vector<ENumber> origCornerFunc(N);
        for (int k=0;k<N;k++)
          origCornerFunc[k]=exactParamFuncsN[cutF(i,j)][k];
        
        /*cout<<"origCornerFunc: ";
        for (int k=0;k<N;k++)
          cout<<origCornerFunc[k].to_double()<<" ";
        cout<<endl;*/
        
        vector<ENumber> embCornerFunc;
        exactDenseMult(embNumMat, embDenMat, origCornerFunc, embCornerFunc);
        /*cout<<"embCornerFunc: ";
        for (int k=0;k<NEmb;k++)
          cout<<embCornerFunc[k].to_double()<<" ";
        cout<<endl;*/
        for (int k=capN;k<capN+capEmbN;k++)
          exactWholeCornerParamFuncsN[i][NFull*j+k] = embCornerFunc[k-capN];
        
        /*cout<<"exactWholeCornerParamFuncsN[i]: ";
        for (int k=0;k<NFull;k++)
          cout<<exactWholeCornerParamFuncsN[i][NFull*j+k].to_double()<<" ";
        cout<<endl;*/
        
        //cout<<"paramCornerFunc: "<<Halfedges[FH(i,j)].paramFuncs<<endl;
      }
    }
    
    for (int i=0;i<FH.rows();i++)
      for (int j=0;j<FH.cols();j++){
        Halfedges[FH(i,j)].exactParamFuncs.resize(NFull);
        for (int k=0;k<NFull;k++)
          Halfedges[FH(i,j)].exactParamFuncs[k] = exactWholeCornerParamFuncsN[i][NFull*j+k];
      }
    
    //sanity check
    double maxError = -32767000.0;
    for (int i=0;i<Halfedges.size();i++){
      for (int j=0;j<NFull;j++){
        double fromExact = Halfedges[i].exactParamFuncs[j].to_double();
        //cout<<"fromExact: "<<fromExact<<endl;
        //cout<<"Halfedges[i].paramFuncs[j]: "<<Halfedges[i].paramFuncs[j]<<endl;
        if (abs(fromExact-Halfedges[i].paramFuncs[j])>maxError)
          maxError =abs(fromExact-Halfedges[i].paramFuncs[j]);
      }
      
    }
    cout<<"double from exact in halfedges maxError: "<<maxError<<endl;
  }
  
  
  //corner angles is per vertex in each F
  void toHedra(Eigen::MatrixXd& generatedV, Eigen::VectorXi& generatedD, Eigen::MatrixXi& generatedF, Eigen::MatrixXi& generatedFfuncNum, Eigen::MatrixXd& cornerAngles){
    generatedV.resize(Vertices.size(),3);
    
    generatedD.resize(Faces.size());
    
    for (int i=0;i<Vertices.size();i++)
      generatedV.row(i)<<Vertices[i].Coordinates.x(), Vertices[i].Coordinates.y(),Vertices[i].Coordinates.z();
    
    for (int i=0;i<Faces.size();i++){
      int hebegin = Faces[i].AdjHalfedge;
      //reseting to first vertex
      int vCount=0;
      int heiterate=hebegin;
      do{
        vCount++;
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
      generatedD(i)=vCount;
    }
    
    generatedF.resize(Faces.size(),generatedD.maxCoeff());
    for (int i=0;i<Faces.size();i++){
      int hebegin = Faces[i].AdjHalfedge;
      int vCount=0;
      int heiterate=hebegin;
      do{
        generatedF(i,vCount++)=Halfedges[heiterate].Origin;
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
      
    }
    
    generatedFfuncNum.resize(Faces.size(),generatedD.maxCoeff());
    cornerAngles=Eigen::MatrixXd::Constant(Faces.size(),generatedD.maxCoeff(),-1.0);
    //prescribedAnglesInt.resize(Faces.size(),generatedD.maxCoeff());
    for (int i=0;i<Faces.size();i++){
      int hebegin = Faces[i].AdjHalfedge;
      int vCount=0;
      int heiterate=hebegin;
      do{
        generatedFfuncNum(i,vCount)=Halfedges[heiterate].OrigParamFunc;
        cornerAngles(i,vCount++)=Halfedges[heiterate].prescribedAngle;
        //prescribedAnglesInt(i,vCount++)=Halfedges[heiterate].prescribedAngleDiff;
        //cout<<"Halfedges[heiterate].prescribedAngleDiff: "<<Halfedges[heiterate].prescribedAngleDiff<<endl;
        heiterate=Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
    }
    
    
  }
  
  Mesh(){}
  ~Mesh(){}
  
};


struct EdgeData{
  int ID;
  bool isFunction;
  int OrigHalfedge;
  bool isBoundary;
  int funcNum;  //in case of function segment
  
  EdgeData():ID(-1), isFunction(false), OrigHalfedge(-1), isBoundary(false), funcNum(-1){}
  ~EdgeData(){};
};

namespace CGAL {

template <class ArrangementA, class ArrangementB, class ArrangementR>
class Arr_function_overlay_traits
{
public:
  
  typedef typename ArrangementA::Face_const_handle    Face_handle_A;
  typedef typename ArrangementB::Face_const_handle    Face_handle_B;
  typedef typename ArrangementR::Face_handle          Face_handle_R;
  
  typedef typename ArrangementA::Vertex_const_handle    Vertex_handle_A;
  typedef typename ArrangementB::Vertex_const_handle    Vertex_handle_B;
  typedef typename ArrangementR::Vertex_handle          Vertex_handle_R;
  
  typedef typename ArrangementA::Halfedge_const_handle    Halfedge_handle_A;
  typedef typename ArrangementB::Halfedge_const_handle    Halfedge_handle_B;
  typedef typename ArrangementR::Halfedge_handle          Halfedge_handle_R;
  
public:
  
  virtual void create_face (Face_handle_A f1,
                            Face_handle_B f2,
                            Face_handle_R f) const
  {
    // Overlay the data objects associated with f1 and f2 and store the result
    // with f.
    f->set_data (f1->data());
    return;
  }
  
  //-1 - triangle vertex (non-hex vertex)
  //-2 - hex vertex
  virtual void	create_vertex ( Vertex_handle_A v1, Vertex_handle_B v2, Vertex_handle_R v)
  {
    v->set_data(-2);
  }
  
  
  virtual void create_vertex ( Vertex_handle_A v1, Halfedge_handle_B e2, Vertex_handle_R v)
  {
    v->set_data(-1);
  }
  
  
  virtual void create_vertex ( Vertex_handle_A v1, Face_handle_B f2, Vertex_handle_R v)
  {
    v->set_data(-1);
  }
  
  virtual void create_vertex ( Halfedge_handle_A e1, Vertex_handle_B v2, Vertex_handle_R v)
  {
    v->set_data(-2);
  }
  
  virtual void create_vertex ( Face_handle_A f1, Vertex_handle_B v2, Vertex_handle_R v)
  {
    v->set_data(-2);
  }
  
  virtual void create_vertex ( Halfedge_handle_A e1, Halfedge_handle_B e2, Vertex_handle_R v)
  {
    v->set_data(-1);
  }
  
  
  virtual void create_edge ( Halfedge_handle_A e1, Halfedge_handle_B e2, Halfedge_handle_R e)
  {
    EdgeData data;
    data.ID=-2;
    data.OrigHalfedge=e1->data().OrigHalfedge;
    data.funcNum = e2->data().funcNum;
    //std::cout<<"e2->data().funcNum: "<<e2->data().funcNum<<std::endl;
    e->set_data(data);
    e->twin()->set_data(data);
  }
  
  
  virtual void create_edge ( Halfedge_handle_A e1, Face_handle_B f2, Halfedge_handle_R e)
  {
    EdgeData data;
    data.ID=-1;
    data.OrigHalfedge=e1->data().OrigHalfedge;
    data.funcNum = -1;  //triangle edge
    //std::cout<<"data.funcNum: "<<data.funcNum<<std::endl;
    e->set_data(data);
    e->twin()->set_data(data);
  }
  
  virtual void create_edge ( Face_handle_A f1, Halfedge_handle_B e2, Halfedge_handle_R e)
  {
    EdgeData data;
    data.ID=-2;
    data.funcNum = e2->data().funcNum;
    // if (f1->data())
    //std::cout<<"e2->data().funcNum: "<<e2->data().funcNum<<std::endl;
    e->set_data(data);
    e->twin()->set_data(data);
  }
  
};




} //namespace CGAL


bool mesh_function(const Eigen::MatrixXd& V,
                   const Eigen::MatrixXi& F,
                   const Eigen::MatrixXi& EV,
                   const Eigen::MatrixXi& EF,
                   const Eigen::MatrixXi& FE,
                   const IntegrationData& intData,
                   const Eigen::MatrixXd NFunction,
                   const bool verbose,
                   Eigen::MatrixXd& V,
                   Eigen::VectorXi& D,
                   Eigen::MatrixXi& F){
  
  
  int N = NFunction.cols();
  FunctionMesh TMesh, FMesh;
  
  Eigen::VectorXi VHPoly, HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, HVPoly,innerEdgesPoly;
  Eigen::MatrixXi EHPoly,EFiPoly, FHPoly, EFPoly,EVPoly,FEPoly;
  Eigen::MatrixXd FEsPoly;
  hedra::polygonal_edge_topology(VectorXi::Constant(FMeshWhole.rows(),3), FMeshWhole,EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly);
  hedra::dcel(VectorXi::Constant(FMeshWhole.rows(),3),FMeshWhole,EVPoly,EFPoly, EFiPoly,innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly);
  
  TMesh.fromHedraDCEL(VectorXi::Constant(FMeshWhole.rows(),3),VMeshWhole, FMeshWhole, EVPoly,FEPoly,EFPoly, EFiPoly, FEsPoly, innerEdgesPoly,VHPoly, EHPoly, FHPoly,  HVPoly,  HEPoly, HFPoly, nextHPoly, prevHPoly, twinHPoly, VMeshCut, FMeshCut, paramFuncsd, pd.d, pd.N, intData.vertexTrans2CutMatInteger*pd.symmMatInteger*intData.singIntSpanMatInteger*intData.intSpanMatInteger,  intData.constraintMatInteger*intData.symmMatInteger*intData.singIntSpanMatInteger*intData.intSpanMatInteger, intData.symmFunc*intData.intFunc, reducedCornerFuncs/*wholeCornerParamFuncsN*/, pd.integerVars, embNumMat, embDenMat, singVertices);
  
  if (verbose);
  cout<<"Generating mesh"<<endl;
  TMesh.GenerateMesh(FMesh);
  cout<<"Done generating!"<<endl;
  
  Eigen::VectorXi genInnerEdges,genTF;
  Eigen::MatrixXi genEV,genEFi, genEF,genFE, genTEdges;
  Eigen::MatrixXd genFEs, genCEdges, genVEdges;
  
  bool success = FMesh.SimplifyMesh(N);
  
  FMesh.toHedra(V,D,  F, simpFfuncNum);

}

} //namespace directional






#endif
