// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef FUNCTION_MESH_CLASS_HEADER_FILE
#define FUNCTION_MESH_CLASS_HEADER_FILE


#include <iosfwd>
#include <vector>
#include <set>
#include <math.h>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <iostream>
#include <fstream>
#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_2.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <igl/PI.h>
#include <directional/FunctionMesh.h>

namespace Directional{

class NFunctionMesher{
public:
  
  
  struct MergeData {
    const bool operator()(const int& v1, const int& v2) const {return v1; }
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
  
  //typedef CGAL::Filtered_kernel< CGAL::Cartesian<CORE::Expr> > CKernel;
   typedef ::CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef ::CGAL::Simple_cartesian<::CGAL::Gmpq>  EKernel;
  typedef ::CGAL::Gmpz  EInt;
  

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
    virtual void  create_vertex ( Vertex_handle_A v1, Vertex_handle_B v2, Vertex_handle_R v)
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

   typedef Kernel::RT Number;
   typedef EKernel::RT ENumber;
   typedef Kernel::Point_2 Point2D;
   typedef Kernel::Point_3 Point3D;
   typedef Kernel::Ray_2 Ray2D;
   typedef Kernel::Line_2 Line2D;
   typedef EKernel::Line_2 ELine2D;
   typedef Kernel::Line_3 Line3D;
   typedef Kernel::Vector_2 Vector2D;
   typedef Kernel::Vector_3 Vector3D;
   typedef Kernel::Segment_2 Segment2D;
   typedef EKernel::Segment_2 ESegment2D;
   typedef Kernel::Segment_3 Segment3D;
   typedef Kernel::Plane_3 Plane3D;
   typedef Kernel::Triangle_3 Triangle3D;
   typedef Kernel::Triangle_2 Triangle2D;
   typedef Kernel::Circle_2 Circle2D;
   typedef EKernel::Point_2 EPoint2D;
   typedef EKernel::Vector_2 EVector2D;
   typedef EKernel::Point_3 EPoint3D;
   typedef EKernel::Triangle_2 ETriangle2D;
   typedef Kernel::Ray_3 Ray3D;
   typedef Kernel::Direction_2 Direction2D;
   typedef EKernel::Direction_2 EDirection2D;
   typedef Kernel::Direction_3 Direction3D;
   typedef Kernel::Sphere_3 Sphere3D;
   typedef Kernel::Aff_transformation_3 Transform3D;
   typedef ::CGAL::Polyhedron_3<Kernel> Polyhedron3D;
   typedef ::CGAL::Polygon_2<Kernel>    Polygon2D;
   typedef ::CGAL::Arr_linear_traits_2<EKernel>                    Traits_2;
   typedef Traits_2::Point_2                                     Point_2;
   typedef Traits_2::Segment_2                                   Segment_2;
   typedef Traits_2::Ray_2                                       Ray_2;
   typedef Traits_2::Line_2                                      Line_2;
   typedef Traits_2::X_monotone_curve_2                          X_monotone_curve_2;
   typedef ::CGAL::Arr_extended_dcel<Traits_2, int,EdgeData,int>   Dcel;
   typedef ::CGAL::Arrangement_with_history_2<Traits_2, Dcel>      Arr_2;
   typedef Arr_2::Face_iterator                                  Face_iterator;
   typedef Arr_2::Face_handle                                    Face_handle;
   typedef Arr_2::Edge_iterator                                  Edge_iterator;
   typedef Arr_2::Halfedge_iterator                              Halfedge_iterator;
   typedef Arr_2::Vertex_iterator                                Vertex_iterator;
   typedef Arr_2::Vertex_handle                                  Vertex_handle;
   typedef Arr_2::Halfedge_handle                                Halfedge_handle;
   typedef Arr_2::Ccb_halfedge_circulator                        Ccb_halfedge_circulator;
   typedef Arr_function_overlay_traits <Arr_2,Arr_2,Arr_2>      Overlay_traits;
   
   typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::undirectedS> Graph;
  

  
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
    std::vector<ENumber> exactParamFuncs;
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

  std::vector<Vertex> Vertices;
  std::vector<Halfedge> Halfedges;
  std::vector<Face> Faces;
  
  std::vector<int> TransVertices;
  std::vector<int> InStrip;
  std::vector<std::set<int> > VertexChains;
  
 
  bool JoinFace(int heindex){
    if (Halfedges[heindex].Twin<0)
      return true;  //there is no joining of boundary faces
    
    int Face1=Halfedges[heindex].AdjFace;
    int Face2=Halfedges[Halfedges[heindex].Twin].AdjFace;
  
    int hebegin=Faces[Face1].AdjHalfedge;
    int heiterate = hebegin;
    do{
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
   
    hebegin=Faces[Face1].AdjHalfedge;
    heiterate = hebegin;
    do{
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
  
    //check if spike edge
    if ((Halfedges[heindex].Prev==Halfedges[heindex].Twin)||(Halfedges[heindex].Next==Halfedges[heindex].Twin)){
      
      
      int CloseEdge=heindex;
      if (Halfedges[heindex].Prev==Halfedges[heindex].Twin)
        CloseEdge=Halfedges[heindex].Twin;
      
      Halfedges[CloseEdge].Valid=Halfedges[Halfedges[CloseEdge].Twin].Valid=false;
      
      Vertices[Halfedges[CloseEdge].Origin].AdjHalfedge=Halfedges[Halfedges[CloseEdge].Twin].Next;
      Faces[Face1].AdjHalfedge=Halfedges[CloseEdge].Prev;
      
      Halfedges[Halfedges[CloseEdge].Prev].Next=Halfedges[Halfedges[CloseEdge].Twin].Next;
      Halfedges[Halfedges[Halfedges[CloseEdge].Twin].Next].Prev=Halfedges[CloseEdge].Prev;
      
      Vertices[Halfedges[Halfedges[CloseEdge].Twin].Origin].Valid=false;
      //Faces[Face1].NumVertices-=2;  //although one vertex should appear twice
      
     
      int hebegin=Faces[Face1].AdjHalfedge;
      int heiterate = hebegin;
      do{
        
        heiterate = Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);

      hebegin=Faces[Face1].AdjHalfedge;
      heiterate = hebegin;
      do{
        heiterate = Halfedges[heiterate].Next;
      }while (heiterate!=hebegin);
      
      
      return true;
    }
    
    if (Face1==Face2)
      return false;  //we don't remove non-spike edges with the same faces to not disconnect a chain
    
    hebegin=Faces[Face2].AdjHalfedge;
    heiterate = hebegin;
    do{
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
    
    Faces[Face1].AdjHalfedge=Halfedges[heindex].Next;
    Faces[Face2].Valid=false;
    
    //Faces[Face2].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
    
    Halfedges[heindex].Valid=Halfedges[Halfedges[heindex].Twin].Valid=false;
    
    Halfedges[Halfedges[heindex].Next].Prev=Halfedges[Halfedges[heindex].Twin].Prev;
    Halfedges[Halfedges[Halfedges[heindex].Twin].Prev].Next=Halfedges[heindex].Next;
    
    Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[heindex].Prev;
    Halfedges[Halfedges[heindex].Prev].Next=Halfedges[Halfedges[heindex].Twin].Next;
    
    Vertices[Halfedges[heindex].Origin].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
    Vertices[Halfedges[Halfedges[heindex].Next].Origin].AdjHalfedge=Halfedges[heindex].Next;
    
    //all other floating halfedges should renounce this one
    for (int i=0;i<Halfedges.size();i++)
      if (Halfedges[i].AdjFace==Face2)
        Halfedges[i].AdjFace=Face1;
    
    //Faces[Face1].NumVertices+=Faces[Face2].NumVertices-2;
    
    //DebugLog<<"Official number of vertices: "<<Faces[Face1].NumVertices;
    
    hebegin=Faces[Face1].AdjHalfedge;
    heiterate = hebegin;
    int currVertex=0;
    do{
      //Faces[Face1].Vertices[currVertex++]=Halfedges[heiterate].Origin;
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);

    return true;
    
    
  }
  void UnifyEdges(int heindex){
    //if (Halfedges[heindex].Twin<0)
     //  return;
     //adjusting source
     Vertices[Halfedges[heindex].Origin].Valid=false;
     Halfedges[heindex].Origin=Halfedges[Halfedges[heindex].Prev].Origin;
     if (Halfedges[heindex].prescribedAngle<0.0)
       Halfedges[heindex].prescribedAngle=Halfedges[Halfedges[heindex].Prev].prescribedAngle;
     Vertices[Halfedges[heindex].Origin].AdjHalfedge=heindex;
     
     Faces[Halfedges[heindex].AdjFace].AdjHalfedge=Halfedges[heindex].Next;
     //Faces[Halfedges[heindex].AdjFace].NumVertices--;
     
     
     
     //adjusting halfedges
     Halfedges[Halfedges[heindex].Prev].Valid=false;
     Halfedges[heindex].Prev=Halfedges[Halfedges[heindex].Prev].Prev;
     Halfedges[Halfedges[heindex].Prev].Next=heindex;
    
     //adjusting twin, if exists
     if (Halfedges[heindex].Twin>=0){
       if (Halfedges[Halfedges[heindex].Twin].prescribedAngle<0.0)
         Halfedges[Halfedges[heindex].Twin].prescribedAngle=Halfedges[Halfedges[Halfedges[heindex].Twin].Next].prescribedAngle;
       Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Valid=false;
       Halfedges[Halfedges[heindex].Twin].Next=Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Next;
       Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[heindex].Twin;
       Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
       //Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices--;
     }
  }
  bool CheckMesh(const bool verbose, const bool checkHalfedgeRepetition, const bool CheckTwinGaps, const bool checkPureBoundary){
    for (int i=0;i<Vertices.size();i++){
       if (!Vertices[i].Valid)
         continue;
       
       if (Vertices[i].AdjHalfedge==-1){
         if (verbose) std::cout<<"Valid Vertex "<<i<<" points to non-valid value -1 "<<std::endl;
         return false;
       }
       
       if (!Halfedges[Vertices[i].AdjHalfedge].Valid){
         if (verbose) std::cout<<"Valid Vertex "<<i<<" points to non-valid AdjHalfedge "<<Vertices[i].AdjHalfedge<<std::endl;
         return false;
       }
       
       
       if (Halfedges[Vertices[i].AdjHalfedge].Origin!=i){
         if (verbose) std::cout<<"Adjacent Halfedge "<<Vertices[i].AdjHalfedge<<" of vertex "<<i<<"does not point back"<<std::endl;
         return false;
       }
       
     }
     
     for (int i=0;i<Halfedges.size();i++){
       if (!Halfedges[i].Valid)
         continue;
       
       
       if (Halfedges[i].Next==-1){
         if (verbose) std::cout<<"Valid Halfedge "<<i<<"points to Next non-valid value -1"<<std::endl;
         return false;
       }
       
       if (Halfedges[i].Prev==-1){
         if (verbose) std::cout<<"Valid Halfedge "<<i<<"points to Prev non-valid value -1"<<std::endl;
         return false;
       }
       
       
       if (Halfedges[i].Origin==-1){
         if (verbose) std::cout<<"Valid Halfedge "<<i<<"points to Origin non-valid value -1"<<std::endl;
         return false;
       }
       
       if (Halfedges[i].AdjFace==-1){
         if (verbose) std::cout<<"Valid Halfedge "<<i<<"points to AdjFace non-valid value -1"<<std::endl;
         return false;
       }
       
       if (Halfedges[Halfedges[i].Next].Prev!=i){
         if (verbose) std::cout<<"Halfedge "<<i<<"Next is "<<Halfedges[i].Next<<" which doesn't point back as Prev"<<std::endl;
         return false;
       }
       
       
       if (Halfedges[Halfedges[i].Prev].Next!=i){
         if (verbose) std::cout<<"Halfedge "<<i<<"Prev is "<<Halfedges[i].Prev<<" which doesn't point back as Next"<<std::endl;
         return false;
       }
       
       if (!Vertices[Halfedges[i].Origin].Valid){
         if (verbose) std::cout<<"The Origin of halfedges "<<i<<", vertex "<<Halfedges[i].Origin<<" is not valid"<<std::endl;
         return false;
       }
       
       if (!Faces[Halfedges[i].AdjFace].Valid){
        if (verbose) std::cout<<"The face of halfedges "<<i<<", face "<<Halfedges[i].AdjFace<<" is not valid"<<std::endl;
         return false;
       }
       
       if (Halfedges[Halfedges[i].Next].Origin==Halfedges[i].Origin){  //a degenerate edge{
         if (verbose) std::cout<<"Halfedge "<<i<<" with twin"<<Halfedges[i].Twin<<" is degenerate with vertex "<<Halfedges[i].Origin<<std::endl;
         return false;
       }
       
       if (Halfedges[i].Twin>=0){
         if (Halfedges[Halfedges[i].Twin].Twin!=i){
           if (verbose) std::cout<<"Halfedge "<<i<<"twin is "<<Halfedges[i].Twin<<" which doesn't point back"<<std::endl;
           return false;
         }
         
         if (!Halfedges[Halfedges[i].Twin].Valid){
           if (verbose) std::cout<<"halfedge "<<i<<" is twin with invalid halfedge"<<Halfedges[i].Twin<<std::endl;
           return false;
         }
       }
       
       if (!Halfedges[Halfedges[i].Next].Valid){
         if (verbose) std::cout<<"halfedge "<<i<<" has next invalid halfedge"<<Halfedges[i].Next<<std::endl;
         return false;
       }
       
       if (!Halfedges[Halfedges[i].Prev].Valid){
         if (verbose) std::cout<<"halfedge "<<i<<" has prev invalid halfedge"<<Halfedges[i].Prev<<std::endl;
         return false;
       }
       
       if (Halfedges[i].isFunction){  //checking that it is not left alone
         if (Halfedges[i].Prev==Halfedges[i].Twin){
           if (verbose) std::cout<<"Hex halfedge "<<i<<" has Halfedge "<<Halfedges[i].Prev<<" and both prev and twin"<<std::endl;
           return false;
         }
         
         
         if (Halfedges[i].Next==Halfedges[i].Twin){
           if (verbose) std::cout<<"Hex halfedge "<<i<<" has Halfedge "<<Halfedges[i].Next<<" and both next and twin"<<std::endl;
           return false;
         }
       }
     }
     
     std::vector<std::set<int>> halfedgesinFace(Faces.size());
     std::vector<std::set<int>> verticesinFace(Faces.size());
     for (int i=0;i<Faces.size();i++){
       if (!Faces[i].Valid)
         continue;
       
       //if (Faces[i].NumVertices<3)  //we never allow this
       //  return false;
       int hebegin=Faces[i].AdjHalfedge;
       int heiterate=hebegin;
       int NumEdges=0;
       int actualNumVertices=0;
       
       do{
         if (verticesinFace[i].find(Halfedges[heiterate].Origin)!=verticesinFace[i].end())
           if (verbose) std::cout<<"Warning: Vertex "<<Halfedges[heiterate].Origin<<" appears more than once in face "<<i<<std::endl;
         
         verticesinFace[i].insert(Halfedges[heiterate].Origin);
         halfedgesinFace[i].insert(heiterate);
         actualNumVertices++;
         if (!Halfedges[heiterate].Valid)
           return false;
         
         if (Halfedges[heiterate].AdjFace!=i){
           if (verbose) std::cout<<"Face "<<i<<" has halfedge "<<heiterate<<" that does not point back"<<std::endl;
           return false;
         }
       
         heiterate=Halfedges[heiterate].Next;
         NumEdges++;
         if (NumEdges>Halfedges.size()){
           if (verbose) std::cout<<"Infinity loop!"<<std::endl;
           return false;
         }
          
         
       }while (heiterate!=hebegin);
       
       /*if (actualNumVertices!=Faces[i].NumVertices){
         DebugLog<<"Faces "<<i<<" lists "<<Faces[i].NumVertices<<" vertices but has a chain of "<<actualNumVertices<<endl;
         return false;
       }
       
       for (int j=0;j<Faces[i].NumVertices;j++){
         if (Faces[i].Vertices[j]<0){
           DebugLog<<"Faces "<<i<<".vertices "<<j<<"is undefined"<<endl;
           return false;
         }
         if (!Vertices[Faces[i].Vertices[j]].Valid){
           DebugLog<<"Faces "<<i<<".vertices "<<j<<"is not valid"<<endl;
           return false;
         }
       }*/
     }
     
     //checking if all halfedges that relate to a face are part of its recognized chain (so no floaters)
     for (int i=0;i<Halfedges.size();i++){
       if (!Halfedges[i].Valid)
         continue;
       int currFace = Halfedges[i].AdjFace;
       if (halfedgesinFace[currFace].find(i)==halfedgesinFace[currFace].end()){
         if (verbose) std::cout<<"Halfedge "<<i<<" is floating in face "<<currFace<<std::endl;
         return false;
       }
     }
     
     //check if mesh is a manifold: every halfedge appears only once
     if (checkHalfedgeRepetition){
       std::set<TwinFinder> HESet;
       for (int i=0;i<Halfedges.size();i++){
         if (!Halfedges[i].Valid)
           continue;
         std::set<TwinFinder>::iterator HESetIterator=HESet.find(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
         if (HESetIterator!=HESet.end()){
           if (verbose) std::cout<<"Warning: the halfedge ("<<Halfedges[i].Origin<<","<<Halfedges[Halfedges[i].Next].Origin<<") appears at least twice in the mesh"<<std::endl;
           if (verbose) std::cout<<"for instance halfedges "<<i<<" and "<<HESetIterator->index<<std::endl;
           return false;
           //return false;
         }else{
           HESet.insert(TwinFinder(i,Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
           //if (verbose) std::cout<<"inserting halfedge "<<i<<" which is "<<Halfedges[i].Origin<<", "<<Halfedges[Halfedges[i].Next].Origin<<endl;
         }
       }
     }
     
     if (CheckTwinGaps){
       std::set<TwinFinder> HESet;
       //checking if there is a gap: two halfedges that share the same opposite vertices but do not have twins
       for (int i=0;i<Halfedges.size();i++){
         if (!Halfedges[i].Valid)
           continue;
         
         std::set<TwinFinder>::iterator HESetIterator=HESet.find(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
         if (HESetIterator==HESet.end()){
           HESet.insert(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
           continue;
         }
         
         HESetIterator=HESet.find(TwinFinder(i, Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
         if (HESetIterator!=HESet.end()){
           
           if (Halfedges[i].Twin==-1){
             if (verbose) std::cout<<"Halfedge "<<i<<"has no twin although halfedge "<<HESetIterator->index<<" can be a twin"<<std::endl;
             return false;
           }
           if (Halfedges[HESetIterator->index].Twin==-1){
             if (verbose) std::cout<<"Halfedge "<<HESetIterator->index<<"has no twin although halfedge "<<i<<" can be a twin"<<std::endl;
             return false;
           }
         }
       }
     }
     
     //checking if there are pure boundary faces (there shouldn't be)
     if (checkPureBoundary){
      for (int i=0;i<Halfedges.size();i++){
        if (!Halfedges[i].Valid)
          continue;
        if ((Halfedges[i].Twin<0)&&(Halfedges[i].isFunction))
          if (verbose) std::cout<<"WARNING: Halfedge "<<i<<" is a hex edge without twin!"<<std::endl;
        
        if (Halfedges[i].Twin>0)
          continue;
        
        bool pureBoundary=true;
        int hebegin=i;
        int heiterate=hebegin;
        do{
          if (Halfedges[heiterate].Twin>0){
            pureBoundary=false;
            break;
          }
          heiterate = Halfedges[heiterate].Next;
        }while (heiterate!=hebegin);
        if (pureBoundary){
          if (verbose) std::cout<<"Face "<<Halfedges[i].AdjFace<<"is a pure boundary face!"<<std::endl;
          return false;
        }
      }
        
       //checking for latent valence 2 faces
       std::vector<int> Valences(Vertices.size());
       for (int i=0;i<Vertices.size();i++)
         Valences[i]=0;
       
       for (int i=0;i<Halfedges.size();i++){
         if (Halfedges[i].Valid){
           Valences[Halfedges[i].Origin]++;
           //Valences[Halfedges[Halfedges[i].Next].Origin]++;
           if (Halfedges[i].Twin<0)  //should account for the target as well
             Valences[Halfedges[Halfedges[i].Next].Origin]++;
         }
       }
        
       int countThree;
       for (int i=0;i<Faces.size();i++){
         if (!Faces[i].Valid)
           continue;
         countThree=0;
         int hebegin = Faces[i].AdjHalfedge;
         int heiterate=hebegin;
         int numEdges=0;
         do{
           if (Valences[Halfedges[heiterate].Origin]>2)
             countThree++;
           heiterate=Halfedges[heiterate].Next;
           numEdges++;
           if (numEdges>Halfedges.size()){
             if (verbose) std::cout<<"Infinity loop in face "<<i<<"!"<<std::endl;
             return false;
           }
         }while (heiterate!=hebegin);
         if (countThree<3){
           if (verbose) std::cout<<"Face "<<i<<" is a latent valence 2 face!"<<std::endl;
           if (verbose) std::cout<<"Its vertices are "<<std::endl;
           do{
             if (verbose) std::cout<<"Vertex "<<Halfedges[heiterate].Origin<<" halfedge "<<heiterate<<" valence "<<Valences[Halfedges[heiterate].Origin]<<std::endl;
             
             if (Valences[Halfedges[heiterate].Origin]>2)
               countThree++;
             heiterate=Halfedges[heiterate].Next;
             numEdges++;
             if (numEdges>Halfedges.size()){
               if (verbose) std::cout<<"Infinity loop in face "<<i<<"!"<<std::endl;
               return false;
             }
           }while (heiterate!=hebegin);
           
           //return false;
         }
       }
     }
     
     if (verbose) std::cout<<"Mesh is clear according to given checks"<<std::endl;
     return true;  //most likely the mesh is solid
    
  }
  void CleanMesh(){
    //removing nonvalid vertices
    std::vector<int> TransVertices(Vertices.size());
    std::vector<Vertex> NewVertices;
    for (int i=0;i<Vertices.size();i++){
      if (!Vertices[i].Valid)
        continue;
      
      Vertex NewVertex=Vertices[i];
      NewVertex.ID=NewVertices.size();
      NewVertices.push_back(NewVertex);
      TransVertices[i]=NewVertex.ID;
    }
    
    
    
    Vertices=NewVertices;
    for (int i=0;i<Halfedges.size();i++)
      Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];
    

    
    /*for (int i=0;i<Faces.size();i++){
      if (!Faces[i].Valid)
        continue;
      for (int j=0;j<Faces[i].NumVertices;j++){
        DebugLog<<"i is "<<i<<", j is "<<j<<", Faces[i].Vertices[j] is "<<Faces[i].Vertices[j]<<endl;
        Faces[i].Vertices[j]=TransVertices[Faces[i].Vertices[j]];
      }
    }*/
    
   
    
    //removing nonvalid faces
    std::vector<Face> NewFaces;
    std::vector<int> TransFaces(Faces.size());
    for (int i=0;i<Faces.size();i++){
      if (!Faces[i].Valid)
        continue;
      
      Face NewFace=Faces[i];
      NewFace.ID=NewFaces.size();
      NewFaces.push_back(NewFace);
      TransFaces[i]=NewFace.ID;
    }
    Faces=NewFaces;
    for (int i=0;i<Halfedges.size();i++)
      Halfedges[i].AdjFace=TransFaces[Halfedges[i].AdjFace];
    
    
    
    //removing nonvalid halfedges
    std::vector<Halfedge> NewHalfedges;
    std::vector<int> TransHalfedges(Halfedges.size());
    for (int i=0;i<Halfedges.size();i++){
      if (!Halfedges[i].Valid)
        continue;
      
      Halfedge NewHalfedge=Halfedges[i];
      NewHalfedge.ID=NewHalfedges.size();
      NewHalfedges.push_back(NewHalfedge);
      TransHalfedges[i]=NewHalfedge.ID;
    }
    
     
    
    Halfedges=NewHalfedges;
    for (int i=0;i<Faces.size();i++)
      Faces[i].AdjHalfedge=TransHalfedges[Faces[i].AdjHalfedge];
    
  
    
    for (int i=0;i<Vertices.size();i++)
      Vertices[i].AdjHalfedge=TransHalfedges[Vertices[i].AdjHalfedge];
    
    
    for (int i=0;i<Halfedges.size();i++){
      if (Halfedges[i].Twin!=-1)
        Halfedges[i].Twin=TransHalfedges[Halfedges[i].Twin];
      Halfedges[i].Next=TransHalfedges[Halfedges[i].Next];
      Halfedges[i].Prev=TransHalfedges[Halfedges[i].Prev];
    }
    
  }
  void ComputeTwins(){
    //twinning up edges
    std::set<TwinFinder> Twinning;
    for (int i=0;i<Halfedges.size();i++){
      if (Halfedges[i].Twin>=0)
        continue;
      
      std::set<TwinFinder>::iterator Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
      if (Twinit!=Twinning.end()){
        Halfedges[Twinit->index].Twin=i;
        Halfedges[i].Twin=Twinit->index;
        Twinning.erase(*Twinit);
      } else {
        Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
      }
    }
  }
  void WalkBoundary(int &CurrEdge){
    do{
       CurrEdge=Halfedges[CurrEdge].Next;
       if (Halfedges[CurrEdge].Twin<0)
         break;  //next boundary over a 2-valence vertex
       CurrEdge=Halfedges[CurrEdge].Twin;
     }while (Halfedges[CurrEdge].Twin>=0);
  }
  void RemoveVertex(int vindex, std::deque<int>& removeVertexQueue){
     int hebegin = Vertices[vindex].AdjHalfedge;
     int heiterate=hebegin;
     do{
       if (heiterate==-1){  //boundary vertex
         return;
       }
       heiterate=Halfedges[Halfedges[heiterate].Prev].Twin;
     }while(heiterate!=hebegin);

     Vertices[vindex].Valid = false;
     
     int remainingFace = Halfedges[hebegin].AdjFace;
    
     
     Faces[remainingFace].AdjHalfedge= Halfedges[hebegin].Next;
     heiterate=hebegin;
     int infinityCounter=0;
     do{
       
       int NextEdge = Halfedges[heiterate].Next;
       int PrevEdge = Halfedges[Halfedges[heiterate].Twin].Prev;
      
       Halfedges[NextEdge].Prev = PrevEdge;
       Halfedges[PrevEdge].Next = NextEdge;
       if (Halfedges[NextEdge].AdjFace!=remainingFace)
         Faces[Halfedges[NextEdge].AdjFace].Valid = false;
       
       if (Halfedges[PrevEdge].AdjFace!=remainingFace)
         Faces[Halfedges[PrevEdge].AdjFace].Valid = false;
       
       
       Halfedges[PrevEdge].AdjFace = Halfedges[NextEdge].AdjFace = remainingFace;
       Halfedges[heiterate].Valid = false;
       Halfedges[Halfedges[heiterate].Twin].Valid = false;
       heiterate=Halfedges[Halfedges[heiterate].Prev].Twin;
       infinityCounter++;
       if (infinityCounter>Halfedges.size())
         return;
       
     }while(heiterate!=hebegin);
     
     //cleaning new face
     hebegin =Faces[remainingFace].AdjHalfedge;
     //Faces[remainingFace].NumVertices=0;
     heiterate=hebegin;
     infinityCounter=0;
     do{
       //Faces[remainingFace].NumVertices++;
       Halfedges[heiterate].AdjFace=remainingFace;
       Vertices[Halfedges[heiterate].Origin].AdjHalfedge=heiterate;
       removeVertexQueue.push_front(Halfedges[heiterate].Origin);
       infinityCounter++;
       if (infinityCounter>Halfedges.size())
         return;
       
       heiterate=Halfedges[heiterate].Next;
     }while (heiterate!=hebegin);
  }
  

  void TestUnmatchedTwins();
  
  
  void GenerateMesh(const Eigen::VectorXd& funcOrientations, NFunctionMesher& funcMesh){
    
    using namespace std;
    using namespace Eigen;
    using namespace ::CGAL;
    
    
    funcMesh.Vertices.clear();
    funcMesh.Halfedges.clear();
    funcMesh.Faces.clear();
    
    int numParamFuncs=Halfedges[0].exactParamFuncs.size();
    
    //DebugLog.open("Debugging.txt");
    
    //resolution is set to 10e-6 of bounding box of mesh
    vector<Point3D> coordList;
    for (int i=0;i<Vertices.size();i++)
      coordList.push_back(Vertices[i].Coordinates);
    
    Bbox_3 boundBox = ::CGAL::bbox_3  ( coordList.begin(), coordList.end());
    
    /*double minRange = 3276700.0;
     for (int i=0;i<2;i++)
     minRange=std::min(minRange, boundBox.max(i)-boundBox.min(i));*/
    
    long Resolution=1e8; //pow(10,ceil(10/log10(minRange)));
    //cout<<"Resolution: "<<Resolution<<endl;
    
    for (int findex=0;findex<Faces.size();findex++){
      
      //building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face
      
      int ebegin=Faces[findex].AdjHalfedge;
      int eiterate=ebegin;
      //vector<Point2D> TriPoints2D;
      
      //basis for triangle
      /*do{
       Point2D Location((Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis1,(Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis2);
       TriPoints2D.push_back(Location);
       eiterate=Halfedges[eiterate].Next;
       }while (eiterate!=ebegin);
       Triangle2D CurrTri(TriPoints2D[0], TriPoints2D[1], TriPoints2D[2]);*/
      
      vector<vector<ENumber> > funcValues(3);
      
      //DebugLog<<"Working on triangle "<<findex<<"\n";
      vector<ENumber> minFuncs(numParamFuncs);
      vector<ENumber> maxFuncs(numParamFuncs);
      for (int k=0;k<numParamFuncs;k++){
        minFuncs[k]=ENumber(327600);
        maxFuncs[k]=ENumber(-327600.0);
      }
      
      Arr_2 ParamArr,TriangleArr, FullArr;
      ebegin=Faces[findex].AdjHalfedge;
      eiterate=ebegin;
      int currVertex=0;
      do{
        for(int i=0;i<numParamFuncs;i++){
          if (Halfedges[eiterate].exactParamFuncs[i]>maxFuncs[i]) maxFuncs[i]=Halfedges[eiterate].exactParamFuncs[i];
          if (Halfedges[eiterate].exactParamFuncs[i]<minFuncs[i]) minFuncs[i]=Halfedges[eiterate].exactParamFuncs[i];
        }
        funcValues[currVertex++]=Halfedges[eiterate].exactParamFuncs;
        eiterate=Halfedges[eiterate].Next;
      }while (eiterate!=ebegin);
      
      //building the one-triangle arrangement
      ebegin=Faces[findex].AdjHalfedge;
      eiterate=ebegin;
      vector<EPoint2D> ETriPoints2D;
      //vector<Point2D> TriPoints;
      vector<EPoint3D> ETriPoints3D;
      vector<EdgeData> EdgeDatas;
      ETriPoints2D.push_back(EPoint2D(0,0));
      ETriPoints2D.push_back(EPoint2D(1,0));
      ETriPoints2D.push_back(EPoint2D(0,1));
      do{
        //cout<<"Halfedges[eiterate].Origin: "<<Halfedges[eiterate].Origin<<endl;
        //Point2D Location((Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis1,(Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis2);
        //cout<<"Location: "<<Location<<endl;
        
        Point3D Position=Vertices[Halfedges[eiterate].Origin].Coordinates;
        //ENumber cx=ENumber((int)(Location.x()*Resolution),Resolution);
        //ENumber cy=ENumber((int)(Location.y()*Resolution),Resolution);
        ENumber x=ENumber((long)round(Position.x()*Resolution),Resolution);
        ENumber y=ENumber((long)round(Position.y()*Resolution),Resolution);
        ENumber z=ENumber((long)round(Position.z()*Resolution),Resolution);
        //ETriPoints.push_back(EPoint2D(cx,cy));
        //TriPoints.push_back(Location);
        ETriPoints3D.push_back(EPoint3D(x,y,z));
        int DomEdge;
        
        if ((Halfedges[eiterate].Twin<0)||(Halfedges[eiterate].Twin>eiterate))
          DomEdge=eiterate;
        else
          DomEdge=Halfedges[eiterate].Twin;
        EdgeData ed; ed.OrigHalfedge=DomEdge;
        ed.isBoundary=(Halfedges[eiterate].Twin<0);
        EdgeDatas.push_back(ed);
        eiterate=Halfedges[eiterate].Next;
      }while(ebegin!=eiterate);
      
      for (int i=0;i<3;i++){
        X_monotone_curve_2 c =ESegment2D(ETriPoints2D[i],ETriPoints2D[(i+1)%3]);
        Halfedge_handle he=CGAL::insert_non_intersecting_curve(TriangleArr,c);
        he->set_data(EdgeDatas[i]);
        if (EdgeDatas[i].isBoundary)
          he->source()->data()=he->target()->data()=0;
        else
          he->source()->data()=he->target()->data()=1;
        
        he->twin()->set_data(EdgeDatas[i]);
      }
      
      for (Face_iterator fi= TriangleArr.faces_begin(); fi != TriangleArr.faces_end(); fi++){
        if (fi->is_unbounded())
          fi->data()=0;
        else
          fi->data()=1;
      }
      
      //creating the primal arrangement of lines
      vector<ELine2D> paramLines;
      vector<EDirection2D> isoDirections(numParamFuncs);
      //int jumps = (numParamFuncs%2==0 ? 2 : 1);
      for (int funcIter=0;funcIter<numParamFuncs/*/jumps*/;funcIter++){
        
        vector<EInt> isoValues;
        //cout<<"isoValues: "<<endl;
        EInt q,r;
        CGAL::div_mod(minFuncs[funcIter].numerator(), minFuncs[funcIter].denominator(), q, r);
        EInt minIsoValue = q + (r<0 ? -1 : 0);
        CGAL::div_mod(maxFuncs[funcIter].numerator(), maxFuncs[funcIter].denominator(), q, r);
        EInt maxIsoValue = q + (r<0  ? 0 : -1);
        for (EInt isoValue=minIsoValue-2;isoValue <=maxIsoValue+2;isoValue++){
          //cout<<"isoValue: "<<isoValue<<endl;
          isoValues.push_back(isoValue);
        }
        
        //computing gradient of function in plane
        EVector2D e01 =ETriPoints2D[1] - ETriPoints2D[0];
        EVector2D e12 =ETriPoints2D[2] - ETriPoints2D[1];
        EVector2D e20 =ETriPoints2D[0] - ETriPoints2D[2];
        
        //a and b values of lines
        EVector2D gradVector = funcValues[2][funcIter]*EVector2D(-e01.y(), e01.x())+
        funcValues[0][funcIter]*EVector2D(-e12.y(), e12.x())+
        funcValues[1][funcIter]*EVector2D(-e20.y(), e20.x());
        
        isoDirections[funcIter]=EDirection2D(gradVector);
        
        //Number avgFuncValue = (funcValues[0](funcIter)+funcValues[1](funcIter)+funcValues[2](funcIter))/3.0;
        //TODO: find c = z1*u+z2 of ax+by+c(u) ad then use it to generate all values between floor and ceil.
        
        //pinv of [a 1;b 1;c 1] is [           2*a - b - c,           2*b - a - c,           2*c - b - a]
        //[ b^2 - a*b + c^2 - a*c, a^2 - b*a + c^2 - b*c, a^2 - c*a + b^2 - c*b]/(2*a^2 - 2*a*b - 2*a*c + 2*b^2 - 2*b*c + 2*c^2)
        
        ENumber a=funcValues[0][funcIter];
        ENumber b=funcValues[1][funcIter];
        ENumber c=funcValues[2][funcIter];
        if ((a==b)&&(b==c))
          continue;  //that means a degenerate function on the triangle
        
        //cout<<"a,b,c: "<<a.to_double()<<","<<b.to_double()<<","<<c.to_double()<<endl;
        
        ENumber rhs[3];
        rhs[0]=-gradVector[0]*ETriPoints2D[0].x()-gradVector[1]*ETriPoints2D[0].y();
        rhs[1]=-gradVector[0]*ETriPoints2D[1].x()-gradVector[1]*ETriPoints2D[1].y();
        rhs[2]=-gradVector[0]*ETriPoints2D[2].x()-gradVector[1]*ETriPoints2D[2].y();
        
        ENumber invM[2][3];
        invM[0][0]= 2*a-b-c;
        invM[0][1]= 2*b-a-c;
        invM[0][2]= 2*c-b-a;
        invM[1][0]=b*b - a*b + c*c - a*c;
        invM[1][1]=a*a - b*a + c*c - b*c;
        invM[1][2]=a*a - c*a + b*b - c*b;
        for (int row=0;row<2;row++)
          for (int col=0;col<3;col++)
            invM[row][col]/=(ENumber(2)*(a*a - a*b - a*c + b*b- b*c + c*c));
        
        //cout<<(ENumber(2)*(a*a - a*b - a*c + b*b- b*c + c*c)).to_double()<<endl;
        
        ENumber x[2];
        x[0] = invM[0][0]*rhs[0]+invM[0][1]*rhs[1]+invM[0][2]*rhs[2];
        x[1] = invM[1][0]*rhs[0]+invM[1][1]*rhs[1]+invM[1][2]*rhs[2];
        
        
        //RowVectorXd x = lhs.colPivHouseholderQr().solve(rhs).transpose();
        
        //sanity check
        ENumber error[3];
        error[0]=x[0]*a+x[1]-rhs[0];
        error[1]=x[0]*b+x[1]-rhs[1];
        error[2]=x[0]*c+x[1]-rhs[2];
        
        
        //cout<<"lhs*x - rhs: "<<error[0].to_double()<<","<<error[1].to_double()<<","<<error[2].to_double()<<endl;
        
        //full sanity check
        /*MatrixXd lhss(3,2);
         lhss<<funcValues[0][funcIter].to_double(), 1.0,
         funcValues[1][funcIter].to_double(), 1.0,
         funcValues[2][funcIter].to_double(), 1.0;
         
         VectorXd rhss(3);
         rhss<<-gradVector[0].to_double()*ETriPoints2D[0].x().to_double()-gradVector[1].to_double()*ETriPoints2D[0].y().to_double(),
         -gradVector[0].to_double()*ETriPoints2D[1].x().to_double()-gradVector[1].to_double()*ETriPoints2D[1].y().to_double(),
         -gradVector[0].to_double()*ETriPoints2D[2].x().to_double()-gradVector[1].to_double()*ETriPoints2D[2].y().to_double();
         
         /*cout<<"invM: "<<endl;
         for (int r=0;r<2;r++)
         for (int c=0;c<3;c++)
         cout<<invM[r][c].to_double()<<","<<endl;
         
         
         
         
         MatrixXd invlhs = (lhss.transpose()*lhss).inverse()*lhss.transpose();
         cout<<"invLhs: "<<invlhs<<endl;
         RowVectorXd xx = lhss.colPivHouseholderQr().solve(rhss).transpose();
         xx =invlhs*rhss;
         cout<<"lhss*xx - rhss"<<lhss*xx.transpose() - rhss<<endl;*/
        
        
        //generating all lines
        for (int isoIndex=0;isoIndex<isoValues.size();isoIndex++){
          //ENumber isoVec[2];
          //isoVec[0]=isoValues[isoIndex];
          //isoVec[1]= ENumber(1);
          ENumber currc = isoValues[isoIndex]*x[0]+x[1];
          // ENumber a=ENumber((int)(gradVector[0]*Resolution),Resolution);
          //ENumber b=ENumber((int)(gradVector[1]*Resolution),Resolution);
          //ENumber c=ENumber((int)(currc(0)*Resolution),Resolution);
          paramLines.push_back(ELine2D(gradVector[0],gradVector[1],currc));
          //cout<<"paramLine: "<<gradVector[0]<<","<<gradVector[1]<<","<<currc<<endl;
        }
      }
      
      //cout<<"paramLines.size() :"<<paramLines.size()<<endl;
      CGAL::insert(ParamArr, paramLines.begin(), paramLines.end());
      
      //giving edge data to curve arrangement
      Arr_2::Edge_iterator                  eit;
      Arr_2::Originating_curve_iterator     ocit;
      for (eit = ParamArr.edges_begin(); eit != ParamArr.edges_end(); ++eit) {
        for (ocit = ParamArr.originating_curves_begin(eit);
             ocit != ParamArr.originating_curves_end(eit); ++ocit){
          EDirection2D thisDirection =  EDirection2D(ocit->supporting_line().a(), ocit->supporting_line().b());
          //cout<<"thisDirection: "<<thisDirection<<endl;
          for (int paramIter = 0;paramIter<numParamFuncs/*/jumps*/;paramIter++){
            //cout<<"isoDirections[paramIter]: "<<isoDirections[paramIter]<<endl;
            if ((thisDirection==isoDirections[paramIter])||(thisDirection==-isoDirections[paramIter])){
              eit->data().funcNum=paramIter;
              eit->twin()->data().funcNum=paramIter;
              //cout<<"assigning "<<paramIter<<endl;
            }
          }
        }
      }
      
      
      //sanity check: all edges are assigned
      /*cout<<"paramarr edges: "<<endl;
       for (eit = ParamArr.edges_begin(); eit != ParamArr.edges_end(); ++eit)
       cout<<"paramarr eit->data().funcNum: "<<eit->data().funcNum<<endl;*/
      
      //
      //creating the overlay
      Overlay_traits ot;
      overlay (TriangleArr, ParamArr, FullArr, ot);
      
      /*cout<<"FullArr edges: "<<endl;
       for (eit = FullArr.edges_begin(); eit != FullArr.edges_end(); ++eit)
       cout<<"FullArr eit->data().funcNum: "<<eit->data().funcNum<<endl;*/
      
      
      for (Face_iterator fi=FullArr.faces_begin();fi!=FullArr.faces_end();fi++){
        if (!fi->data())
          continue;  //not participating
        
        Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
        Ccb_halfedge_circulator heiterate=hebegin;
        do{
          
          if (heiterate->source()->data()<0){  //new vertex
            Vertex NewVertex;
            NewVertex.ID=funcMesh.Vertices.size();
            NewVertex.isFunction=(heiterate->source()->data()==-2);
            funcMesh.Vertices.push_back(NewVertex);
            heiterate->source()->data()=NewVertex.ID;
          }
          
          if (heiterate->data().ID<0){  //new halfedge
            Halfedge NewHalfedge;
            NewHalfedge.ID=funcMesh.Halfedges.size();
            NewHalfedge.isFunction=(heiterate->data().ID==-2);
            NewHalfedge.Origin=heiterate->source()->data();
            NewHalfedge.OrigHalfedge=heiterate->data().OrigHalfedge;
            NewHalfedge.OrigParamFunc=heiterate->data().funcNum;
            //cout<<"NewHalfedge.OrigParamFunc :"<<NewHalfedge.OrigParamFunc<<endl;
            funcMesh.Vertices[heiterate->source()->data()].AdjHalfedge=NewHalfedge.ID;
            funcMesh.Halfedges.push_back(NewHalfedge);
            heiterate->data().ID=NewHalfedge.ID;
          }
          heiterate++;
        }while(heiterate!=hebegin);
        
        //now assigning nexts and prevs
        do{
          funcMesh.Halfedges[heiterate->data().ID].Next=heiterate->next()->data().ID;
          funcMesh.Halfedges[heiterate->data().ID].Prev=heiterate->prev()->data().ID;
          funcMesh.Halfedges[heiterate->data().ID].Twin=heiterate->twin()->data().ID;
          if (heiterate->twin()->data().ID>=0)
            funcMesh.Halfedges[heiterate->twin()->data().ID].Twin=heiterate->data().ID;
          
          heiterate++;
        }while (heiterate!=hebegin);
      }
      
      //constructing the actual vertices
      for (Vertex_iterator vi=FullArr.vertices_begin();vi!=FullArr.vertices_end();vi++){
        if (vi->data()<0)
          continue;
        
        //finding out barycentric coordinates
        ENumber BaryValues[3];
        ENumber Sum=0;
        for (int i=0;i<3;i++){
          ETriangle2D t(vi->point(), ETriPoints2D[(i+1)%3], ETriPoints2D[(i+2)%3]);
          BaryValues[i]=t.area();
          Sum+=BaryValues[i];
        }
        for (int i=0;i<3;i++)
          BaryValues[i]/=Sum;
        
        EPoint3D ENewPosition(0,0,0);
        for (int i=0;i<3;i++)
          ENewPosition=ENewPosition+(ETriPoints3D[i]-CGAL::ORIGIN)*BaryValues[i];
        
        Point3D NewPosition(to_double(ENewPosition.x()), to_double(ENewPosition.y()), to_double(ENewPosition.z()));
        funcMesh.Vertices[vi->data()].Coordinates=NewPosition;
        funcMesh.Vertices[vi->data()].ECoordinates=ENewPosition;
        
        //DebugLog<<"Creating Vertex "<<vi->data()<<" with 2D coordinates ("<<vi->point().x()<<","<<vi->point().y()<<") "<<" and 3D Coordinates ("<<std::setprecision(10) <<NewPosition.x()<<","<<NewPosition.y()<<","<<NewPosition.z()<<")\n";
      }
      
      for (Face_iterator fi=FullArr.faces_begin();fi!=FullArr.faces_end();fi++){
        if (!fi->data())
          continue;
        
        int FaceSize=0;
        Ccb_halfedge_circulator hebegin=fi->outer_ccb ();
        Ccb_halfedge_circulator heiterate=hebegin;
        do{ FaceSize++;  heiterate++; }while(heiterate!=hebegin);
        int CurrPlace=0;
        
        Face NewFace;
        NewFace.ID=funcMesh.Faces.size();
        //NewFace.NumVertices=FaceSize;
        NewFace.AdjHalfedge=hebegin->data().ID;
        
        do{
          //NewFace.Vertices[CurrPlace++]=heiterate->source()->data();
          funcMesh.Halfedges[heiterate->data().ID].AdjFace=NewFace.ID;
          heiterate++;
        }while(heiterate!=hebegin);
        funcMesh.Faces.push_back(NewFace);
      }
      
    }
    
    //devising angles from differences in functions
    //int ratio = (numParamFuncs%2==0 ? 1 : 2);
    for (int hi=0;hi<funcMesh.Halfedges.size();hi++){
      //cout<<"funcMesh.Halfedges[hi].OrigParamFunc: "<<funcMesh.Halfedges[hi].OrigParamFunc<<endl;
      //cout<<"funcMesh.Halfedges[Halfedges[hi].Prev].OrigParamFunc: "<<funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigParamFunc<<endl;
      if ((funcMesh.Halfedges[hi].OrigParamFunc==-1)||(funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigParamFunc==-1))
        funcMesh.Halfedges[hi].prescribedAngle=-1.0;  //one of the edges is a triangle edge, and it will be devised later.
      else{
        //int func1 =(ratio*(funcMesh.Halfedges[hi].OrigParamFunc)) % (numParamFuncs/(3-ratio));
        //int func2 =(ratio*(funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigParamFunc)) % (numParamFuncs/(3-ratio));
        
        double funcOrient1 = funcOrientations(funcMesh.Halfedges[hi].OrigParamFunc);
        double funcOrient2 = funcOrientations(funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigParamFunc);
        //cout<<"funcOrient1: "<<funcOrient1<<endl;
        //cout<<"funcOrient2: "<<funcOrient2<<endl;
        funcMesh.Halfedges[hi].prescribedAngle=funcOrient2-funcOrient1;
        //getting difference between [-pi,pi]
        while (funcMesh.Halfedges[hi].prescribedAngle>igl::PI)
          funcMesh.Halfedges[hi].prescribedAngle-=2*igl::PI;
        while (funcMesh.Halfedges[hi].prescribedAngle<-igl::PI)
          funcMesh.Halfedges[hi].prescribedAngle+=2*igl::PI;
        //cout<<"After 2pi correction: "<<funcMesh.Halfedges[hi].prescribedAngle<<endl;
        if (funcMesh.Halfedges[hi].prescribedAngle<0)
          funcMesh.Halfedges[hi].prescribedAngle+=igl::PI;//+funcMesh.Halfedges[hi].prescribedAngle;
      
      }
    }
  }
  
  struct PointPair{
    int Index1, Index2;
    ENumber Distance;
    
    PointPair(int i1, int i2, ENumber d):Index1(i1), Index2(i2), Distance(d){}
    ~PointPair(){}
    
    const bool operator<(const PointPair& pp) const {
      if (Distance>pp.Distance) return false;
      if (Distance<pp.Distance) return true;
      
      if (Index1>pp.Index1) return false;
      if (Index1<pp.Index1) return true;
      
      if (Index2>pp.Index2) return false;
      if (Index2<pp.Index2) return true;
      
      return false;
      
    }
  };

  std::vector<std::pair<int,int>> FindVertexMatch(const bool verbose, std::vector<EPoint3D>& Set1, std::vector<EPoint3D>& Set2)
  {
    std::set<PointPair> PairSet;
    for (int i=0;i<Set1.size();i++)
      for (int j=0;j<Set2.size();j++)
        PairSet.insert(PointPair(i,j,squared_distance(Set1[i],Set2[j])));
    
    if (Set1.size()!=Set2.size())  //should not happen anymore
      std::cout<<"FindVertexMatch(): The two sets are of different sizes!! "<<std::endl;
    
    //adding greedily legal connections until graph is full
    std::vector<bool> Set1Connect(Set1.size());
    std::vector<bool> Set2Connect(Set2.size());
    
    std::vector<std::pair<int, int> > Result;
    
    for (int i=0;i<Set1.size();i++)
      Set1Connect[i]=false;
    
    for (int i=0;i<Set2.size();i++)
      Set2Connect[i]=false;
    
    /*if (Set1.size()!=Set2.size())
     int kaka=9;*/
    
    int NumConnected=0;
    
    //categorically match both ends
    
    Result.push_back(std::pair<int, int>(0,0));
    Result.push_back(std::pair<int, int>(Set1.size()-1,Set2.size()-1));
    for (std::set<PointPair>::iterator ppi=PairSet.begin();ppi!=PairSet.end();ppi++)
    {
      PointPair CurrPair=*ppi;
      //checking legality - if any of one's former are connected to ones latters or vice versa
      bool FoundConflict=false;
      for (int i=0;i<Result.size();i++){
        if (((Result[i].first>CurrPair.Index1)&&(Result[i].second<CurrPair.Index2))||
            ((Result[i].first<CurrPair.Index1)&&(Result[i].second>CurrPair.Index2))){
          FoundConflict=true;
          break;
        }
      }
      
      if (FoundConflict)
        continue;
      
      //if both are already matched, this matching is redundant
      if ((Set1Connect[CurrPair.Index1])&&(Set2Connect[CurrPair.Index2]))
        continue;  //there is no reason for this matching
      
      //otherwise this edge is legal, so add it
      Result.push_back(std::pair<int, int>(CurrPair.Index1, CurrPair.Index2));
      if (!Set1Connect[CurrPair.Index1]) NumConnected++;
      if (!Set2Connect[CurrPair.Index2]) NumConnected++;
      Set1Connect[CurrPair.Index1]=Set2Connect[CurrPair.Index2]=true;
      /*if (NumConnected==Set1.size()+Set2.size())
       break;  //all nodes are connected*/
    }
    
    for (int i=0;i<Set1.size();i++)
      if ((!Set1Connect[i])&&(verbose))
        std::cout<<"Relative Vertex "<<i<<" in Set1 is unmatched!"<<std::endl;
    
    for (int i=0;i<Set2.size();i++)
      if ((!Set2Connect[i])&&(verbose))
        std::cout<<"Relative Vertex "<<i<<" in Set2 is unmatched!"<<std::endl;
    
    /*if (NumConnected!=Set1.size()+Set2.size())
     int kaka=9;*/
    
    if (verbose){
      for (int i=0;i<Result.size();i++){
        if (squared_distance(Set1[Result[i].first],Set2[Result[i].second])>0){
          std::cout<<"("<<Result[i].first<<","<<Result[i].second<<") with dist "<<squared_distance(Set1[Result[i].first],Set2[Result[i].second])<<std::endl;
          std::cout<<"Distance is abnormally not zero!"<<std::endl;
        }
      }
    }
    
    
    return Result;
    
  }


  struct TwinFinder{
    int index;
    int v1,v2;
    
    TwinFinder(int i, int vv1, int vv2):index(i), v1(vv1), v2(vv2){}
    ~TwinFinder(){}
    
    const bool operator<(const TwinFinder& tf) const
    {
      if (v1<tf.v1) return false;
      if (v1>tf.v1) return true;
      
      if (v2<tf.v2) return false;
      if (v2>tf.v2) return true;
      
      return false;
    }
    
    
  };

  
  bool SimplifyFuncMesh(const bool verbose, int N){
     //unifying vertices which are similar
     
    using namespace std;
    using namespace Eigen;
    using namespace boost;
    
     if (!CheckMesh(verbose, false, false, false))
       return false;
     
     int MaxOrigHE=-3276700.0;
     for (int i=0;i<Halfedges.size();i++)
       MaxOrigHE=std::max(MaxOrigHE, Halfedges[i].OrigHalfedge);
     
     vector<bool> visitedOrig(MaxOrigHE+1);
     for (int i=0;i<MaxOrigHE+1;i++) visitedOrig[i]=false;
     for (int i=0;i<Halfedges.size();i++){
       if (Halfedges[i].OrigHalfedge<0)
         continue;
       if (visitedOrig[Halfedges[i].OrigHalfedge])
         continue;
       
       int hebegin = i;
       int heiterate = hebegin;
       do{
         visitedOrig[Halfedges[heiterate].OrigHalfedge]=true;
         WalkBoundary(heiterate);
       }while (heiterate!=hebegin);
       
     }
     
     vector< vector<int> > BoundEdgeCollect1(MaxOrigHE+1);
     vector< vector<int> > BoundEdgeCollect2(MaxOrigHE+1);
     vector<bool> Marked(Halfedges.size());
     for (int i=0;i<Halfedges.size();i++) Marked[i]=false;
     //finding out vertex correspondence along twin edges of the original mesh by walking on boundaries
     for (int i=0;i<Halfedges.size();i++){
       if ((Halfedges[i].OrigHalfedge<0)||(Marked[i]))
         continue;
       
       //find the next beginning of a boundary
       int PrevOrig;
       int CurrEdge=i;
       do{
         PrevOrig=Halfedges[CurrEdge].OrigHalfedge;
         WalkBoundary(CurrEdge);
       }while(PrevOrig==Halfedges[CurrEdge].OrigHalfedge);
       
       //filling out strips of boundary with the respective attached original halfedges
       int BeginEdge=CurrEdge;
       vector<pair<int,int> > CurrEdgeCollect;
       do{
         CurrEdgeCollect.push_back(pair<int, int> (Halfedges[CurrEdge].OrigHalfedge, CurrEdge));
         Marked[CurrEdge]=true;
         WalkBoundary(CurrEdge);
       }while (CurrEdge!=BeginEdge);
       
       PrevOrig=-1000;
       bool In1;
       for (int j=0;j<CurrEdgeCollect.size();j++){
         if (CurrEdgeCollect[j].first!=PrevOrig)
           In1=BoundEdgeCollect1[CurrEdgeCollect[j].first].empty();
         
         if (In1)
           BoundEdgeCollect1[CurrEdgeCollect[j].first].push_back(CurrEdgeCollect[j].second);
         else
           BoundEdgeCollect2[CurrEdgeCollect[j].first].push_back(CurrEdgeCollect[j].second);
         PrevOrig=CurrEdgeCollect[j].first;
       }
     }
     
     //editing the edges into two vector lists per associated original edge
     vector< vector<int> > VertexSets1(MaxOrigHE+1), VertexSets2(MaxOrigHE+1);
     for (int i=0;i<MaxOrigHE+1;i++){
       for (int j=0;j<BoundEdgeCollect1[i].size();j++)
         VertexSets1[i].push_back(Halfedges[BoundEdgeCollect1[i][j]].Origin);
       
       if (BoundEdgeCollect1[i].size()>0)
         VertexSets1[i].push_back(Halfedges[Halfedges[BoundEdgeCollect1[i][BoundEdgeCollect1[i].size()-1]].Next].Origin);
       
       for (int j=0;j<BoundEdgeCollect2[i].size();j++)
         VertexSets2[i].push_back(Halfedges[BoundEdgeCollect2[i][j]].Origin);
       
       if (BoundEdgeCollect2[i].size()>0)
         VertexSets2[i].push_back(Halfedges[Halfedges[BoundEdgeCollect2[i][BoundEdgeCollect2[i].size()-1]].Next].Origin);
       
       std::reverse(VertexSets2[i].begin(),VertexSets2[i].end());
     }
     
     //finding out vertex matches
     vector<pair<int, int> > VertexMatches;
     for (int i=0;i<MaxOrigHE+1;i++){
       vector<EPoint3D> PointSet1(VertexSets1[i].size());
       vector<EPoint3D> PointSet2(VertexSets2[i].size());
       for (int j=0;j<PointSet1.size();j++)
         PointSet1[j]=Vertices[VertexSets1[i][j]].ECoordinates;
       
       for (int j=0;j<PointSet2.size();j++)
         PointSet2[j]=Vertices[VertexSets2[i][j]].ECoordinates;
       
       vector<pair<int, int> > CurrMatches;
       if ((!PointSet1.empty())&&(!PointSet2.empty()))
         CurrMatches=FindVertexMatch(verbose, PointSet1, PointSet2);
       
       for (int j=0;j<CurrMatches.size();j++){
         CurrMatches[j].first =VertexSets1[i][CurrMatches[j].first];
         CurrMatches[j].second=VertexSets2[i][CurrMatches[j].second];
       }
       
       VertexMatches.insert(VertexMatches.end(), CurrMatches.begin(), CurrMatches.end() );
     }
     
     //finding connected components, and uniting every component into a random single vertex in it (it comes out the last mentioned)
     Graph MatchGraph;
     for (int i=0;i<Vertices.size();i++)
       add_vertex(MatchGraph);
     for (int i=0;i<VertexMatches.size();i++)
       add_edge(VertexMatches[i].first, VertexMatches[i].second, MatchGraph);
     
     double MaxDist=-327670000.0;
     for (int i=0;i<VertexMatches.size();i++)
       MaxDist=std::max(MaxDist, (Vertices[VertexMatches[i].first].Coordinates-Vertices[VertexMatches[i].second].Coordinates).squared_length());
     
    if (verbose)
      std::cout<<"Max matching distance: "<<MaxDist<<endl;
     
     //vector<int> TransVertices(Vertices.size());
    TransVertices.resize(Vertices.size());
    int NumNewVertices = connected_components(MatchGraph, &TransVertices[0]);
    
    if (!CheckMesh(verbose, false, false, false))
       return false;
     
     vector<bool> transClaimed(NumNewVertices);
     for (int i=0;i<NumNewVertices;i++)
       transClaimed[i]=false;
     //unifying all vertices into the TransVertices
     vector<Vertex> NewVertices(NumNewVertices);
     for (int i=0;i<Vertices.size();i++){  //redundant, but not terrible
       if (!Vertices[i].Valid)
         continue;
       Vertex NewVertex=Vertices[i];
       NewVertex.ID=TransVertices[i];
       transClaimed[TransVertices[i]]=true;
       NewVertices[TransVertices[i]]=NewVertex;
     }
     
     for (int i=0;i<NumNewVertices;i++)
       if (!transClaimed[i])
         NewVertices[i].Valid=false;  //this vertex is dead to begin with
       
     Vertices=NewVertices;
     
     for (int i=0;i<Halfedges.size();i++){
       if (!Halfedges[i].Valid)
         continue;
       Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];
       Vertices[Halfedges[i].Origin].AdjHalfedge=i;
     }
     
     
     if (!CheckMesh(verbose, true, false, false))
       return false;
     
     //twinning up edges
     set<TwinFinder> Twinning;
     for (int i=0;i<Halfedges.size();i++){
       if ((Halfedges[i].Twin>=0)||(!Halfedges[i].Valid))
         continue;
       
       set<TwinFinder>::iterator Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
       if (Twinit!=Twinning.end()){
         if ((Halfedges[Twinit->index].Twin!=-1)&&(verbose))
           std::cout<<"warning: halfedge "<<Twinit->index<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<std::endl;
         if ((Halfedges[i].Twin!=-1)&&(verbose))
           std::cout<<"warning: halfedge "<<i<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<std::endl;
         Halfedges[Twinit->index].Twin=i;
         Halfedges[i].Twin=Twinit->index;
         
         if (Halfedges[i].isFunction){
           Halfedges[Twinit->index].isFunction = true;
         } else if (Halfedges[Twinit->index].isFunction){
           Halfedges[i].isFunction = true;
         }
         Twinning.erase(*Twinit);
       } else {
         Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
       }
     }
     
     //check if there are any non-twinned edge which shouldn't be in a closed mesh
    if (verbose){
     for (int i=0;i<Halfedges.size();i++){
       if (Halfedges[i].Twin==-1)
         std::cout<<"Halfedge "<<i<<" does not have a twin!"<<std::endl;
     }
    }
     

     if (!CheckMesh(verbose, true, true, true))
       return false;
     
     //removing triangle components
     
     //starting with pure triangle vertices
     std::vector<bool> isPureTriangle(Vertices.size());
     std::vector<bool> isBoundary(Vertices.size());
     for (int i=0;i<Vertices.size();i++){
       isPureTriangle[i]=true;
       isBoundary[i]=false;
     }
     for (int i=0;i<Halfedges.size();i++){
       if ((Halfedges[i].isFunction)&&(Halfedges[i].Valid)){
         isPureTriangle[Halfedges[i].Origin]=isPureTriangle[Halfedges[Halfedges[i].Next].Origin]=false;  //adjacent to at least one hex edge
       }
       if (Halfedges[i].Twin==-1){
         isBoundary[Halfedges[i].Origin]=true;
         isPureTriangle[Halfedges[i].Origin]=false;  //this shouldn't be removed
       }
     }
     
     std::vector<bool> isEar(Vertices.size());
     for (int i=0;i<Vertices.size();i++){
       isEar[i] = (Halfedges[Vertices[i].AdjHalfedge].Twin==-1)&&(Halfedges[Halfedges[Vertices[i].AdjHalfedge].Prev].Twin==-1);
       if (isEar[i]) isPureTriangle[i]=false;
     }
     
     //realigning halfedges in hex vertices to only follow other hex edges
     for (int i=0;i<Vertices.size();i++){
       if ((isPureTriangle[i])||(!Vertices[i].Valid))
         continue;
       
       vector<int> hexHEorder;
       int hebegin = Vertices[i].AdjHalfedge;
       if (isBoundary[i]){
         //finding the first hex halfedge
         while (Halfedges[Halfedges[hebegin].Prev].Twin!=-1)
           hebegin =Halfedges[Halfedges[hebegin].Prev].Twin;
       }
       
       int heiterate=hebegin;
       do{
         if ((Halfedges[heiterate].isFunction)||(Halfedges[heiterate].Twin==-1))
           hexHEorder.push_back(heiterate);
         if (Halfedges[heiterate].Twin==-1)
           break;
         heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
       }while(heiterate!=hebegin);
       
       
       for (int j=0;j<hexHEorder.size();j++){
         if ((isBoundary[i])&&(j==hexHEorder.size()-1))
           continue;
         Halfedges[hexHEorder[(j+1)%hexHEorder.size()]].Prev =Halfedges[hexHEorder[j]].Twin;
         Halfedges[Halfedges[hexHEorder[j]].Twin].Next =hexHEorder[(j+1)%hexHEorder.size()];
         Vertices[Halfedges[hexHEorder[j]].Origin].AdjHalfedge=hexHEorder[j];
       }
       
       if (isBoundary[i]){ //connect first to the prev
         Halfedges[hexHEorder[0]].Prev = Halfedges[hebegin].Prev;
         Halfedges[Halfedges[hebegin].Prev].Next =hexHEorder[0];
         Vertices[Halfedges[hexHEorder[0]].Origin].AdjHalfedge=hexHEorder[0];
       }
     }

     //invalidating all triangle vertices and edges
     for (int i=0;i<Vertices.size();i++)
       if (isPureTriangle[i])
         Vertices[i].Valid=false;
     
     for (int i=0;i<Halfedges.size();i++)
       if ((!Halfedges[i].isFunction)&&(Halfedges[i].Twin!=-1))
         Halfedges[i].Valid=false;
     
     //realigning faces
     VectorXi visitedHE=VectorXi::Zero(Halfedges.size());
     VectorXi usedFace=VectorXi::Zero(Faces.size());
     for (int i=0;i<Halfedges.size();i++){
       if ((!Halfedges[i].Valid)||(visitedHE[i]!=0))
         continue;
       
       //following the loop and reassigning face
       int currFace=Halfedges[i].AdjFace;
       Faces[currFace].AdjHalfedge=i;
       usedFace[currFace]=1;
       int hebegin=i;
       int heiterate=hebegin;
       int infinityCounter=0;
       do{
         infinityCounter++;
         if (infinityCounter>Halfedges.size()){
           std::cout<<"Infinity loop in realigning faces on halfedge "<<i<<std::endl;
           return false;
         }
         Halfedges[heiterate].AdjFace=currFace;
         heiterate=Halfedges[heiterate].Next;
       }while (heiterate!=hebegin);
     }
     
     int countThree=0;
     for (int i=0;i<Faces.size();i++)
       if (!usedFace[i])
         Faces[i].Valid=false;
     
     
     //killing perfect ear faces (not doing corners atm)
     //counting valences
     vector<int> Valences(Vertices.size());
     for (int i=0;i<Vertices.size();i++)
       Valences[i]=0;
     
     for (int i=0;i<Halfedges.size();i++){
       if (Halfedges[i].Valid){
         Valences[Halfedges[i].Origin]++;
         //Valences[Halfedges[Halfedges[i].Next].Origin]++;
         if (Halfedges[i].Twin<0)  //should account for the target as well
           Valences[Halfedges[Halfedges[i].Next].Origin]++;
       }
     }
        
     for (int i=0;i<Faces.size();i++){
       if (!Faces[i].Valid)
         continue;
       countThree=0;
       int hebegin = Faces[i].AdjHalfedge;
       int heiterate=hebegin;
       do{
         if (Valences[Halfedges[heiterate].Origin]>2)
           countThree++;
         heiterate=Halfedges[heiterate].Next;
       }while (heiterate!=hebegin);
       if (countThree<3){
         do{
           /*DebugLog<<"Invalidating Vertex "<<Halfedges[heiterate].Origin<<"and  halfedge "<<heiterate<<" of valence "<<Valences[Halfedges[heiterate].Origin]<<endl;*/
           
           Halfedges[heiterate].Valid=false;
           if (Halfedges[heiterate].Twin!=-1)
             Halfedges[Halfedges[heiterate].Twin].Twin=-1;
           if ((Halfedges[heiterate].Twin==-1)&&(Halfedges[Halfedges[heiterate].Prev].Twin==-1))  //origin is a boundary vertex
             Vertices[Halfedges[heiterate].Origin].Valid=false;
           
           heiterate=Halfedges[heiterate].Next;
           
           
         }while (heiterate!=hebegin);
         Faces[i].Valid=false;
         //return false;
       }
     }
     
     //need to realign all vertices pointing
     for (int i=0;i<Halfedges.size();i++)
       if (Halfedges[i].Valid)
         Vertices[Halfedges[i].Origin].AdjHalfedge=i;
       
        
     if (!CheckMesh(verbose, true, true, true))
       return false;
     
     for (int i=0;i<Valences.size();i++)
       if ((Vertices[i].Valid)&&(Valences[i]<2))
         Vertices[i].Valid=false;
     
     for (int i=0;i<Vertices.size();i++){
       if ((Vertices[i].Valid)&&(Valences[i]<=2)&&(!isEar[i]))
         UnifyEdges(Vertices[i].AdjHalfedge);
     }
     
     if (!CheckMesh(verbose, true, true, true))
       return false;
     
     //remove non-valid components
     CleanMesh();
     
     //checking if mesh is valid
     if (!CheckMesh(verbose, true, true, true))
       return false;
     
     return true;
     
  }
  
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
                     const Eigen::MatrixXd& nFunction,
                     const int n,
                     const int N,
                     const Eigen::SparseMatrix<int>& n2NMat,
                     const Eigen::VectorXi& integerVars){
    
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
    
    /*vector<ENumber> constraintError;
    exactSparseMult(constraintMatInteger, exactParamFuncsd,constraintError);
    
    ENumber MaxError(0);
    
    for (int i=0;i<constraintError.size();i++)
      if (abs(constraintError[i])>MaxError)
        MaxError=abs(constraintError[i]);
    
    if (verbose)
      cout<<"constraintMatInteger*exactParamFuncsd MaxError: "<<MaxError<<endl;*/
    
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
  
  NFunctionMesher(){}
  ~NFunctionMesher(){}
  
};





} //namespace directional


#endif
