#include "Definitions.h"
#include "Mesh.h"
#include <Eigen/Dense>
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <igl/PI.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <CGAL/Arrangement_with_history_2.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/number_utils.h>
#include <math.h>

using namespace boost;
using namespace CGAL;
using namespace Eigen;
using namespace std;
 


void Mesh::CleanMesh()
{
  DebugLog<<"Cleaning mesh (removing invalid components"<<endl;
  //removing nonvalid vertices
  vector<int> TransVertices(Vertices.size());
  vector<Vertex> NewVertices;
  for (int i=0;i<Vertices.size();i++){
    if (!Vertices[i].Valid)
      continue;
    
    Vertex NewVertex=Vertices[i];
    NewVertex.ID=NewVertices.size();
    NewVertices.push_back(NewVertex);
    TransVertices[i]=NewVertex.ID;
  }
  
   DebugLog<<"ok here 1"<<endl;
  
  Vertices=NewVertices;
  for (int i=0;i<Halfedges.size();i++)
    Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];
  
   DebugLog<<"ok here 2"<<endl;
  
  /*for (int i=0;i<Faces.size();i++){
    if (!Faces[i].Valid)
      continue;
    for (int j=0;j<Faces[i].NumVertices;j++){
      DebugLog<<"i is "<<i<<", j is "<<j<<", Faces[i].Vertices[j] is "<<Faces[i].Vertices[j]<<endl;
      Faces[i].Vertices[j]=TransVertices[Faces[i].Vertices[j]];
    }
  }*/
  
  DebugLog<<"ok here 3"<<endl;
  
  //removing nonvalid faces
  vector<Face> NewFaces;
  vector<int> TransFaces(Faces.size());
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
  
   DebugLog<<"ok here 4"<<endl;
  
  //removing nonvalid halfedges
  vector<Halfedge> NewHalfedges;
  vector<int> TransHalfedges(Halfedges.size());
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    
    Halfedge NewHalfedge=Halfedges[i];
    NewHalfedge.ID=NewHalfedges.size();
    NewHalfedges.push_back(NewHalfedge);
    TransHalfedges[i]=NewHalfedge.ID;
  }
  
   DebugLog<<"ok here 5"<<endl;
  
  Halfedges=NewHalfedges;
  for (int i=0;i<Faces.size();i++)
    Faces[i].AdjHalfedge=TransHalfedges[Faces[i].AdjHalfedge];
  
   DebugLog<<"ok here 6"<<endl;
  
  for (int i=0;i<Vertices.size();i++)
    Vertices[i].AdjHalfedge=TransHalfedges[Vertices[i].AdjHalfedge];
  
   DebugLog<<"ok here 7"<<endl;
  
  for (int i=0;i<Halfedges.size();i++){
    if (Halfedges[i].Twin!=-1)
      Halfedges[i].Twin=TransHalfedges[Halfedges[i].Twin];
    Halfedges[i].Next=TransHalfedges[Halfedges[i].Next];
    Halfedges[i].Prev=TransHalfedges[Halfedges[i].Prev];
  }
  
}


void Mesh::ComputeTwins()
{
  //twinning up edges
  set<TwinFinder> Twinning;
  for (int i=0;i<Halfedges.size();i++){
    if (Halfedges[i].Twin>=0)
      continue;
    
    set<TwinFinder>::iterator Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
    if (Twinit!=Twinning.end()){
      if (Halfedges[Twinit->index].Twin!=-1)
        DebugLog<<"Warning: halfedge "<<Twinit->index<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<endl;
      Halfedges[Twinit->index].Twin=i;
      Halfedges[i].Twin=Twinit->index;
      Twinning.erase(*Twinit);
    } else {
      Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
    }
  }
  
}

void Mesh::WalkBoundary(int &CurrEdge)
{
  do{
    CurrEdge=Halfedges[CurrEdge].Next;
    if (Halfedges[CurrEdge].Twin<0)
      break;  //next boundary over a 2-valence vertex
    CurrEdge=Halfedges[CurrEdge].Twin;
  }while (Halfedges[CurrEdge].Twin>=0);
  
}

//typedef pair< vector<int>, vector<int> > VertexWalk;

typedef adjacency_list <vecS, vecS, undirectedS> Graph;

/*void Mesh::TestUnmatchedTwins()
 {
 vector<int> Untwinned;
 for (int i=0;i<Halfedges.size();i++)
 if ((Halfedges[i].Twin<0)&&(Halfedges[i].Valid))
 Untwinned.push_back(i);
 
 /*for (int i=0;i<Untwinned.size();i++){
 for (int j=0;j<Untwinned.size();j++){
 Vector3D diff1=Vertices[Halfedges[Untwinned[i]].Origin].Coordinates-Vertices[Halfedges[Halfedges[Untwinned[j]].Next].Origin].Coordinates;
 Vector3D diff2=Vertices[Halfedges[Untwinned[j]].Origin].Coordinates-Vertices[Halfedges[Halfedges[Untwinned[i]].Next].Origin].Coordinates;
 if ((sqrt(diff1.squared_length())<10e-6)&&(sqrt(diff2.squared_length())<10e-6)){
 DebugLog<<"Halfedge "<<Untwinned[i]<<":("<<Halfedges[Untwinned[i]].Origin<<","<<Halfedges[Halfedges[Untwinned[i]].Next].Origin<<") is untwinned to ";
 DebugLog<<"Halfedge "<<Untwinned[j]<<":("<<Halfedges[Untwinned[j]].Origin<<","<<Halfedges[Halfedges[Untwinned[j]].Next].Origin<<")\n";
 DebugLog<<Vertices[Halfedges[Untwinned[i]].Origin].Coordinates<<" and "<<Vertices[Halfedges[Untwinned[j]].Origin].Coordinates<<"\n";
 }
 }
 }*/
//}


//currently assuming it's not on the boundary
void Mesh::RemoveVertex(int vindex, std::deque<int>& removeVertexQueue){
  
  DebugLog<<"Trying to removing triangle Vertex "<<vindex<<endl;
  int hebegin = Vertices[vindex].AdjHalfedge;
  int heiterate=hebegin;
  do{
    if (heiterate==-1){
      DebugLog<<"This is a boundary vertex so doing nothing"<<endl;
      return;
    }
    heiterate=Halfedges[Halfedges[heiterate].Prev].Twin;
  }while(heiterate!=hebegin);

  Vertices[vindex].Valid = false;
  
  int remainingFace = Halfedges[hebegin].AdjFace;
  DebugLog<<"Keeping face "<<remainingFace<<endl;
  
  Faces[remainingFace].AdjHalfedge= Halfedges[hebegin].Next;
  heiterate=hebegin;
  int infinityCounter=0;
  do{
    DebugLog<<"Removing halfedge "<<heiterate<<" and its twin "<<Halfedges[heiterate].Twin<<endl;
    int NextEdge = Halfedges[heiterate].Next;
    int PrevEdge = Halfedges[Halfedges[heiterate].Twin].Prev;
    DebugLog<<"Connecting halfedges "<<PrevEdge<<"->"<<NextEdge<<endl;
    Halfedges[NextEdge].Prev = PrevEdge;
    Halfedges[PrevEdge].Next = NextEdge;
    if (Halfedges[NextEdge].AdjFace!=remainingFace){
      DebugLog<<"Invalidating Face "<<Halfedges[NextEdge].AdjFace<<" for NextEdge "<<endl;
      Faces[Halfedges[NextEdge].AdjFace].Valid = false;
    }
    if (Halfedges[PrevEdge].AdjFace!=remainingFace){
      DebugLog<<"Invalidating Face "<<Halfedges[PrevEdge].AdjFace<<" for PrevEdge "<<endl;
      Faces[Halfedges[PrevEdge].AdjFace].Valid = false;
    }
    
    Halfedges[PrevEdge].AdjFace = Halfedges[NextEdge].AdjFace = remainingFace;
    Halfedges[heiterate].Valid = false;
    Halfedges[Halfedges[heiterate].Twin].Valid = false;
    heiterate=Halfedges[Halfedges[heiterate].Prev].Twin;
    infinityCounter++;
    if (infinityCounter>Halfedges.size()){
      DebugLog<<"Created an infinite loop!"<<endl;
      return;
    }
  }while(heiterate!=hebegin);
  
  //cleaning new face
  hebegin =Faces[remainingFace].AdjHalfedge;
  //Faces[remainingFace].NumVertices=0;
  heiterate=hebegin;
  DebugLog<<"New Face edges are "<<endl;
  infinityCounter=0;
  do{
    DebugLog<<heiterate<<",";
    //Faces[remainingFace].NumVertices++;
    Halfedges[heiterate].AdjFace=remainingFace;
    Vertices[Halfedges[heiterate].Origin].AdjHalfedge=heiterate;
    removeVertexQueue.push_front(Halfedges[heiterate].Origin);
    infinityCounter++;
    if (infinityCounter>Halfedges.size()){
      DebugLog<<"Created an infinite loop!"<<endl;
      return;
    }
    heiterate=Halfedges[heiterate].Next;
  }while (heiterate!=hebegin);
  DebugLog<<endl;
  
}

/*void Mesh::RemoveEdge(int heindex)
{
  if ((Halfedges[heindex].Twin>=0)&&(Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices<=2)){
    RemoveFace(Halfedges[Halfedges[heindex].Twin].AdjFace,Halfedges[heindex].Twin);
    return;
  }
  
  DebugLog<<"Removing Edge "<<heindex<<" with Twin "<<Halfedges[heindex].Twin<<" and face "<<Halfedges[heindex].AdjFace<<"\n";
  int face1=Halfedges[heindex].AdjFace;
  DebugLog<<"Face "<<face1<<" vertices before: ";
  int hebegin=Faces[face1].AdjHalfedge;
  int heiterate = hebegin;
  do{
    DebugLog<<Halfedges[heiterate].Origin<<", ";
    heiterate = Halfedges[heiterate].Next;
  }while (heiterate!=hebegin);
  DebugLog<<endl;
  Halfedges[heindex].Valid=false;
  Halfedges[Halfedges[heindex].Next].Prev=Halfedges[heindex].Prev;
  Halfedges[Halfedges[heindex].Prev].Next=Halfedges[heindex].Next;
  DebugLog<<"Connecting halfedges "<<Halfedges[heindex].Prev<<"->"<<Halfedges[heindex].Next<<"\n";
  Vertices[Halfedges[heindex].Origin].AdjHalfedge=Halfedges[heindex].Next;
  DebugLog<<"Vertex "<<Halfedges[heindex].Origin<<"("<<TransVertices[Halfedges[heindex].Origin]<<") points to "<<Halfedges[heindex].Next<<"\n";
  int LeftVertex=Halfedges[heindex].Origin;
  int RemoveVertex=Halfedges[Halfedges[heindex].Next].Origin;
  DebugLog<<"vertex "<<RemoveVertex<<"("<<TransVertices[RemoveVertex]<<")"<<" is removed"<<"\n";
  Vertices[RemoveVertex].Valid=false;
  Halfedges[Halfedges[heindex].Next].Origin=Halfedges[heindex].Origin;
  DebugLog<<"halfedge "<<Halfedges[heindex].Next<<" has vertex"<<Halfedges[heindex].Origin<<"("<<TransVertices[Halfedges[heindex].Origin]<<") as origin\n";
  Faces[Halfedges[heindex].AdjFace].AdjHalfedge=Halfedges[heindex].Next;
  //Faces[Halfedges[heindex].AdjFace].NumVertices--;
  
  
  if (Faces[Halfedges[heindex].AdjFace].NumVertices==0)
    Faces[Halfedges[heindex].AdjFace].Valid=false;
  
  if (Halfedges[heindex].Twin>=0){
    Halfedges[Halfedges[heindex].Twin].Valid=false;
    Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev=Halfedges[Halfedges[heindex].Twin].Prev;
    Halfedges[Halfedges[Halfedges[heindex].Twin].Prev].Next=Halfedges[Halfedges[heindex].Twin].Next;
    DebugLog<<"Connecting "<<Halfedges[Halfedges[heindex].Twin].Prev<<"->"<<Halfedges[Halfedges[heindex].Twin].Next<<"\n";
    //Vertices[Halfedges[Halfedges[heindex].Twin].Origin].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
    Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Origin=LeftVertex;
    Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;
    Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices--;
    if (Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices==0)
      Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].Valid=false;
  }
  
  for (int i=0;i<Halfedges.size();i++)
    if (Halfedges[i].Origin==RemoveVertex){
      Halfedges[i].Origin=LeftVertex;
      DebugLog<<"Other halfedge "<<i<<" has vertex"<<LeftVertex<<"("<<TransVertices[LeftVertex]<<") as origin\n";
    }
  
  if (Faces[face1].NumVertices>0){
    DebugLog<<"Face "<<face1<<" vertices after: ";
    int hebegin=Faces[Halfedges[heindex].AdjFace].AdjHalfedge;
    int heiterate = hebegin;
    do{
      DebugLog<<Halfedges[heiterate].Origin<<", ";
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
  } else
    DebugLog<<"Face "<<face1<<" now has zero vertices"<<endl;
  DebugLog<<"Done removing edge"<<endl;
  
}*/

/*void Mesh::RemoveFace(int findex, int heindex)
{
  //return;
  DebugLog<<"Removing Face "<<findex<<" with initial edge "<<heindex<<" and # vertices "<<Faces[findex].NumVertices<<"\n";
  
  DebugLog<<"Face state before: \n";
  int hebegin=heindex;
  int heiterate=hebegin;
  do{
    DebugLog<<"Halfedge "<<heiterate<<" with origin "<<Halfedges[heiterate].Origin<<"("<<TransVertices[Halfedges[heiterate].Origin]<<")\n";
    heiterate=Halfedges[heiterate].Next;
  }while (heiterate!=hebegin);
  
  int LeftVertex=Halfedges[heindex].Origin;
  
  DebugLog<<"Leaving Vertex "<<LeftVertex<<"\n";
  Faces[findex].Valid=false;
  vector<int> ReplaceOrigins;
  do{
    DebugLog<<"Removing halfedge "<<heiterate<<" with twin "<<Halfedges[heiterate].Twin<<" and origin "<<Halfedges[heiterate].Origin<<"\n";
    
    Halfedges[heiterate].Valid=false;
    if (Halfedges[heiterate].Origin!=LeftVertex){
      Vertices[Halfedges[heiterate].Origin].Valid=false;
      ReplaceOrigins.push_back(Halfedges[heiterate].Origin);
    }
    if (Halfedges[heiterate].Twin<0){
      heiterate=Halfedges[heiterate].Next;
      continue;
    }
    int ReduceEdge=Halfedges[heiterate].Twin;
    Halfedges[ReduceEdge].Valid=false;
    Halfedges[Halfedges[ReduceEdge].Next].Prev=Halfedges[ReduceEdge].Prev;
    Halfedges[Halfedges[ReduceEdge].Prev].Next=Halfedges[ReduceEdge].Next;
    Faces[Halfedges[ReduceEdge].AdjFace].NumVertices--;
    
    DebugLog<<"Connecting "<<Halfedges[ReduceEdge].Prev<<"->"<<Halfedges[ReduceEdge].Next<<"\n";
    
    Faces[Halfedges[ReduceEdge].AdjFace].AdjHalfedge=Halfedges[ReduceEdge].Next;
    Vertices[LeftVertex].AdjHalfedge=Halfedges[ReduceEdge].Next;
    DebugLog<<"Vertex "<<LeftVertex<<"("<<TransVertices[LeftVertex]<<") now points to "<<Halfedges[ReduceEdge].Next<<"\n";
    
    heiterate=Halfedges[heiterate].Next;
  }while (heiterate!=hebegin);
  
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    for (int j=0;j<ReplaceOrigins.size();j++)
      if (Halfedges[i].Origin==ReplaceOrigins[j]){
        DebugLog<<"Now halfedge "<<i<<" has vertex "<<LeftVertex<<"("<<TransVertices[LeftVertex]<<") as origin instead of "<<ReplaceOrigins[j]<<"("<<TransVertices[ReplaceOrigins[j]]<<")\n";
        Halfedges[i].Origin=LeftVertex;
        Vertices[LeftVertex].AdjHalfedge=i;
      }
  }
  
  DebugLog<<"Finished Removing face\n";
}

void Mesh::RemoveDegree2Faces()
{
  for (int i=0;i<Faces.size();i++){
    if ((Faces[i].NumVertices!=2)||(!Faces[i].Valid))
      continue;
    
    Faces[i].Valid=false;
    int he1 = Faces[i].AdjHalfedge;
    int he2 = Halfedges[he1].Next;
    
    DebugLog<<"Removing Face "<<i<<" with halfedges "<<he1<<" and "<<he2<<endl;
    DebugLog<<"Twinning halfedges "<<Halfedges[he1].Twin<<" and "<<Halfedges[he2].Twin<<endl;
    
    Halfedges[he1].Valid=false;
    Halfedges[he2].Valid=false;
    
    if (Halfedges[he1].Twin!=-1)
      Halfedges[Halfedges[he1].Twin].Twin = Halfedges[he2].Twin;
    
    if (Halfedges[he2].Twin!=-1)
      Halfedges[Halfedges[he2].Twin].Twin = Halfedges[he1].Twin;
    
    if (Halfedges[Halfedges[he1].Next].Twin!=-1)
      Vertices[Halfedges[he1].Origin].AdjHalfedge=Halfedges[Halfedges[he1].Next].Twin;
    else if (Halfedges[he1].Twin!=-1)
      Vertices[Halfedges[he1].Origin].AdjHalfedge=Halfedges[Halfedges[he1].Twin].Next;
    else //both are -1 which means isolated vertex
      Vertices[Halfedges[he1].Origin].Valid=false;
    
    if (Halfedges[Halfedges[he2].Next].Twin!=-1)
      Vertices[Halfedges[he2].Origin].AdjHalfedge=Halfedges[Halfedges[he2].Next].Twin;
    else if (Halfedges[he2].Twin!=-1)
      Vertices[Halfedges[he2].Origin].AdjHalfedge=Halfedges[Halfedges[he2].Twin].Next;
    else //both are -1 which means isolated vertex
      Vertices[Halfedges[he2].Origin].Valid=false;
    
    if ((Halfedges[Halfedges[he1].Next].Twin!=-1)&&(Halfedges[Halfedges[he1].Next].Twin!=-1))
    {
      int spikeLeadingHE=-1;
      if (Halfedges[Halfedges[he2].Twin].Next == Halfedges[he1].Twin)
        spikeLeadingHE = Halfedges[he1].Twin;
      else if (Halfedges[Halfedges[he1].Twin].Next == Halfedges[he2].Twin)
        spikeLeadingHE =Halfedges[he2].Twin;
      
      if (spikeLeadingHE==-1)
        continue;
      
      DebugLog<<"Removing spike "<<Halfedges[spikeLeadingHE].Prev<<"->"<<spikeLeadingHE<<endl;
      
      //removing spike
      Vertices[Halfedges[spikeLeadingHE].Origin].Valid=false;
      Halfedges[spikeLeadingHE].Valid=false;
      Halfedges[Halfedges[spikeLeadingHE].Prev].Valid=false;
      
      
      Halfedges[Halfedges[Halfedges[spikeLeadingHE].Prev].Prev].Next = Halfedges[spikeLeadingHE].Next;
      Halfedges[Halfedges[spikeLeadingHE].Next].Prev =Halfedges[Halfedges[spikeLeadingHE].Prev].Prev;
      Vertices[Halfedges[Halfedges[spikeLeadingHE].Prev].Origin].AdjHalfedge=Halfedges[spikeLeadingHE].Next;
      Faces[Halfedges[spikeLeadingHE].AdjFace].NumVertices-=2;  //the base vertex for the spike appears twice
      
      DebugLog<<"Connecting halfedges "<<Halfedges[Halfedges[spikeLeadingHE].Prev].Prev<<"->"<<Halfedges[spikeLeadingHE].Next<<endl;
      DebugLog<<"Vertex "<<Halfedges[Halfedges[spikeLeadingHE].Prev].Origin<<" now points to halfedge "<<Halfedges[spikeLeadingHE].Next<<endl;
      DebugLog<<"Face "<<Halfedges[spikeLeadingHE].AdjFace<<" now has "<<Faces[Halfedges[spikeLeadingHE].AdjFace].NumVertices<<" vertices"<<endl;
    }
    
    //checking if twinned edges are not consecutive. If they are, remove them as this is a spike
    
    
    if (!CheckMesh(false, false))
      return;
  }
}*/



bool Mesh::SimplifyHexMesh(int N)
{
  //unifying vertices which are similar
  
  DebugLog.open("Debugging.txt");
  
  if (!CheckMesh(false, false, false))
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
    DebugLog<<"Walking original triangle boundary"<<endl;
    do{
      visitedOrig[Halfedges[heiterate].OrigHalfedge]=true;
      DebugLog<<"Walking boundary halfedge "<<heiterate<<" with vertex "<<Halfedges[heiterate].Origin<<" on original halfedge "<<Halfedges[heiterate].OrigHalfedge<<endl;
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
      CurrMatches=FindVertexMatch(DebugLog, PointSet1, PointSet2);
    
    for (int j=0;j<CurrMatches.size();j++){
      CurrMatches[j].first =VertexSets1[i][CurrMatches[j].first];
      CurrMatches[j].second=VertexSets2[i][CurrMatches[j].second];
      DebugLog<<"Vertex "<<CurrMatches[j].first<<" is matched with vertex "<<CurrMatches[j].second<<endl;
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
  
  cout<<"Max matching distance: "<<MaxDist<<endl;
  
  //vector<int> TransVertices(Vertices.size());
  TransVertices.resize(Vertices.size());
  int NumNewVertices = connected_components(MatchGraph, &TransVertices[0]);
  for (int i=0;i<NumNewVertices;i++){
    DebugLog<<"TransVertices group "<<i<<": ";
    for (int j=0;j<TransVertices.size();j++)
      if (TransVertices[j]==i)
        DebugLog<<j<<", ";
    DebugLog<<endl;
  }
  
  for (int i=0;i<Faces.size();i++){
    int hebegin = Faces[i].AdjHalfedge;
    int heiterate = hebegin;
    DebugLog<<"Face "<<i<<" is initially of ";
    do{
      DebugLog<<heiterate<<"["<<Halfedges[heiterate].Origin<<"], ";
      heiterate = Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
    DebugLog<<endl;
  }
  
  for (int i=0;i<Halfedges.size();i++){
    DebugLog<<"Halfedge "<<i<<" is initially a ";
    if (Halfedges[i].isHex)
      DebugLog<<"Hex edge with ";
    else
      DebugLog<<"Triangle edge with ";
    
    DebugLog<<"Origin: "<<Halfedges[i].Origin<<"("<<TransVertices[Halfedges[i].Origin]<<")\n";
    DebugLog<<"Prev: "<<Halfedges[i].Prev<<"\n";
    DebugLog<<"Next: "<<Halfedges[i].Next<<"\n";
    DebugLog<<"Twin: "<<Halfedges[i].Twin<<"\n";
    DebugLog<<"Face: "<<Halfedges[i].AdjFace<<"\n";
  }
  
  //seeing if there are faces that are going to degenerate
  /*vector<int> TransVertices2(NewVertices.size());
   for (int i=0;i<NewVertices.size();i++)
   TransVertices2[i]=i;*/
  
  //adding other vertex to the degeneration if needed
  /*bool ThereisChange;
   do{
   ThereisChange=false;
   for (int i=0;i<Halfedges.size();i++){
   if ((!Halfedges[i].Valid)||(Halfedges[i].Twin>=0))
   continue;
   
   if (TransVertices[Halfedges[i].Origin]!=TransVertices[Halfedges[Halfedges[i].Next].Origin])
   continue;  //this edge is not to be collapsed
   
   if (Faces[Halfedges[i].AdjFace].NumVertices<=3) {
   for (int j=0;j<Faces[Halfedges[i].AdjFace].NumVertices;j++){
   if (TransVertices[Faces[Halfedges[i].AdjFace].Vertices[j]]!=TransVertices[Halfedges[i].Origin]){
   add_edge(Faces[Halfedges[i].AdjFace].Vertices[j], Halfedges[i].Origin, MatchGraph);
   ThereisChange=true;
   }
   }
   }
   }
   
   if (ThereisChange){
   TransVertices.clear();
   TransVertices.resize(Vertices.size());
   NumNewVertices = connected_components(MatchGraph, &TransVertices[0]);
   }
   }while (ThereisChange);*/
  
  if (!CheckMesh(false, false, false))
    return false;
  
  //removing edges (and consequent faces) which will degenerate
  /*for (int i=0;i<Halfedges.size();i++){
    if ((!Halfedges[i].Valid)||(Halfedges[i].Twin>=0))
      continue;
    if (TransVertices[Halfedges[i].Origin]!=TransVertices[Halfedges[Halfedges[i].Next].Origin])
      continue;  //this edge is OK
    
    if (Faces[Halfedges[i].AdjFace].NumVertices<=2)
      RemoveFace(Halfedges[i].AdjFace, i);
    else
      RemoveEdge(i);
    
    if (!CheckMesh(false, false))
      return false;
  }
  
  //Removing degree 2 faces
  RemoveDegree2Faces();
  
  if (!CheckMesh(false, false))
    return false;*/
  
  
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
    if (!transClaimed[i]){
      DebugLog<<"TransVertex "<<i<<" not claimed!"<<endl;
      DebugLog<<"Group of vertices: ";
      for (int j=0;j<TransVertices.size();j++)
        if (TransVertices[j]==i)
          DebugLog<<j<<", ";
      DebugLog<<endl;
      //return false;
      NewVertices[i].Valid=false;  //this vertex is dead to begin with
    }
  
  Vertices=NewVertices;
  /*for (int i=0;i<Faces.size();i++){
    if (!Faces[i].Valid)
      continue;
    for (int j=0;j<Faces[i].NumVertices;j++)
      Faces[i].Vertices[j]=TransVertices[Faces[i].Vertices[j]];
  }*/
  
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    Halfedges[i].Origin=TransVertices[Halfedges[i].Origin];
    Vertices[Halfedges[i].Origin].AdjHalfedge=i;
  }
  
  //sanity check: that every halfedge has a single potential twin (only for closed meshes)
  /* set<TwinFinder> checkTwinning;
   for (int i=0;i<Halfedges.size();i++){
   if (!Halfedges[i].Valid)
   continue;
   
   set<TwinFinder>::iterator Twinit=checkTwinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
   if (Twinit!=checkTwinning.end())
   checkTwinning.erase(*Twinit);
   else
   checkTwinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
   
   }
   
   //checkTwinning should be empty
   for (set<TwinFinder>::iterator it = checkTwinning.begin();it!=checkTwinning.end();it++)
   DebugLog<<"Halfedge "<<it->index<<" with vertices " <<it->v1<<","<<it->v2<<" doesn't have a twin."<<endl;
   
   if (!checkTwinning.empty())
   return false;*/
  
  
  if (!CheckMesh(true, false, false))
    return false;
  
  //twinning up edges
  set<TwinFinder> Twinning;
  for (int i=0;i<Halfedges.size();i++){
    if ((Halfedges[i].Twin>=0)||(!Halfedges[i].Valid))
      continue;
    
    set<TwinFinder>::iterator Twinit=Twinning.find(TwinFinder(0,Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
    if (Twinit!=Twinning.end()){
      DebugLog<<"Twinning halfedge "<<i<<" and halfedge "<<Twinit->index<<"\n";
      if (Halfedges[Twinit->index].Twin!=-1)
        DebugLog<<"warning: halfedge "<<Twinit->index<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<endl;
      if (Halfedges[i].Twin!=-1)
        DebugLog<<"warning: halfedge "<<i<<" is already twinned to halfedge "<<Halfedges[Twinit->index].Twin<<endl;
      Halfedges[Twinit->index].Twin=i;
      Halfedges[i].Twin=Twinit->index;
      
      if (Halfedges[i].isHex){
        DebugLog<<"Halfedge "<<i<<" is hex, infecting the other\n";
        Halfedges[Twinit->index].isHex = true;
      } else if (Halfedges[Twinit->index].isHex){
        DebugLog<<"Halfedge "<<Twinit->index<<" is hex, infecting the other\n";
        Halfedges[i].isHex = true;
      }
      Twinning.erase(*Twinit);
    } else {
      Twinning.insert(TwinFinder(i,Halfedges[i].Origin,Halfedges[Halfedges[i].Next].Origin));
    }
  }
  
  //check if there are any non-twinned edge which shouldn't be in a closed mesh
  for (int i=0;i<Halfedges.size();i++){
    if (Halfedges[i].Twin==-1)
      DebugLog<<"Halfedge "<<i<<" does not have a twin!"<<endl;
  }
  

  
  for (int i=0;i<Halfedges.size();i++){
    if (!Halfedges[i].Valid)
      continue;
    if (Halfedges[i].isHex)
      DebugLog<<"Hex edge "<<i<<"\n";
    else
      DebugLog<<"Triangle edge "<<i<<"\n";
    
    DebugLog<<"Origin: "<<Halfedges[i].Origin<<"\n";
    DebugLog<<"Prev: "<<Halfedges[i].Prev<<"\n";
    DebugLog<<"Next: "<<Halfedges[i].Next<<"\n";
    DebugLog<<"Twin: "<<Halfedges[i].Twin<<"\n";
    DebugLog<<"Face: "<<Halfedges[i].AdjFace<<"\n";
  }
  
  if (!CheckMesh(true, true, true))
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
    if ((Halfedges[i].isHex)&&(Halfedges[i].Valid)){
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
  DebugLog<<"Realigning edges around vertices"<<endl;
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
      if ((Halfedges[heiterate].isHex)||(Halfedges[heiterate].Twin==-1))
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
  DebugLog<<"Invalidating triangle vertices and edges"<<endl;
  for (int i=0;i<Vertices.size();i++)
    if (isPureTriangle[i])
      Vertices[i].Valid=false;
  
  for (int i=0;i<Halfedges.size();i++)
    if ((!Halfedges[i].isHex)&&(Halfedges[i].Twin!=-1))
      Halfedges[i].Valid=false;
  
  //realigning faces
  DebugLog<<"Realigning faces"<<endl;
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
        DebugLog<<"Infinity loop in realigning faces on halfedge "<<i<<endl;
        return false;
      }
      Halfedges[heiterate].AdjFace=currFace;
      heiterate=Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
  }
  
  int countThree=0;
  DebugLog<<"Invalidating remainder faces"<<endl;
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
     
  DebugLog<<"Invalidating ear (latent valence 2) faces"<<endl;
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
      DebugLog<<"Invalidating ear Face "<<i<<endl;
      DebugLog<<"Its vertices are "<<endl;
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
    
  
 /* std::deque<int> removeVertexQueue;
  for (int i=0;i<Vertices.size();i++){
    if ((isPureTriangle[i])&&(Vertices[i].Valid)){
      removeVertexQueue.push_back(i);
    }
  }
  
  while (!removeVertexQueue.empty()){
    int currVertex=removeVertexQueue.front();
    removeVertexQueue.pop_front();
    if (!Vertices[currVertex].Valid)
      continue;
    
    cout<<"Removing vertex "<<currVertex<<endl;
    RemoveVertex(currVertex, removeVertexQueue);
    if (!CheckMesh(false, true, true))
        return false;
    }
  
  cout<<"Done removing vertices!"<<endl;*/
  
  if (!CheckMesh(true, true, true))
    return false;
  
  
  
  
  /*std::queue<int> removeVertices;
  
  //removing the rest of triangle halfedges
  std::queue<int> removeHalfedges;
  for (int i=0;i<Halfedges.size();i++)
    if ((!Halfedges[i].isHex)&&(Halfedges[i].Valid))
      removeHalfedges.push(i);
  
  while (!removeHalfedges.empty()){
    int currhe = removeHalfedges.front();
    removeHalfedges.pop();
    if (!Halfedges[currhe].Valid)
      continue;
    if (!JoinFace(currhe))
      removeHalfedges.push(currhe);
  }*/
  
  
  //if (!CheckMesh(true, true, true))
   // return false;
  
  //REMOVE!!!!!
  /*cout<<"Cleaning mesh cheatingly!"<<endl;
  CleanMesh();
  cout<<"Mesh is cheatingly clean!"<<endl;
  return true;*/
  
  //unifying chains of edges
  

  
  for (int i=0;i<Valences.size();i++)
    if ((Vertices[i].Valid)&&(Valences[i]<2))
      Vertices[i].Valid=false;
  
  DebugLog<<"Starting unifying edges" <<endl;
  for (int i=0;i<Vertices.size();i++){
    
    if ((Vertices[i].Valid)&&(Valences[i]<=2)&&(!isEar[i]))
      UnifyEdges(Vertices[i].AdjHalfedge);
    
      /*if (Halfedges[Vertices[i].AdjHalfedge].Twin<0)
        if (!CheckMesh(true, true, true))
              return false;*/
    
  }
  
  
  /*if (!CheckMesh(true, true, true))
    return false;
  
  //re-assigning vertices on every face according to halfedges
  /*for (int i=0;i<Faces.size();i++)
  {
    if (!Faces[i].Valid){
      Faces[i].NumVertices=0;
      continue;
    }
    int hebegin=Faces[i].AdjHalfedge;
    int heiterate=hebegin;
    Faces[i].NumVertices=0;
    do{
      Faces[i].Vertices[Faces[i].NumVertices]=Halfedges[heiterate].Origin;
      Faces[i].NumVertices++;
      heiterate=Halfedges[heiterate].Next;
    }while (heiterate!=hebegin);
  }*/
  
  if (!CheckMesh(true, true, true))
    return false;
  
  //remove non-valid components
  CleanMesh();
  
  //checking if mesh is valid
  if (!CheckMesh(true, true, true))
    return false;
  
  //computing centroids
  /*for (int i=0;i<Faces.size();i++){
    Faces[i].Centroid=Point3D(0.0,0.0,0.0);
    for (int j=0;j<Faces[i].NumVertices;j++){
      Faces[i].Centroid=Faces[i].Centroid+(Vertices[Faces[i].Vertices[j]].Coordinates-CGAL::ORIGIN);
    }
    Faces[i].Centroid=CGAL::ORIGIN+(Faces[i].Centroid-CGAL::ORIGIN)/(double)Faces[i].NumVertices;
  }*/
  
  //completing angle diffs - currently not working with a boundary
  /*for (int i=0;i<Vertices.size();i++){
    int hebegin = Vertices[i].AdjHalfedge;
    //weeding out boundaries - WHAT TO DO WITH THEM?
    int heiterate = hebegin;
    bool isBoundary=false;
    do{
      if (Halfedges[heiterate].Twin==-1){
        isBoundary=true;
        break;
      }
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);
    
    if (isBoundary)
      continue;
    
    heiterate = hebegin;
    double sumPrescribedDiff=0;
    int missingPrescribed=0;
    do{
      if (Halfedges[heiterate].prescribedAngle>0.0)
        sumPrescribedDiff+=Halfedges[heiterate].prescribedAngle;
      else
        missingPrescribed++;
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);
    
    //TOMMOROW - USE ACTUAL DOUBLE ANGLES BC IT"S NOT POSSIBLE OTHERWISE
    //sanity check
    //if (missingPrescribed==0)
      //cout<<"sumPrescribedDiff: "<<sumPrescribedDiff<<endl;
    
    double completionAngle = 2*igl::PI-sumPrescribedDiff;
    do{
      if (Halfedges[heiterate].prescribedAngle<=0)
        //Halfedges[heiterate].prescribedAngle=2*igl::PI*(double)Halfedges[heiterate].prescribedAngleDiff/N;*/
      //else
        /*Halfedges[heiterate].prescribedAngle = completionAngle;
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);
    
    //sanity check
   /* double shouldBeTwoPi=0;
    do{
      shouldBeTwoPi+=Halfedges[heiterate].prescribedAngle;
      heiterate = Halfedges[Halfedges[heiterate].Twin].Next;
    }while(heiterate!=hebegin);
    cout<<"shouldBeTwoPi: "<<shouldBeTwoPi<<endl;*/
    
  //}
  
  
  return true;
  
}


