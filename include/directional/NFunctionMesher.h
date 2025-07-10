//
// Created by Amir Vaxman on 20.04.24.
//

#ifndef DIRECTIONAL_N_FUNCTION_MESHER
#define DIRECTIONAL_N_FUNCTION_MESHER

#include <set>
#include <math.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <directional/exact_geometric_definitions.h>
#include <directional/dcel.h>
#include <directional/setup_mesher.h>

namespace directional{

class NFunctionMesher {
public:
    
    const TriMesh& origMesh;
    const MesherData& mData;
    
    struct SegmentData{
        bool isFunction;
        int origHalfedge;
        int origNFunctionIndex;  //the original parameteric function assoicated with this edge
        int lineInPencil;
        std::set<ENumber> intParams;
        //double prescribedAngle;  //the actual prescribed angle
        
        SegmentData():isFunction(false), origHalfedge(-1), origNFunctionIndex(-1),  lineInPencil(-1), intParams(){} //prescribedAngle(-1.0){}
    };
    
    struct VData{
        Eigen::RowVector3d coords;
        EVector3 eCoords;
    };
    
    typedef DCEL<VData,  SegmentData, bool, bool> FunctionDCEL;
    FunctionDCEL genDcel;
    
    //halfedge quantities
    Eigen::MatrixXd NFunction;
    std::vector<std::vector<ENumber>> exactNFunction;
    
    
    //mesh generation functions found in generate_mesh.h
    void arrange_on_triangle(const std::vector<EVector2>& triangle,
                             const std::vector<std::pair<int, bool>>& triangleData,
                             const std::vector<LinePencil>& linePencils,
                             const std::vector<int>& linePencilData,
                             std::vector<EVector2>& V,
                             FunctionDCEL& triDcel);
    
    void segment_arrangement(const std::vector<Segment2>& segments,
                             const std::vector<SegmentData>& data,
                             const Eigen::Matrix<ENumber, Eigen::Dynamic ,2> I2dts,
                             const Eigen::Matrix<ENumber, Eigen::Dynamic ,1> t00s,
                             std::vector<EVector2>& V,
                             FunctionDCEL& triDcel);
    
    void generate_mesh(const unsigned long Resolution);
    
    std::vector<int> TransVertices;
    std::vector<int> InStrip;
    std::vector<std::set<int> > VertexChains;
    
    
    struct MergeData {
        const bool operator()(const int &v1, const int &v2) const { return v1; }
    };
    
    
    void TestUnmatchedTwins();
    
    struct PointPair{
        int Index1, Index2;
        ENumber Distance;
        
        PointPair(int i1, int i2, const EVector3& d):Index1(i1), Index2(i2){
            Distance=d.max_abs();
        }
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
    
    std::vector<std::pair<int,int>> FindVertexMatch(const bool verbose, std::vector<EVector3>& Set1, std::vector<EVector3>& Set2)
    {
        std::set<PointPair> PairSet;
        for (int i=0;i<Set1.size();i++)
            for (int j=0;j<Set2.size();j++)
                PairSet.insert(PointPair(i,j,Set1[i]-Set2[j]));
        
        assert (Set1.size()==Set2.size() && "NFunctionMesher::FindVertexMatch(): The two sets are of different sizes!! ");
        
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
                if (squaredDistance(Set1[Result[i].first],Set2[Result[i].second])>ENumber(0)){
                    std::cout<<"("<<Result[i].first<<","<<Result[i].second<<") with dist "<<squaredDistance(Set1[Result[i].first],Set2[Result[i].second]).to_double()<<std::endl;
                    std::cout<<"Distance is abnormally not zero!"<<std::endl;
                }
            }
        }
        
        
        return Result;
        
    }
    
    bool simplify_mesh(){
        //unifying vertices which are similar
        
        using namespace std;
        using namespace Eigen;
        
        if (!genDcel.check_consistency(mData.verbose, false, false, false))
            return false;
        
        int MaxOrigHE=-3276700.0;
        for (int i=0;i<genDcel.halfedges.size();i++)
            MaxOrigHE=std::max(MaxOrigHE, genDcel.halfedges[i].data.origHalfedge);
        
        vector<bool> visitedOrig(MaxOrigHE+1);
        for (int i=0;i<MaxOrigHE+1;i++) visitedOrig[i]=false;
        for (int i=0;i<genDcel.halfedges.size();i++){
            if (genDcel.halfedges[i].data.origHalfedge<0)
                continue;
            if (visitedOrig[genDcel.halfedges[i].data.origHalfedge])
                continue;
            
            int hebegin = i;
            int heiterate = hebegin;
            do{
                visitedOrig[genDcel.halfedges[heiterate].data.origHalfedge]=true;
                genDcel.walk_boundary(heiterate);
            }while (heiterate!=hebegin);
            
        }
        
        vector< vector<int> > BoundEdgeCollect1(MaxOrigHE+1);
        vector< vector<int> > BoundEdgeCollect2(MaxOrigHE+1);
        vector<bool> Marked(genDcel.halfedges.size());
        for (int i=0;i<genDcel.halfedges.size();i++) Marked[i]=false;
        //finding out vertex correspondence along twin edges of the original mesh by walking on boundaries
        for (int i=0;i<genDcel.halfedges.size();i++){
            if ((genDcel.halfedges[i].data.origHalfedge<0)||(Marked[i]))
                continue;
            
            //find the next beginning of a boundary
            int PrevOrig;
            int CurrEdge=i;
            do{
                PrevOrig=genDcel.halfedges[CurrEdge].data.origHalfedge;
                genDcel.walk_boundary(CurrEdge);
            }while(PrevOrig==genDcel.halfedges[CurrEdge].data.origHalfedge);
            
            //filling out strips of boundary with the respective attached original halfedges
            int BeginEdge=CurrEdge;
            vector<pair<int,int> > CurrEdgeCollect;
            do{
                CurrEdgeCollect.push_back(pair<int, int> (genDcel.halfedges[CurrEdge].data.origHalfedge, CurrEdge));
                Marked[CurrEdge]=true;
                genDcel.walk_boundary(CurrEdge);
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
                VertexSets1[i].push_back(genDcel.halfedges[BoundEdgeCollect1[i][j]].vertex);
            
            if (BoundEdgeCollect1[i].size()>0)
                VertexSets1[i].push_back(genDcel.halfedges[genDcel.halfedges[BoundEdgeCollect1[i][BoundEdgeCollect1[i].size()-1]].next].vertex);
            
            for (int j=0;j<BoundEdgeCollect2[i].size();j++)
                VertexSets2[i].push_back(genDcel.halfedges[BoundEdgeCollect2[i][j]].vertex);
            
            if (BoundEdgeCollect2[i].size()>0)
                VertexSets2[i].push_back(genDcel.halfedges[genDcel.halfedges[BoundEdgeCollect2[i][BoundEdgeCollect2[i].size()-1]].next].vertex);
            
            std::reverse(VertexSets2[i].begin(),VertexSets2[i].end());
        }
        
        //finding out vertex matches
        vector<pair<int, int> > VertexMatches;
        for (int i=0;i<MaxOrigHE+1;i++){
            vector<EVector3> PointSet1(VertexSets1[i].size());
            vector<EVector3> PointSet2(VertexSets2[i].size());
            for (int j=0;j<PointSet1.size();j++)
                PointSet1[j]=genDcel.vertices[VertexSets1[i][j]].data.eCoords;
            
            for (int j=0;j<PointSet2.size();j++)
                PointSet2[j]=genDcel.vertices[VertexSets2[i][j]].data.eCoords;
            
            vector<pair<int, int> > CurrMatches;
            if ((!PointSet1.empty())&&(!PointSet2.empty()))
                CurrMatches=FindVertexMatch(mData.verbose, PointSet1, PointSet2);
            
            for (int j=0;j<CurrMatches.size();j++){
                CurrMatches[j].first =VertexSets1[i][CurrMatches[j].first];
                CurrMatches[j].second=VertexSets2[i][CurrMatches[j].second];
            }
            
            VertexMatches.insert(VertexMatches.end(), CurrMatches.begin(), CurrMatches.end() );
        }
        
        //finding connected components, and uniting every component into a random single vertex in it (it comes out the last mentioned)
        /*Graph MatchGraph;
         for (int i=0;i<vertices.size();i++)
         add_vertex(MatchGraph);
         for (int i=0;i<VertexMatches.size();i++)
         add_edge(VertexMatches[i].first, VertexMatches[i].second, MatchGraph);*/
        
        double MaxDist=-327670000.0;
        for (int i=0;i<VertexMatches.size();i++)
            MaxDist=std::max(MaxDist, (genDcel.vertices[VertexMatches[i].first].data.coords-genDcel.vertices[VertexMatches[i].second].data.coords).squaredNorm());
        
        if (mData.verbose)
            std::cout<<"Max matching distance: "<<MaxDist<<endl;
        
        //vector<int> Transvertices(vertices.size());
        TransVertices.resize(genDcel.vertices.size());
        int NumNewVertices = connectedComponents(VertexMatches, TransVertices);
        
        if (!genDcel.check_consistency(mData.verbose, false, false, false))
            return false;
        
        vector<bool> transClaimed(NumNewVertices);
        for (int i=0;i<NumNewVertices;i++)
            transClaimed[i]=false;
        //unifying all vertices into the TransVertices
        vector<FunctionDCEL::Vertex> NewVertices(NumNewVertices);
        for (int i=0;i<genDcel.vertices.size();i++){  //redundant, but not terrible
            if (!genDcel.vertices[i].valid)
                continue;
            FunctionDCEL::Vertex NewVertex=genDcel.vertices[i];
            NewVertex.ID=TransVertices[i];
            transClaimed[TransVertices[i]]=true;
            NewVertices[TransVertices[i]]=NewVertex;
        }
        
        for (int i=0;i<NumNewVertices;i++)
            if (!transClaimed[i])
                NewVertices[i].valid=false;  //this vertex is dead to begin with
        
        genDcel.vertices=NewVertices;
        
        for (int i=0;i<genDcel.halfedges.size();i++){
            if (!genDcel.halfedges[i].valid)
                continue;
            genDcel.halfedges[i].vertex=TransVertices[genDcel.halfedges[i].vertex];
            genDcel.vertices[genDcel.halfedges[i].vertex].halfedge=i;
        }
        
        
        if (!genDcel.check_consistency(mData.verbose, true, false, false))
            return false;
        
        //twinning up halfedges
        set<FunctionDCEL::TwinFinder> Twinning;
        for (int i=0;i<genDcel.halfedges.size();i++){
            if ((genDcel.halfedges[i].twin>=0)||(!genDcel.halfedges[i].valid))
                continue;
            
            set<FunctionDCEL::TwinFinder>::iterator Twinit=Twinning.find(FunctionDCEL::TwinFinder(0,genDcel.halfedges[genDcel.halfedges[i].next].vertex,
                                                                                                  genDcel.halfedges[i].vertex));
            if (Twinit!=Twinning.end()){
                if ((genDcel.halfedges[Twinit->index].twin!=-1)&&(mData.verbose))
                    std::cout<<"warning: halfedge "<<Twinit->index<<" is already twinned to halfedge "<<genDcel.halfedges[Twinit->index].twin<<std::endl;
                if ((genDcel.halfedges[i].twin!=-1)&&(mData.verbose))
                    std::cout<<"warning: halfedge "<<i<<" is already twinned to halfedge "<<genDcel.halfedges[Twinit->index].twin<<std::endl;
                genDcel.halfedges[Twinit->index].twin=i;
                genDcel.halfedges[i].twin=Twinit->index;
                
                //assigning a single edge out of them and invalidaing the other one.
                genDcel.edges[genDcel.halfedges[Twinit->index].edge].valid=false;
                genDcel.halfedges[Twinit->index].edge = genDcel.halfedges[i].edge;
                
                //std::cout<<"Twinning halfedge "<<i<<" to halfedge "<<Twinit->index<<std::endl;
                
                if (genDcel.halfedges[i].data.isFunction){
                    genDcel.halfedges[Twinit->index].data.isFunction = true;
                } else if (genDcel.halfedges[Twinit->index].data.isFunction){
                    genDcel.halfedges[i].data.isFunction = true;
                }
                Twinning.erase(*Twinit);
            } else {
                Twinning.insert(FunctionDCEL::TwinFinder(i,genDcel.halfedges[i].vertex,genDcel.halfedges[genDcel.halfedges[i].next].vertex));
            }
        }
        
        //check if there are any non-twinned edge which shouldn't be in a closed mesh
        /*if (verbose){
         for (int i=0;i<Halfedges.size();i++){
         if (Halfedges[i].twin==-1)
         std::cout<<"Halfedge "<<i<<" does not have a twin!"<<std::endl;
         }
         }*/
        
        
        if (!genDcel.check_consistency(mData.verbose, true, true, true))
            return false;
        
        //removing triangle components
        
        //starting with pure triangle vertices
        std::vector<bool> isPureTriangle(genDcel.vertices.size());
        std::vector<bool> isBoundary(genDcel.vertices.size());
        for (int i=0;i<genDcel.vertices.size();i++){
            isPureTriangle[i]=true;
            isBoundary[i]=false;
        }
        for (int i=0;i<genDcel.halfedges.size();i++){
            if ((genDcel.halfedges[i].data.isFunction)&&(genDcel.halfedges[i].valid)){
                isPureTriangle[genDcel.halfedges[i].vertex]=isPureTriangle[genDcel.halfedges[genDcel.halfedges[i].next].vertex]=false;  //adjacent to at least one hex edge
            }
            if (genDcel.halfedges[i].twin==-1){
                isBoundary[genDcel.halfedges[i].vertex]=true;
                isPureTriangle[genDcel.halfedges[i].vertex]=false;  //this shouldn't be removed
            }
        }
        
        std::vector<bool> isEar(genDcel.vertices.size());
        for (int i=0;i<genDcel.vertices.size();i++){
            isEar[i] = (genDcel.halfedges[genDcel.vertices[i].halfedge].twin==-1)&&(genDcel.halfedges[genDcel.halfedges[genDcel.vertices[i].halfedge].prev].twin==-1);
            if (isEar[i]) isPureTriangle[i]=false;
        }
        
        //realigning halfedges in hex vertices to only follow other hex edges
        for (int i=0;i<genDcel.vertices.size();i++){
            if ((isPureTriangle[i])||(!genDcel.vertices[i].valid))
                continue;
            
            vector<int> hexHEorder;
            int hebegin = genDcel.vertices[i].halfedge;
            if (isBoundary[i]){
                //finding the first hex halfedge
                while (genDcel.halfedges[genDcel.halfedges[hebegin].prev].twin!=-1)
                    hebegin =genDcel.halfedges[genDcel.halfedges[hebegin].prev].twin;
            }
            
            int heiterate=hebegin;
            do{
                if ((genDcel.halfedges[heiterate].data.isFunction)||(genDcel.halfedges[heiterate].twin==-1))
                    hexHEorder.push_back(heiterate);
                if (genDcel.halfedges[heiterate].twin==-1)
                    break;
                heiterate = genDcel.halfedges[genDcel.halfedges[heiterate].twin].next;
            }while(heiterate!=hebegin);
            
            
            for (int j=0;j<hexHEorder.size();j++){
                if ((isBoundary[i])&&(j==hexHEorder.size()-1))
                    continue;
                genDcel.halfedges[hexHEorder[(j+1)%hexHEorder.size()]].prev =genDcel.halfedges[hexHEorder[j]].twin;
                genDcel.halfedges[genDcel.halfedges[hexHEorder[j]].twin].next =hexHEorder[(j+1)%hexHEorder.size()];
                genDcel.vertices[genDcel.halfedges[hexHEorder[j]].vertex].halfedge=hexHEorder[j];
            }
            
            if (isBoundary[i]){ //connect first to the prev
                genDcel.halfedges[hexHEorder[0]].prev = genDcel.halfedges[hebegin].prev;
                genDcel.halfedges[genDcel.halfedges[hebegin].prev].next =hexHEorder[0];
                genDcel.vertices[genDcel.halfedges[hexHEorder[0]].vertex].halfedge=hexHEorder[0];
            }
        }
        
        //invalidating all triangle vertices and edges
        for (int i=0;i<genDcel.vertices.size();i++)
            if (isPureTriangle[i])
                genDcel.vertices[i].valid=false;
        
        for (int i=0;i<genDcel.halfedges.size();i++)
            if ((!genDcel.halfedges[i].data.isFunction)&&(genDcel.halfedges[i].twin!=-1))
                genDcel.halfedges[i].valid=genDcel.edges[genDcel.halfedges[i].edge].valid = false;
        
        //realigning faces
        VectorXi visitedHE=VectorXi::Zero(genDcel.halfedges.size());
        VectorXi usedFace=VectorXi::Zero(genDcel.faces.size());
        for (int i=0;i<genDcel.halfedges.size();i++){
            if ((!genDcel.halfedges[i].valid)||(visitedHE[i]!=0))
                continue;
            
            //following the loop and reassigning face
            int currFace=genDcel.halfedges[i].face;
            genDcel.faces[currFace].halfedge=i;
            usedFace[currFace]=1;
            int hebegin=i;
            int heiterate=hebegin;
            int infinityCounter=0;
            do{
                infinityCounter++;
                if (infinityCounter>genDcel.halfedges.size()){
                    std::cout<<"Infinity loop in realigning faces on halfedge "<<i<<std::endl;
                    return false;
                }
                genDcel.halfedges[heiterate].face=currFace;
                heiterate=genDcel.halfedges[heiterate].next;
            }while (heiterate!=hebegin);
        }
        
        int countThree=0;
        for (int i=0;i<genDcel.faces.size();i++)
            if (!usedFace[i])
                genDcel.faces[i].valid=false;
        
        
        //killing perfect ear faces (not doing corners atm)
        //counting valences
        vector<int> Valences(genDcel.vertices.size());
        for (int i=0;i<genDcel.vertices.size();i++)
            Valences[i]=0;
        
        for (int i=0;i<genDcel.halfedges.size();i++){
            if (genDcel.halfedges[i].valid){
                Valences[genDcel.halfedges[i].vertex]++;
                //Valences[Halfedges[Halfedges[i].next].vertex]++;
                if (genDcel.halfedges[i].twin<0)  //should account for the target as well
                    Valences[genDcel.halfedges[genDcel.halfedges[i].next].vertex]++;
            }
        }
        
        for (int i=0;i<genDcel.faces.size();i++){
            if (!genDcel.faces[i].valid)
                continue;
            countThree=0;
            int hebegin = genDcel.faces[i].halfedge;
            int heiterate=hebegin;
            do{
                if (Valences[genDcel.halfedges[heiterate].vertex]>2)
                    countThree++;
                heiterate=genDcel.halfedges[heiterate].next;
            }while (heiterate!=hebegin);
            if (countThree<3){
                do{
                    //std::cout<<"Invalidating Vertex "<<genDcel.halfedges[heiterate].vertex<<" and  halfedge "<<heiterate<<" of valence "<<Valences[genDcel.halfedges[heiterate].vertex]<<std::endl;
                    
                    genDcel.halfedges[heiterate].valid=false;
                    
                    //invalidating edge or assigning it to the twin
                    if (genDcel.halfedges[heiterate].twin!=-1){
                        if (!genDcel.halfedges[genDcel.halfedges[heiterate].twin].valid)
                            genDcel.edges[genDcel.halfedges[heiterate].edge].valid=false;
                        else
                            genDcel.edges[genDcel.halfedges[heiterate].edge].halfedge = genDcel.halfedges[heiterate].twin;
                    } else genDcel.edges[genDcel.halfedges[heiterate].edge].valid=false;
                    
                    if (genDcel.halfedges[heiterate].twin!=-1)
                        genDcel.halfedges[genDcel.halfedges[heiterate].twin].twin=-1;
                    if ((genDcel.halfedges[heiterate].twin==-1)&&(genDcel.halfedges[genDcel.halfedges[heiterate].prev].twin==-1))  //origin is a boundary vertex
                        genDcel.vertices[genDcel.halfedges[heiterate].vertex].valid=false;
                    
                    heiterate=genDcel.halfedges[heiterate].next;
                    
                    
                }while (heiterate!=hebegin);
                genDcel.faces[i].valid=false;
                //return false;
            }
        }
        
        //need to realign all vertices pointing
        for (int i=0;i<genDcel.halfedges.size();i++)
            if (genDcel.halfedges[i].valid)
                genDcel.vertices[genDcel.halfedges[i].vertex].halfedge=i;
        
        
        if (!genDcel.check_consistency(mData.verbose, true, true, true))
            return false;
        
        for (int i=0;i<Valences.size();i++)
            if ((genDcel.vertices[i].valid)&&(Valences[i]<2))
                genDcel.vertices[i].valid=false;
        
        for (int i=0;i<genDcel.vertices.size();i++){
            if ((genDcel.vertices[i].valid)&&(Valences[i]<=2)&&(!isEar[i])) {
                genDcel.unify_edges(genDcel.vertices[i].halfedge);
                /*if (!genDcel.check_consistency(verbose, true, true, true))
                 return false;*/
            }
        }
        
        if (!genDcel.check_consistency(mData.verbose, true, true, true))
            return false;
        
        //remove non-valid components
        genDcel.clean_mesh();
        
        //checking if mesh is valid
        if (!genDcel.check_consistency(mData.verbose, true, true, true))
            return false;
        
        return true;
        
    }
    
    void RemoveDegree2Faces();
    
    /*void Allocate(int NumofVertices, int NumofFaces, int NumofHEdges)
     {
     Vertices.resize(NumofVertices);
     faces.resize(NumofFaces);
     Halfedges.resize(NumofHEdges);
     }*/
    
    /*void init(const TriMesh& origMesh,
     const Eigen::MatrixXd& cutV,
     const Eigen::MatrixXi& cutF,
     const Eigen::VectorXd& vertexNFunction,
     const int N,
     const Eigen::SparseMatrix<double>& vertexToCornerMat,
     const Eigen::SparseMatrix<int>& exactVertexToCornerMat,
     const Eigen::VectorXi& integerVars,
     const unsigned long resolution=1e7)*/
    void init(const unsigned long resolution=1e7){
        
        using namespace std;
        using namespace Eigen;
        
        //computing exact rational corner values by quantizing the free variables d and then manually performing the sparse matrix multiplication
        vector<ENumber> exactVertexNFunction(mData.vertexNFunction.size());
        double tol = 1.0/(double)resolution;
        for (int i=0;i<mData.vertexNFunction.size();i++){
            //exactVertexNFunction[i]=ENumber((long long)round((long double)(mData.vertexNFunction(i)*resolution)),(long long)resolution);
            exactVertexNFunction[i]=ENumber(mData.vertexNFunction(i),tol);
            
            /*if (abs(exactVertexNFunction[i].to_double() - mData.vertexNFunction(i))>2.0/(double)resolution) {
             cout << "exactVertexNFunction[i].to_double(): " << exactVertexNFunction[i].to_double() << endl;
             cout << "vertexNFunction(i): " << mData.vertexNFunction(i) << endl;
             cout << "(long double)(vertexNFunction(i)*resolution): " << (long double)(mData.vertexNFunction(i) * resolution) << endl;
             }*/
        }
        
        for (int i=0;i<mData.integerVars.size();i++){
            exactVertexNFunction[mData.integerVars(i)]=ENumber((long)round(mData.vertexNFunction(mData.integerVars(i))));
            //cout<<"rounding diff of integer var "<<mData.integerVars(i)<<" is "<<exactVertexNFunction[mData.integerVars(i)].to_double()-mData.vertexNFunction(mData.integerVars(i))<<endl;
        }
        
        VectorXd cutNFunctionVec = mData.orig2CutMat*mData.vertexNFunction;
        vector<ENumber> exactCutNFunctionVec;
        exactSparseMult(mData.exactOrig2CutMat, exactVertexNFunction,exactCutNFunctionVec);
        
        //sanity check - comparing exact to double
        double maxError2 = -32767000.0;
        for (int i=0;i<exactCutNFunctionVec.size();i++){
            double fromExact = exactCutNFunctionVec[i].to_double();
            if (abs(fromExact-cutNFunctionVec[i])>maxError2){
                maxError2 =abs(fromExact-cutNFunctionVec[i]);
                //cout<<"i, fromExact, cutNFunctionVec[i]: "<<i<<","<<fromExact<<","<<cutNFunctionVec[i]<<endl;
            }
        }
        
        if (mData.verbose)
            cout<<"double from exact in halfedges maxError2: "<<maxError2<<endl;
        
        exactNFunction.resize(origMesh.F.size());
        NFunction.resize(origMesh.F.size(), 3*mData.N);
        
        for (int i=0;i<origMesh.F.rows();i++){
            exactNFunction[i].resize(3*mData.N);
            for (int j=0;j<3;j++){
                //Halfedges[FH(i,j)].exactNFunction.resize(N);
                NFunction.block(i, mData.N*j, 1, mData.N) = cutNFunctionVec.segment(mData.N*mData.cutF(i,j), mData.N).transpose();
                for (int k=0;k<mData.N;k++)
                    exactNFunction[i][j*mData.N+k] = exactCutNFunctionVec[mData.N*mData.cutF(i,j)+k];
            }
        }
        
    }
    
    
    //corner angles is per vertex in each F
    void to_polygonal(Eigen::MatrixXd& generatedV,
                      Eigen::VectorXi& generatedD,
                      Eigen::MatrixXi& generatedF){
        generatedV.resize(genDcel.vertices.size(),3);
        
        generatedD.resize(genDcel.faces.size());
        
        for (int i=0;i<genDcel.vertices.size();i++)
            generatedV.row(i)=genDcel.vertices[i].data.coords;
        
        
        for (int i=0;i<genDcel.faces.size();i++){
            int hebegin = genDcel.faces[i].halfedge;
            //reseting to first vertex
            int vCount=0;
            int heiterate=hebegin;
            do{
                vCount++;
                heiterate=genDcel.halfedges[heiterate].next;
            }while (heiterate!=hebegin);
            generatedD(i)=vCount;
        }
        
        generatedF.resize(genDcel.faces.size(),generatedD.maxCoeff());
        for (int i=0;i<genDcel.faces.size();i++){
            int hebegin = genDcel.faces[i].halfedge;
            int vCount=0;
            int heiterate=hebegin;
            do{
                generatedF(i,vCount++)=genDcel.halfedges[heiterate].vertex;
                heiterate=genDcel.halfedges[heiterate].next;
            }while (heiterate!=hebegin);
            
        }
        
        
    }
    
    NFunctionMesher(const TriMesh& _origMesh, const MesherData& _mData ):origMesh(_origMesh), mData(_mData){}
    ~NFunctionMesher(){}
    
private:
    void arrange_on_triangle(const std::vector<EVector2>& triangle,
                             const std::vector<std::pair<EVector2, EVector2>>& lines,
                             const Eigen::VectorXi& lineData,
                             std::vector<EVector2>& V,
                             FunctionDCEL & dcel,
                             Eigen::VectorXi& dataH);
    
    void segment_arrangement(const std::vector<std::pair<EVector2, EVector2>>& segments,
                             const std::vector<int>& data,
                             std::vector<EVector2>& V,
                             FunctionDCEL& dcel,
                             Eigen::VectorXi& dataH);
    
};

}

#endif //DIRECTIONAL_N_FUNCTION_MESHER
