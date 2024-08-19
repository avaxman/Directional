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
#include <gmp.h>
#include <gmpxx.h>
#include <directional/GMP_definitions.h>
#include <directional/dcel.h>
#include <directional/setup_mesh_function_isolines.h>

namespace directional{

    class NFunctionMesher {
    public:

        const TriMesh& origMesh;
        const MeshFunctionIsolinesData& mfiData;

        struct HEData{
            bool isFunction;
            int origHalfedge;
            int origNFunctionIndex;  //the original parameteric function assoicated with this edge
            double prescribedAngle;  //the actual prescribed angle

            HEData():isFunction(false), origHalfedge(-1), origNFunctionIndex(-1),  prescribedAngle(-1.0){}
        };

        struct VData{
            Eigen::RowVector3d coords;
            EVector3 eCoords;
        };

        typedef DCEL<VData,  HEData, bool, bool> FunctionDCEL;
        FunctionDCEL genDcel;

        //vertex quantities
        //Eigen::MatrixXd coordinates;
        //std::vector<EVector3> ECoordinates;

        //halfedge quantities
        Eigen::MatrixXd NFunction;
        std::vector<std::vector<ENumber>> exactNFunction;
        /*std::vector<bool> isHalfedgeFunction;
        std::vector<int> origHalfedge;
        std::vector<int> origNFunctionIndex;  //the original parameteric function assoicated with this edge
        //int prescribedAngleDiff;
        std::vector<double> prescribedAngle;  //the actual prescribed angle

        //face quantities
        std::vector<int> origFace;  //in triangle mesh*/

        /*struct EdgeData{
            int ID;
            bool isFunction;
            int OrigHalfedge;
            bool isBoundary;
            int funcNum;  //in case of function segment

            EdgeData():ID(-1), isFunction(false), OrigHalfedge(-1), isBoundary(false), funcNum(-1){}
            ~EdgeData(){};
        };*/

        /*struct VertexCompare{
            const bool operator<(const std::pair<EVector2, int>& v1, const std::pair<EVector2,int> & v2) const {return (v1.first < v2.first);}
        };*/

        //mesh generation functions found in generate_mesh.h
        void arrange_on_triangle(const std::vector<EVector2>& triangle,
                                 const std::vector<Line2>& lines,
                                 const std::vector<int>& lineData,
                                 std::vector<EVector2>& V,
                                 FunctionDCEL& triDcel);

        void segment_arrangement(const std::vector<Segment2>& segments,
                                 const std::vector<int>& data,
                                 std::vector<EVector2>& V,
                                 FunctionDCEL& triDcel);

        void generate_mesh(const unsigned long Resolution);




        /*class Vertex{
        public:
            int ID;
            Eigen::RowVector3d Coordinates;
            EVector3 ECoordinates;
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
            Eigen::VectorXd NFunction;
            std::vector<ENumber> exactNFunction;
            bool isFunction;
            bool Valid;

            //Parametric function values
            int OrigHalfedge;
            int OrigNFunctionIndex;  //the original parameteric function assoicated with this edge
            //int prescribedAngleDiff;
            double prescribedAngle;  //the actual prescribed angle


            Halfedge():ID(-1), Origin(-1), Next(-1), Prev(-1), Twin(-1), AdjFace(-1), isFunction(false), Valid(true), OrigHalfedge(-1), OrigNFunctionIndex(-1),  prescribedAngle(-1.0){}
            ~Halfedge(){}
        };


        class Face{
        public:
            int ID;
            int AdjHalfedge;
            Eigen::RowVector3d Normal;
            Eigen::RowVector3d Centroid;
            Eigen::RowVector3d Basis1, Basis2;
            bool Valid;

            Face():ID(-1), AdjHalfedge(-1), Valid(true){}
            ~Face(){}
        };

        std::vector<Vertex> Vertices;
        std::vector<Halfedge> Halfedges;
        std::vector<Face> Faces;*/

        std::vector<int> TransVertices;
        std::vector<int> InStrip;
        std::vector<std::set<int> > VertexChains;


        struct MergeData {
            const bool operator()(const int &v1, const int &v2) const { return v1; }
        };

        /*bool JoinFace(const int heindex) {
            if (Halfedges[heindex].twin < 0)
                return true;  //there is no joining of boundary faces

            int Face1 = Halfedges[heindex].AdjFace;
            int Face2 = Halfedges[Halfedges[heindex].twin].AdjFace;

            int hebegin = faces[Face1].halfedge;
            int heiterate = hebegin;
            do {
                heiterate = Halfedges[heiterate].next;
            } while (heiterate != hebegin);

            hebegin = faces[Face1].halfedge;
            heiterate = hebegin;
            do {
                heiterate = Halfedges[heiterate].next;
            } while (heiterate != hebegin);

            //check if spike edge
            if ((Halfedges[heindex].prev == Halfedges[heindex].twin) ||
                (Halfedges[heindex].next == Halfedges[heindex].twin)) {


                int CloseEdge = heindex;
                if (Halfedges[heindex].prev == Halfedges[heindex].twin)
                    CloseEdge = Halfedges[heindex].twin;

                Halfedges[CloseEdge].valid = Halfedges[Halfedges[CloseEdge].twin].valid = false;

                Vertices[Halfedges[CloseEdge].vertex].halfedge = Halfedges[Halfedges[CloseEdge].twin].next;
                faces[Face1].halfedge = Halfedges[CloseEdge].prev;

                Halfedges[Halfedges[CloseEdge].prev].next = Halfedges[Halfedges[CloseEdge].twin].next;
                Halfedges[Halfedges[Halfedges[CloseEdge].twin].next].prev = Halfedges[CloseEdge].prev;

                Vertices[Halfedges[Halfedges[CloseEdge].twin].vertex].valid = false;
                //faces[Face1].NumVertices-=2;  //although one vertex should appear twice


                int hebegin = faces[Face1].halfedge;
                int heiterate = hebegin;
                do {

                    heiterate = Halfedges[heiterate].next;
                } while (heiterate != hebegin);

                hebegin = faces[Face1].halfedge;
                heiterate = hebegin;
                do {
                    heiterate = Halfedges[heiterate].next;
                } while (heiterate != hebegin);


                return true;
            }

            if (Face1 == Face2)
                return false;  //we don't remove non-spike edges with the same faces to not disconnect a chain

            hebegin = faces[Face2].halfedge;
            heiterate = hebegin;
            do {
                heiterate = Halfedges[heiterate].next;
            } while (heiterate != hebegin);

            faces[Face1].halfedge = Halfedges[heindex].next;
            faces[Face2].valid = false;

            //faces[Face2].halfedge=Halfedges[Halfedges[heindex].twin].next;

            Halfedges[heindex].valid = Halfedges[Halfedges[heindex].twin].valid = false;

            Halfedges[Halfedges[heindex].next].prev = Halfedges[Halfedges[heindex].twin].prev;
            Halfedges[Halfedges[Halfedges[heindex].twin].prev].next = Halfedges[heindex].next;

            Halfedges[Halfedges[Halfedges[heindex].twin].next].prev = Halfedges[heindex].prev;
            Halfedges[Halfedges[heindex].prev].next = Halfedges[Halfedges[heindex].twin].next;

            Vertices[Halfedges[heindex].vertex].halfedge = Halfedges[Halfedges[heindex].twin].next;
            Vertices[Halfedges[Halfedges[heindex].next].vertex].halfedge = Halfedges[heindex].next;

            //all other floating halfedges should renounce this one
            for (int i = 0; i < Halfedges.size(); i++)
                if (Halfedges[i].AdjFace == Face2)
                    Halfedges[i].AdjFace = Face1;

            //faces[Face1].NumVertices+=faces[Face2].NumVertices-2;

            //DebugLog<<"Official number of vertices: "<<faces[Face1].NumVertices;

            hebegin = faces[Face1].halfedge;
            heiterate = hebegin;
            int currVertex = 0;
            do {
                //faces[Face1].Vertices[currVertex++]=Halfedges[heiterate].vertex;
                heiterate = Halfedges[heiterate].next;
            } while (heiterate != hebegin);

            return true;
        }*/

        /*void UnifyEdges(int heindex) {
            //if (Halfedges[heindex].twin<0)
            //  return;
            //adjusting source
            Vertices[Halfedges[heindex].vertex].valid = false;
            Halfedges[heindex].vertex = Halfedges[Halfedges[heindex].prev].vertex;
            if (Halfedges[heindex].prescribedAngle < 0.0)
                Halfedges[heindex].prescribedAngle = Halfedges[Halfedges[heindex].prev].prescribedAngle;
            Vertices[Halfedges[heindex].vertex].halfedge = heindex;

            faces[Halfedges[heindex].AdjFace].halfedge = Halfedges[heindex].next;
            //faces[Halfedges[heindex].AdjFace].NumVertices--;



            //adjusting halfedges
            Halfedges[Halfedges[heindex].prev].valid = false;
            Halfedges[heindex].prev = Halfedges[Halfedges[heindex].prev].prev;
            Halfedges[Halfedges[heindex].prev].next = heindex;

            //adjusting twin, if exists
            if (Halfedges[heindex].twin >= 0) {
                if (Halfedges[Halfedges[heindex].twin].prescribedAngle < 0.0)
                    Halfedges[Halfedges[heindex].twin].prescribedAngle = Halfedges[Halfedges[Halfedges[heindex].twin].next].prescribedAngle;
                Halfedges[Halfedges[Halfedges[heindex].twin].next].valid = false;
                Halfedges[Halfedges[heindex].twin].next = Halfedges[Halfedges[Halfedges[heindex].twin].next].next;
                Halfedges[Halfedges[Halfedges[heindex].twin].next].prev = Halfedges[heindex].twin;
                faces[Halfedges[Halfedges[heindex].twin].AdjFace].halfedge = Halfedges[Halfedges[heindex].twin].next;
                //faces[Halfedges[Halfedges[heindex].twin].AdjFace].NumVertices--;
            }
        }

        bool CheckMesh(const bool verbose, const bool checkHalfedgeRepetition, const bool CheckTwinGaps,
                       const bool checkPureBoundary) {
            for (int i = 0; i < Vertices.size(); i++) {
                if (!Vertices[i].valid)
                    continue;

                if (Vertices[i].halfedge == -1) {
                    if (verbose) std::cout << "Valid Vertex " << i << " points to non-valid value -1 " << std::endl;
                    return false;
                }

                if (!Halfedges[Vertices[i].halfedge].valid) {
                    if (verbose)
                        std::cout << "Valid Vertex " << i << " points to non-valid halfedge "
                                  << Vertices[i].halfedge << std::endl;
                    return false;
                }


                if (Halfedges[Vertices[i].halfedge].vertex != i) {
                    if (verbose)
                        std::cout << "Adjacent Halfedge " << Vertices[i].halfedge << " of vertex " << i
                                  << "does not point back" << std::endl;
                    return false;
                }

            }

            for (int i = 0; i < Halfedges.size(); i++) {
                if (!Halfedges[i].valid)
                    continue;


                if (Halfedges[i].next == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Next non-valid value -1" << std::endl;
                    return false;
                }

                if (Halfedges[i].prev == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Prev non-valid value -1" << std::endl;
                    return false;
                }


                if (Halfedges[i].vertex == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Origin non-valid value -1" << std::endl;
                    return false;
                }

                if (Halfedges[i].AdjFace == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to AdjFace non-valid value -1" << std::endl;
                    return false;
                }

                if (Halfedges[Halfedges[i].next].prev != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "Next is " << Halfedges[i].next
                                  << " which doesn't point back as Prev" << std::endl;
                    return false;
                }


                if (Halfedges[Halfedges[i].prev].next != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "Prev is " << Halfedges[i].prev
                                  << " which doesn't point back as Next" << std::endl;
                    return false;
                }

                if (!Vertices[Halfedges[i].vertex].valid) {
                    if (verbose)
                        std::cout << "The Origin of halfedges " << i << ", vertex " << Halfedges[i].vertex
                                  << " is not valid" << std::endl;
                    return false;
                }

                if (!faces[Halfedges[i].AdjFace].valid) {
                    if (verbose)
                        std::cout << "The face of halfedges " << i << ", face " << Halfedges[i].AdjFace
                                  << " is not valid" << std::endl;
                    return false;
                }

                if (Halfedges[Halfedges[i].next].vertex == Halfedges[i].vertex) {  //a degenerate edge{
                    if (verbose)
                        std::cout << "Halfedge " << i << " with twin" << Halfedges[i].twin
                                  << " is degenerate with vertex " << Halfedges[i].vertex << std::endl;
                    return false;
                }

                if (Halfedges[i].twin >= 0) {
                    if (Halfedges[Halfedges[i].twin].twin != i) {
                        if (verbose)
                            std::cout << "Halfedge " << i << "twin is " << Halfedges[i].twin
                                      << " which doesn't point back" << std::endl;
                        return false;
                    }

                    if (!Halfedges[Halfedges[i].twin].valid) {
                        if (verbose)
                            std::cout << "halfedge " << i << " is twin with invalid halfedge" << Halfedges[i].twin
                                      << std::endl;
                        return false;
                    }
                }

                if (!Halfedges[Halfedges[i].next].valid) {
                    if (verbose)
                        std::cout << "halfedge " << i << " has next invalid halfedge" << Halfedges[i].next << std::endl;
                    return false;
                }

                if (!Halfedges[Halfedges[i].prev].valid) {
                    if (verbose)
                        std::cout << "halfedge " << i << " has prev invalid halfedge" << Halfedges[i].prev << std::endl;
                    return false;
                }

                if (Halfedges[i].isFunction) {  //checking that it is not left alone
                    if (Halfedges[i].prev == Halfedges[i].twin) {
                        if (verbose)
                            std::cout << "Hex halfedge " << i << " has Halfedge " << Halfedges[i].prev
                                      << " and both prev and twin" << std::endl;
                        return false;
                    }


                    if (Halfedges[i].next == Halfedges[i].twin) {
                        if (verbose)
                            std::cout << "Hex halfedge " << i << " has Halfedge " << Halfedges[i].next
                                      << " and both next and twin" << std::endl;
                        return false;
                    }
                }
            }

            std::vector <std::set<int>> halfedgesinFace(faces.size());
            std::vector <std::set<int>> verticesinFace(faces.size());
            for (int i = 0; i < faces.size(); i++) {
                if (!faces[i].valid)
                    continue;

                //if (faces[i].NumVertices<3)  //we never allow this
                //  return false;
                int hebegin = faces[i].halfedge;
                int heiterate = hebegin;
                int NumEdges = 0;
                int actualNumVertices = 0;

                do {
                    if (verticesinFace[i].find(Halfedges[heiterate].vertex) != verticesinFace[i].end())
                        if (verbose)
                            std::cout << "Warning: Vertex " << Halfedges[heiterate].vertex
                                      << " appears more than once in face " << i << std::endl;

                    verticesinFace[i].insert(Halfedges[heiterate].vertex);
                    halfedgesinFace[i].insert(heiterate);
                    actualNumVertices++;
                    if (!Halfedges[heiterate].valid)
                        return false;

                    if (Halfedges[heiterate].AdjFace != i) {
                        if (verbose)
                            std::cout << "Face " << i << " has halfedge " << heiterate << " that does not point back"
                                      << std::endl;
                        return false;
                    }

                    heiterate = Halfedges[heiterate].next;
                    NumEdges++;
                    if (NumEdges > Halfedges.size()) {
                        if (verbose) std::cout << "Infinity loop!" << std::endl;
                        return false;
                    }


                } while (heiterate != hebegin);

                /*if (actualNumVertices!=faces[i].NumVertices){
                  DebugLog<<"faces "<<i<<" lists "<<faces[i].NumVertices<<" vertices but has a chain of "<<actualNumVertices<<endl;
                  return false;
                }

                for (int j=0;j<faces[i].NumVertices;j++){
                  if (faces[i].Vertices[j]<0){
                    DebugLog<<"faces "<<i<<".vertices "<<j<<"is undefined"<<endl;
                    return false;
                  }
                  if (!Vertices[faces[i].Vertices[j]].valid){
                    DebugLog<<"faces "<<i<<".vertices "<<j<<"is not valid"<<endl;
                    return false;
                  }
                }*/
            /*}/*

            //checking if all halfedges that relate to a face are part of its recognized chain (so no floaters)
            for (int i = 0; i < Halfedges.size(); i++) {
                if (!Halfedges[i].valid)
                    continue;
                int currFace = Halfedges[i].AdjFace;
                if (halfedgesinFace[currFace].find(i) == halfedgesinFace[currFace].end()) {
                    if (verbose) std::cout << "Halfedge " << i << " is floating in face " << currFace << std::endl;
                    return false;
                }
            }

            //check if mesh is a manifold: every halfedge appears only once
            if (checkHalfedgeRepetition) {
                std::set <TwinFinder> HESet;
                for (int i = 0; i < Halfedges.size(); i++) {
                    if (!Halfedges[i].valid)
                        continue;
                    std::set<TwinFinder>::iterator HESetIterator = HESet.find(
                            TwinFinder(i, Halfedges[i].vertex, Halfedges[Halfedges[i].next].vertex));
                    if (HESetIterator != HESet.end()) {
                        if (verbose)
                            std::cout << "Warning: the halfedge (" << Halfedges[i].vertex << ","
                                      << Halfedges[Halfedges[i].next].vertex << ") appears at least twice in the mesh"
                                      << std::endl;
                        if (verbose)
                            std::cout << "for instance halfedges " << i << " and " << HESetIterator->index << std::endl;
                        return false;
                        //return false;
                    } else {
                        HESet.insert(TwinFinder(i, Halfedges[i].vertex, Halfedges[Halfedges[i].next].vertex));
                        //if (verbose) std::cout<<"inserting halfedge "<<i<<" which is "<<Halfedges[i].vertex<<", "<<Halfedges[Halfedges[i].next].vertex<<endl;
                    }
                }
            }

            if (CheckTwinGaps) {
                std::set <TwinFinder> HESet;
                //checking if there is a gap: two halfedges that share the same opposite vertices but do not have twins
                for (int i = 0; i < Halfedges.size(); i++) {
                    if (!Halfedges[i].valid)
                        continue;

                    std::set<TwinFinder>::iterator HESetIterator = HESet.find(
                            TwinFinder(i, Halfedges[i].vertex, Halfedges[Halfedges[i].next].vertex));
                    if (HESetIterator == HESet.end()) {
                        HESet.insert(TwinFinder(i, Halfedges[i].vertex, Halfedges[Halfedges[i].next].vertex));
                        continue;
                    }

                    HESetIterator = HESet.find(TwinFinder(i, Halfedges[Halfedges[i].next].vertex, Halfedges[i].vertex));
                    if (HESetIterator != HESet.end()) {

                        if (Halfedges[i].twin == -1) {
                            if (verbose)
                                std::cout << "Halfedge " << i << "has no twin although halfedge "
                                          << HESetIterator->index << " can be a twin" << std::endl;
                            return false;
                        }
                        if (Halfedges[HESetIterator->index].twin == -1) {
                            if (verbose)
                                std::cout << "Halfedge " << HESetIterator->index << "has no twin although halfedge "
                                          << i << " can be a twin" << std::endl;
                            return false;
                        }
                    }
                }
            }

            //checking if there are pure boundary faces (there shouldn't be)
            if (checkPureBoundary) {
                for (int i = 0; i < Halfedges.size(); i++) {
                    if (!Halfedges[i].valid)
                        continue;
                    if ((Halfedges[i].twin < 0) && (Halfedges[i].isFunction))
                        if (verbose)
                            std::cout << "WARNING: Halfedge " << i << " is a hex edge without twin!" << std::endl;

                    if (Halfedges[i].twin > 0)
                        continue;

                    bool pureBoundary = true;
                    int hebegin = i;
                    int heiterate = hebegin;
                    do {
                        if (Halfedges[heiterate].twin > 0) {
                            pureBoundary = false;
                            break;
                        }
                        heiterate = Halfedges[heiterate].next;
                    } while (heiterate != hebegin);
                    if (pureBoundary) {
                        if (verbose)
                            std::cout << "Face " << Halfedges[i].AdjFace << "is a pure boundary face!" << std::endl;
                        return false;
                    }
                }

                //checking for latent valence 2 faces
                std::vector<int> Valences(Vertices.size());
                for (int i = 0; i < Vertices.size(); i++)
                    Valences[i] = 0;

                for (int i = 0; i < Halfedges.size(); i++) {
                    if (Halfedges[i].valid) {
                        Valences[Halfedges[i].vertex]++;
                        //Valences[Halfedges[Halfedges[i].next].vertex]++;
                        if (Halfedges[i].twin < 0)  //should account for the target as well
                            Valences[Halfedges[Halfedges[i].next].vertex]++;
                    }
                }

                int countThree;
                for (int i = 0; i < faces.size(); i++) {
                    if (!faces[i].valid)
                        continue;
                    countThree = 0;
                    int hebegin = faces[i].halfedge;
                    int heiterate = hebegin;
                    int numEdges = 0;
                    do {
                        if (Valences[Halfedges[heiterate].vertex] > 2)
                            countThree++;
                        heiterate = Halfedges[heiterate].next;
                        numEdges++;
                        if (numEdges > Halfedges.size()) {
                            if (verbose) std::cout << "Infinity loop in face " << i << "!" << std::endl;
                            return false;
                        }
                    } while (heiterate != hebegin);
                    if (countThree < 3) {
                        if (verbose) std::cout << "Face " << i << " is a latent valence 2 face!" << std::endl;
                        if (verbose) std::cout << "Its vertices are " << std::endl;
                        do {
                            if (verbose)
                                std::cout << "Vertex " << Halfedges[heiterate].vertex << " halfedge " << heiterate
                                          << " valence " << Valences[Halfedges[heiterate].vertex] << std::endl;

                            if (Valences[Halfedges[heiterate].vertex] > 2)
                                countThree++;
                            heiterate = Halfedges[heiterate].next;
                            numEdges++;
                            if (numEdges > Halfedges.size()) {
                                if (verbose) std::cout << "Infinity loop in face " << i << "!" << std::endl;
                                return false;
                            }
                        } while (heiterate != hebegin);

                        //return false;
                    }
                }
            }

            if (verbose) std::cout << "Mesh is clear according to given checks" << std::endl;
            return true;  //most likely the mesh is solid

        }*/

        /*void CleanMesh() {
            //removing nonvalid vertices
            std::vector<int> TransVertices(Vertices.size());
            std::vector <Vertex> NewVertices;
            for (int i = 0; i < Vertices.size(); i++) {
                if (!Vertices[i].valid)
                    continue;

                Vertex NewVertex = Vertices[i];
                NewVertex.ID = NewVertices.size();
                NewVertices.push_back(NewVertex);
                TransVertices[i] = NewVertex.ID;
            }


            Vertices = NewVertices;
            for (int i = 0; i < Halfedges.size(); i++)
                Halfedges[i].vertex = TransVertices[Halfedges[i].vertex];



            /*for (int i=0;i<faces.size();i++){
              if (!faces[i].valid)
                continue;
              for (int j=0;j<faces[i].NumVertices;j++){
                DebugLog<<"i is "<<i<<", j is "<<j<<", faces[i].Vertices[j] is "<<faces[i].Vertices[j]<<endl;
                faces[i].Vertices[j]=TransVertices[faces[i].Vertices[j]];
              }
            }*/
/*


            //removing nonvalid faces
            std::vector <Face> NewFaces;
            std::vector<int> TransFaces(faces.size());
            for (int i = 0; i < faces.size(); i++) {
                if (!faces[i].valid)
                    continue;

                Face NewFace = faces[i];
                NewFace.ID = NewFaces.size();
                NewFaces.push_back(NewFace);
                TransFaces[i] = NewFace.ID;
            }
            faces = NewFaces;
            for (int i = 0; i < Halfedges.size(); i++)
                Halfedges[i].AdjFace = TransFaces[Halfedges[i].AdjFace];



            //removing nonvalid halfedges
            std::vector <Halfedge> NewHalfedges;
            std::vector<int> TransHalfedges(Halfedges.size());
            for (int i = 0; i < Halfedges.size(); i++) {
                if (!Halfedges[i].valid)
                    continue;

                Halfedge NewHalfedge = Halfedges[i];
                NewHalfedge.ID = NewHalfedges.size();
                NewHalfedges.push_back(NewHalfedge);
                TransHalfedges[i] = NewHalfedge.ID;
            }


            Halfedges = NewHalfedges;
            for (int i = 0; i < faces.size(); i++)
                faces[i].halfedge = TransHalfedges[faces[i].halfedge];


            for (int i = 0; i < Vertices.size(); i++)
                Vertices[i].halfedge = TransHalfedges[Vertices[i].halfedge];


            for (int i = 0; i < Halfedges.size(); i++) {
                if (Halfedges[i].twin != -1)
                    Halfedges[i].twin = TransHalfedges[Halfedges[i].twin];
                Halfedges[i].next = TransHalfedges[Halfedges[i].next];
                Halfedges[i].prev = TransHalfedges[Halfedges[i].prev];
            }

        }*/


        void TestUnmatchedTwins();

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

        std::vector<std::pair<int,int>> FindVertexMatch(const bool verbose, std::vector<EVector3>& Set1, std::vector<EVector3>& Set2)
        {
            std::set<PointPair> PairSet;
            for (int i=0;i<Set1.size();i++)
                for (int j=0;j<Set2.size();j++)
                    PairSet.insert(PointPair(i,j,squaredDistance(Set1[i],Set2[j])));

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
                    if (squaredDistance(Set1[Result[i].first],Set2[Result[i].second])>0){
                        std::cout<<"("<<Result[i].first<<","<<Result[i].second<<") with dist "<<squaredDistance(Set1[Result[i].first],Set2[Result[i].second])<<std::endl;
                        std::cout<<"Distance is abnormally not zero!"<<std::endl;
                    }
                }
            }


            return Result;

        }

        bool simplify_mesh(const bool verbose, int N){
            //unifying vertices which are similar

            using namespace std;
            using namespace Eigen;

            if (!genDcel.check_consistency(verbose, false, false, false))
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
                    CurrMatches=FindVertexMatch(verbose, PointSet1, PointSet2);

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

            if (verbose)
                std::cout<<"Max matching distance: "<<MaxDist<<endl;

            //vector<int> Transvertices(vertices.size());
            TransVertices.resize(genDcel.vertices.size());
            int NumNewVertices = connectedComponents(VertexMatches, TransVertices);

            if (!genDcel.check_consistency(verbose, false, false, false))
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


            if (!genDcel.check_consistency(verbose, true, false, false))
                return false;

            //twinning up edges
            set<FunctionDCEL::TwinFinder> Twinning;
            for (int i=0;i<genDcel.halfedges.size();i++){
                if ((genDcel.halfedges[i].twin>=0)||(!genDcel.halfedges[i].valid))
                    continue;

                set<FunctionDCEL::TwinFinder>::iterator Twinit=Twinning.find(FunctionDCEL::TwinFinder(0,genDcel.halfedges[genDcel.halfedges[i].next].vertex,
                                                                                                      genDcel.halfedges[i].vertex));
                if (Twinit!=Twinning.end()){
                    if ((genDcel.halfedges[Twinit->index].twin!=-1)&&(verbose))
                        std::cout<<"warning: halfedge "<<Twinit->index<<" is already twinned to halfedge "<<genDcel.halfedges[Twinit->index].twin<<std::endl;
                    if ((genDcel.halfedges[i].twin!=-1)&&(verbose))
                        std::cout<<"warning: halfedge "<<i<<" is already twinned to halfedge "<<genDcel.halfedges[Twinit->index].twin<<std::endl;
                    genDcel.halfedges[Twinit->index].twin=i;
                    genDcel.halfedges[i].twin=Twinit->index;

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


            if (!genDcel.check_consistency(verbose, true, true, true))
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
                    genDcel.halfedges[i].valid=false;

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
                        /*DebugLog<<"Invalidating Vertex "<<Halfedges[heiterate].vertex<<"and  halfedge "<<heiterate<<" of valence "<<Valences[Halfedges[heiterate].vertex]<<endl;*/

                        genDcel.halfedges[heiterate].valid=false;
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


            if (!genDcel.check_consistency(verbose, true, true, true))
                return false;

            for (int i=0;i<Valences.size();i++)
                if ((genDcel.vertices[i].valid)&&(Valences[i]<2))
                    genDcel.vertices[i].valid=false;

            for (int i=0;i<genDcel.vertices.size();i++){
                if ((genDcel.vertices[i].valid)&&(Valences[i]<=2)&&(!isEar[i]))
                    genDcel.unify_edges(genDcel.vertices[i].halfedge);
            }

            if (!genDcel.check_consistency(verbose, true, true, true))
                return false;

            //remove non-valid components
            genDcel.clean_mesh();

            //checking if mesh is valid
            if (!genDcel.check_consistency(verbose, true, true, true))
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
            vector<ENumber> exactVertexNFunction(mfiData.vertexNFunction.size());
            for (int i=0;i<mfiData.vertexNFunction.size();i++){
                exactVertexNFunction[i]=ENumber((signed long)round((long double)(mfiData.vertexNFunction(i)*resolution)),(unsigned long)resolution);
                /*if (abs(exactVertexNFunction[i].to_double() - vertexNFunction(i))>10e-8) {
                    cout << "exactVertexNFunction[i].to_double(): " << exactVertexNFunction[i].to_double() << endl;
                    cout << "vertexNFunction(i): " << vertexNFunction(i) << endl;
                    cout << "(long double)(vertexNFunction(i)*resolution): " << (long double)(vertexNFunction(i) * resolution) << endl;
                }*/
            }

            for (int i=0;i<mfiData.integerVars.size();i++){
                exactVertexNFunction[mfiData.integerVars(i)]=ENumber((long)round(mfiData.vertexNFunction(mfiData.integerVars(i))));
                //cout<<"rounding diff of integer var "<<integerVars(i)<<" is "<<exactVertexNFunction[integerVars(i)].to_double()-vertexNFunction(integerVars(i))<<endl;
            }

            VectorXd cutNFunctionVec = mfiData.orig2CutMat*mfiData.vertexNFunction;
            vector<ENumber> exactCutNFunctionVec;
            exactSparseMult(mfiData.exactOrig2CutMat, exactVertexNFunction,exactCutNFunctionVec);

            //sanity check - comparing exact to double
            double maxError2 = -32767000.0;
            for (int i=0;i<exactCutNFunctionVec.size();i++){
                double fromExact = exactCutNFunctionVec[i].get_d();
                if (abs(fromExact-cutNFunctionVec[i])>maxError2)
                    maxError2 =abs(fromExact-cutNFunctionVec[i]);
            }

            cout<<"double from exact in halfedges maxError2: "<<maxError2<<endl;

            exactNFunction.resize(origMesh.F.size());
            NFunction.resize(origMesh.F.size(), 3*mfiData.N);

            for (int i=0;i<origMesh.F.rows();i++){
                exactNFunction[i].resize(3*mfiData.N);
                for (int j=0;j<3;j++){
                    //Halfedges[FH(i,j)].exactNFunction.resize(N);
                    NFunction.block(i, mfiData.N*j, 1, mfiData.N) = cutNFunctionVec.segment(mfiData.N*mfiData.cutF(i,j), mfiData.N).transpose();
                    for (int k=0;k<mfiData.N;k++)
                        exactNFunction[i][j*mfiData.N+k] = exactCutNFunctionVec[mfiData.N*mfiData.cutF(i,j)+k];
                }
            }

            //sanity check
            /*double maxError = -32767000.0;
            for (int i=0;i<Halfedges.size();i++){
                for (int j=0;j<N;j++){
                    double fromExact = Halfedges[i].exactNFunction[j].to_double();
                    //cout<<"fromExact: "<<fromExact<<endl;
                    //cout<<"Halfedges[i].NFunction[j]: "<<Halfedges[i].NFunction[j]<<endl;
                    if (abs(fromExact-Halfedges[i].NFunction[j])>maxError)
                        maxError =abs(fromExact-Halfedges[i].NFunction[j]);
                }

            }*/
            //cout<<"double from exact in halfedges maxError: "<<maxError<<endl;
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

            /*generatedFfuncNum.resize(faces.size(),generatedD.maxCoeff());
            cornerAngles=Eigen::MatrixXd::Constant(faces.size(),generatedD.maxCoeff(),-1.0);
            //prescribedAnglesInt.resize(faces.size(),generatedD.maxCoeff());
            for (int i=0;i<faces.size();i++){
              int hebegin = faces[i].halfedge;
              int vCount=0;
              int heiterate=hebegin;
              do{
                generatedFfuncNum(i,vCount)=halfedges[heiterate].OrigNFunctionIndex;
                cornerAngles(i,vCount++)=halfedges[heiterate].prescribedAngle;
                //prescribedAnglesInt(i,vCount++)=halfedges[heiterate].prescribedAngleDiff;
                //cout<<"halfedges[heiterate].prescribedAngleDiff: "<<halfedges[heiterate].prescribedAngleDiff<<endl;
                heiterate=halfedges[heiterate].next;
              }while (heiterate!=hebegin);
            }*/


        }

        NFunctionMesher(const TriMesh& _origMesh, const MeshFunctionIsolinesData& _mfiData ):origMesh(_origMesh), mfiData(_mfiData){}
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
