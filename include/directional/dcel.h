// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_DCEL_H
#define DIRECTIONAL_DCEL_H

#include <Eigen/Core>
#include <vector>


namespace directional
{
    template<typename VertexData, typename HalfedgeData, typename EdgeData, typename FaceData>
    class DCEL{
    public:

        struct Vertex{
            int ID;
            bool valid;
            int halfedge;
            VertexData data;

            Vertex():valid(true), halfedge(-1), ID(-1){};
        };

        struct Halfedge{
            int ID;
            bool valid;
            int vertex, face, edge;
            int next, prev, twin;
            HalfedgeData data;

            Halfedge():valid(true), vertex(-1), face(-1), edge(-1), next(-1), prev(-1), twin(-1){}
        };

        struct Edge{
            int ID;
            bool valid;
            int halfedge;
            EdgeData data;

            Edge():valid(true), halfedge(-1){}
        };

        struct Face{
            int ID;
            bool valid;
            int halfedge;
            FaceData data;

            Face():valid(true), halfedge(-1){}
        };

        std::vector<Vertex> vertices;
        std::vector<Halfedge> halfedges;
        std::vector<Edge> edges;
        std::vector<Face> faces;

        DCEL(){}
        ~DCEL(){}

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

        void walk_boundary(int &CurrHalfedge) {
            do {
                CurrHalfedge = halfedges[CurrHalfedge].next;
                if (halfedges[CurrHalfedge].twin < 0)
                    break;  //next boundary over a 2-valence vertex
                CurrHalfedge = halfedges[CurrHalfedge].twin;
            } while (halfedges[CurrHalfedge].twin >= 0);
        }

        void stitch_twins() {
            //twinning up edges
            std::set <TwinFinder> Twinning;
            for (int i = 0; i < halfedges.size(); i++) {
                if (halfedges[i].twin >= 0)
                    continue;

                typename std::set<TwinFinder>::iterator Twinit = Twinning.find(
                        TwinFinder(0, halfedges[halfedges[i].next].vertex, halfedges[i].vertex));
                if (Twinit != Twinning.end()) {
                    halfedges[Twinit->index].twin = i;
                    halfedges[i].twin = Twinit->index;
                    Twinning.erase(*Twinit);
                } else {
                    Twinning.insert(TwinFinder(i, halfedges[i].origin, halfedges[halfedges[i].next].origin));
                }
            }
        }

        bool check_consistency(const bool verbose, const bool checkHalfedgeRepetition, const bool CheckTwinGaps,
                       const bool checkPureBoundary) {
            for (int i = 0; i < vertices.size(); i++) {
                if (!vertices[i].valid)
                    continue;

                if (vertices[i].halfedge == -1) {
                    if (verbose) std::cout << "Valid Vertex " << i << " points to non-valid value -1 " << std::endl;
                    return false;
                }

                if (!halfedges[vertices[i].halfedge].valid) {
                    if (verbose)
                        std::cout << "Valid Vertex " << i << " points to non-valid AdjHalfedge "
                                  << vertices[i].halfedge << std::endl;
                    return false;
                }


                if (halfedges[vertices[i].halfedge].valid != i) {
                    if (verbose)
                        std::cout << "Adjacent Halfedge " << vertices[i].halfedge << " of vertex " << i
                                  << "does not point back" << std::endl;
                    return false;
                }

            }

            for (int i = 0; i < halfedges.size(); i++) {
                if (!halfedges[i].valid)
                    continue;


                if (halfedges[i].next == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Next non-valid value -1" << std::endl;
                    return false;
                }

                if (halfedges[i].prev == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Prev non-valid value -1" << std::endl;
                    return false;
                }


                if (halfedges[i].vertex == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Origin non-valid value -1" << std::endl;
                    return false;
                }

                if (halfedges[i].face == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to AdjFace non-valid value -1" << std::endl;
                    return false;
                }

                if (halfedges[halfedges[i].next].prev != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "Next is " << halfedges[i].next
                                  << " which doesn't point back as Prev" << std::endl;
                    return false;
                }


                if (halfedges[halfedges[i].prev].next != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "Prev is " << halfedges[i].prev
                                  << " which doesn't point back as Next" << std::endl;
                    return false;
                }

                if (!vertices[halfedges[i].vertex].valid) {
                    if (verbose)
                        std::cout << "The Origin of halfedges " << i << ", vertex " << halfedges[i].vertex
                                  << " is not valid" << std::endl;
                    return false;
                }

                if (!faces[halfedges[i].face].valid) {
                    if (verbose)
                        std::cout << "The face of halfedges " << i << ", face " << halfedges[i].face
                                  << " is not valid" << std::endl;
                    return false;
                }

                if (halfedges[halfedges[i].next].vertex == halfedges[i].vertex) {  //a degenerate edge{
                    if (verbose)
                        std::cout << "Halfedge " << i << " with twin" << halfedges[i].twin
                                  << " is degenerate with vertex " << halfedges[i].vertex << std::endl;
                    return false;
                }

                if (halfedges[i].twin >= 0) {
                    if (halfedges[halfedges[i].twin].twin != i) {
                        if (verbose)
                            std::cout << "Halfedge " << i << "twin is " << halfedges[i].twin
                                      << " which doesn't point back" << std::endl;
                        return false;
                    }

                    if (!halfedges[halfedges[i].twin].valid) {
                        if (verbose)
                            std::cout << "halfedge " << i << " is twin with invalid halfedge" << halfedges[i].twin
                                      << std::endl;
                        return false;
                    }
                }

                if (!halfedges[halfedges[i].next].valid) {
                    if (verbose)
                        std::cout << "halfedge " << i << " has next invalid halfedge" << halfedges[i].next << std::endl;
                    return false;
                }

                if (!halfedges[halfedges[i].prev].valid) {
                    if (verbose)
                        std::cout << "halfedge " << i << " has prev invalid halfedge" << halfedges[i].prev << std::endl;
                    return false;
                }

                //if (Halfedges[i].isFunction) {  //checking that it is not left alone
                if (halfedges[i].prev == halfedges[i].twin) {
                    if (verbose)
                        std::cout << "Hex halfedge " << i << " has Halfedge " << halfedges[i].prev
                                  << " and both prev and twin" << std::endl;
                    return false;
                }


                if (halfedges[i].next == halfedges[i].twin) {
                    if (verbose)
                        std::cout << "Hex halfedge " << i << " has Halfedge " << halfedges[i].next
                                  << " and both next and twin" << std::endl;
                    return false;
                }
                //}
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
                    if (verticesinFace[i].find(halfedges[heiterate].vertex) != verticesinFace[i].end())
                        if (verbose)
                            std::cout << "Warning: Vertex " << halfedges[heiterate].vertex
                                      << " appears more than once in face " << i << std::endl;

                    verticesinFace[i].insert(halfedges[heiterate].vertex);
                    halfedgesinFace[i].insert(heiterate);
                    actualNumVertices++;
                    if (!halfedges[heiterate].valid)
                        return false;

                    if (halfedges[heiterate].face != i) {
                        if (verbose)
                            std::cout << "Face " << i << " has halfedge " << heiterate << " that does not point back"
                                      << std::endl;
                        return false;
                    }

                    heiterate = halfedges[heiterate].next;
                    NumEdges++;
                    if (NumEdges > halfedges.size()) {
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
                  if (!Vertices[faces[i].Vertices[j]].Valid){
                    DebugLog<<"faces "<<i<<".vertices "<<j<<"is not valid"<<endl;
                    return false;
                  }
                }*/
            }

            //checking if all halfedges that relate to a face are part of its recognized chain (so no floaters)
            for (int i = 0; i < halfedges.size(); i++) {
                if (!halfedges[i].valid)
                    continue;
                int currFace = halfedges[i].face;
                if (halfedgesinFace[currFace].find(i) == halfedgesinFace[currFace].end()) {
                    if (verbose) std::cout << "Halfedge " << i << " is floating in face " << currFace << std::endl;
                    return false;
                }
            }

            //check if mesh is a manifold: every halfedge appears only once
            if (checkHalfedgeRepetition) {
                std::set <TwinFinder> HESet;
                for (int i = 0; i < halfedges.size(); i++) {
                    if (!halfedges[i].valid.Valid)
                        continue;
                    typename std::set<TwinFinder>::iterator HESetIterator = HESet.find(
                            TwinFinder(i, halfedges[i].vertex, halfedges[halfedges[i].next].vertex));
                    if (HESetIterator != HESet.end()) {
                        if (verbose)
                            std::cout << "Warning: the halfedge (" << halfedges[i].vertex << ","
                                      << halfedges[halfedges[i].next].vertex << ") appears at least twice in the mesh"
                                      << std::endl;
                        if (verbose)
                            std::cout << "for instance halfedges " << i << " and " << HESetIterator->index << std::endl;
                        return false;
                        //return false;
                    } else {
                        HESet.insert(TwinFinder(i, halfedges[i].vertex, halfedges[halfedges[i].next].vertex));
                        //if (verbose) std::cout<<"inserting halfedge "<<i<<" which is "<<Halfedges[i].vertex<<", "<<Halfedges[halfedge[i].next].vertex<<endl;
                    }
                }
            }

            if (CheckTwinGaps) {
                std::set <TwinFinder> HESet;
                //checking if there is a gap: two halfedges that share the same opposite vertices but do not have twins
                for (int i = 0; i < halfedges.size(); i++) {
                    if (!halfedges[i].valid)
                        continue;

                    typename std::set<TwinFinder>::iterator HESetIterator = HESet.find(
                            TwinFinder(i, halfedges[i].vertex, halfedges[halfedges[i].next].vertex));
                    if (HESetIterator == HESet.end()) {
                        HESet.insert(TwinFinder(i, halfedges[i].vertex, halfedges[halfedges[i].next].vertex));
                        continue;
                    }

                    HESetIterator = HESet.find(TwinFinder(i, halfedges[halfedges[i].next].vertex, halfedges[i].vertex));
                    if (HESetIterator != HESet.end()) {

                        if (halfedges[i].twin == -1) {
                            if (verbose)
                                std::cout << "Halfedge " << i << "has no twin although halfedge "
                                          << HESetIterator->index << " can be a twin" << std::endl;
                            return false;
                        }
                        if (halfedges[HESetIterator->index].twin == -1) {
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
                for (int i = 0; i < halfedges.size(); i++) {
                    if (!halfedges[i].valid)
                        continue;

                    //Is this necessary?
                    /*if ((halfedge[i].twin < 0) && (Halfedges[i].isFunction))
                        if (verbose)
                            std::cout << "WARNING: Halfedge " << i << " is a hex edge without twin!" << std::endl;
*/
                    if (halfedges[i].twin > 0)
                        continue;

                    bool pureBoundary = true;
                    int hebegin = i;
                    int heiterate = hebegin;
                    do {
                        if (halfedges[heiterate].twin > 0) {
                            pureBoundary = false;
                            break;
                        }
                        heiterate = halfedges[heiterate].next;
                    } while (heiterate != hebegin);
                    if (pureBoundary) {
                        if (verbose)
                            std::cout << "Face " << halfedges[i].face << "is a pure boundary face!" << std::endl;
                        return false;
                    }
                }

                //checking for latent valence 2 faces
                std::vector<int> Valences(vertices.size());
                for (int i = 0; i < vertices.size(); i++)
                    Valences[i] = 0;

                for (int i = 0; i < halfedges.size(); i++) {
                    if (halfedges[i].valid) {
                        Valences[halfedges[i].vertex]++;
                        //Valences[Halfedges[halfedge[i].next].vertex]++;
                        if (halfedges[i].twin < 0)  //should account for the target as well
                            Valences[halfedges[halfedges[i].next].vertex]++;
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
                        if (Valences[halfedges[heiterate].vertex] > 2)
                            countThree++;
                        heiterate = halfedges[heiterate].next;
                        numEdges++;
                        if (numEdges > halfedges.size()) {
                            if (verbose) std::cout << "Infinity loop in face " << i << "!" << std::endl;
                            return false;
                        }
                    } while (heiterate != hebegin);
                    if (countThree < 3) {
                        if (verbose) std::cout << "Face " << i << " is a latent valence 2 face!" << std::endl;
                        if (verbose) std::cout << "Its vertices are " << std::endl;
                        do {
                            if (verbose)
                                std::cout << "Vertex " << halfedges[heiterate].vertex << " halfedge " << heiterate
                                          << " valence " << Valences[halfedges[heiterate].vertex] << std::endl;

                            if (Valences[halfedges[heiterate].vertex] > 2)
                                countThree++;
                            heiterate = halfedges[heiterate].next;
                            numEdges++;
                            if (numEdges > halfedges.size()) {
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

        }

        //Only used after having checked for consistency!
        void clean_mesh() {

            //Cleaning nonvalid vertices
            std::vector<int> transVertices(vertices.size());
            std::vector<Vertex> newVertices;  //TODO: from here
            //std::vector <Vertex> NewVertices;
            int counter=0;
            for (int i = 0; i < vertices.size(); i++) {
                if (!vertices[i].valid)
                    continue;

                newVertices.push_back(vertices[i]);
                newVertices[newVertices.size()-1].ID = newVertices.size()-1;
                transVertices[i]=newVertices.size()-1;
            }

            vertices=newVertices;
            //updating references to these vertices
            for (int i=0;i<halfedges.size();i++){
                halfedges[i].vertex=transVertices[halfedges[i].vertex];
            }

            //Cleaning nonvalid faces
            std::vector<int> transFaces(faces.size());
            std::vector<Face> newFaces;
            for (int i = 0; i < faces.size(); i++) {
                if (!faces[i].valid)
                    continue;

                newFaces.push_back(faces[i]);
                newFaces[newFaces.size()-1].ID = newFaces.size()-1;
                transFaces[i]=newFaces.size()-1;
            }
            faces = newFaces;
            for (int i = 0; i < halfedges.size(); i++)
                halfedges[i].face = transFaces[halfedges[i].face];

            //Cleaning nonvalid halfedges
            std::vector <Halfedge> newHalfedges;
            std::vector<int> transHalfedges(halfedges.size());
            for (int i = 0; i < halfedges.size(); i++) {
                if (!halfedges[i].Valid)
                    continue;

                Halfedge NewHalfedge = halfedges[i];
                NewHalfedge.ID = newHalfedges.size();
                newHalfedges.push_back(NewHalfedge);
                transHalfedges[i] = NewHalfedge.ID;
            }

            halfedges = newHalfedges;
            for (int i = 0; i < faces.size(); i++)
                faces[i].halfedges = transHalfedges[faces[i].halfedges];

            for (int i = 0; i < vertices.size(); i++)
                vertices[i].halfedge = transHalfedges[vertices[i].halfedge];

            for (int i = 0; i < edges.size(); i++)
                edges[i].halfedge = transHalfedges[edges[i].halfedge];

            for (int i = 0; i < halfedges.size(); i++) {
                if (halfedges[i].twin != -1)
                    halfedges[i].twin = transHalfedges[halfedges[i].twin];
                halfedges[i].next = transHalfedges[halfedges[i].next];
                halfedges[i].prev = transHalfedges[halfedges[i].prev];
            }

            //cleaing non-valid edges
            std::vector<int> transEdges(edges.size());
            std::vector<Vertex> newEdges;
            for (int i = 0; i < edges.size(); i++) {
                if (!edges[i].valid)
                    continue;

                newEdges.push_back(edges[i]);
                newEdges[newEdges.size()-1].ID = newEdges.size()-1;
                transEdges[i]=newEdges.size()-1;
            }

            edges=newEdges;
            //updating references to these vertices
            for (int i=0;i<halfedges.size();i++){
                halfedges[i].edge=transEdges[halfedges[i].edge];
            }
        }

        void ComputeTwins() {
            //twinning up edges
            std::set <TwinFinder> Twinning;
            for (int i = 0; i < halfedges.size(); i++) {
                if (halfedges[i].twin >= 0)
                    continue;

                typename std::set<TwinFinder>::iterator Twinit = Twinning.find(
                        TwinFinder(0, halfedges[halfedges[i].Next].vertex, halfedges[i].vertex));
                if (Twinit != Twinning.end()) {
                    halfedges[Twinit->index].twin = i;
                    halfedges[i].twin = Twinit->index;
                    Twinning.erase(*Twinit);
                } else {
                    Twinning.insert(TwinFinder(i, halfedges[i].vertex, halfedges[halfedges[i].next].vertex));
                }
            }
        }

        bool JoinFace(const int heindex) {
            if (halfedges[heindex].Twin < 0)
                return true;  //there is no joining of boundary faces

            int Face1 = halfedges[heindex].AdjFace;
            int Face2 = halfedges[halfedges[heindex].Twin].AdjFace;

            int hebegin = faces[Face1].halfedge;
            int heiterate = hebegin;
            do {
                heiterate = halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            hebegin = faces[Face1].halfedge;
            heiterate = hebegin;
            do {
                heiterate = halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            //check if spike edge
            if ((halfedges[heindex].Prev == halfedges[heindex].Twin) ||
                (halfedges[heindex].Next == halfedges[heindex].Twin)) {


                int CloseEdge = heindex;
                if (halfedges[heindex].Prev == halfedges[heindex].Twin)
                    CloseEdge = halfedges[heindex].Twin;

                halfedges[CloseEdge].Valid = halfedges[halfedges[CloseEdge].Twin].Valid = false;

                vertices[halfedges[CloseEdge].vertex].halfedge = halfedges[halfedges[CloseEdge].Twin].Next;
                faces[Face1].halfedge = halfedges[CloseEdge].Prev;

                halfedges[halfedges[CloseEdge].Prev].Next = halfedges[halfedges[CloseEdge].Twin].Next;
                halfedges[halfedges[halfedges[CloseEdge].Twin].Next].Prev = halfedges[CloseEdge].Prev;

                vertices[halfedges[halfedges[CloseEdge].Twin].vertex].Valid = false;
                //faces[Face1].Numvertices-=2;  //although one vertex should appear twice


                int hebegin = faces[Face1].halfedge;
                int heiterate = hebegin;
                do {

                    heiterate = halfedges[heiterate].Next;
                } while (heiterate != hebegin);

                hebegin = faces[Face1].halfedge;
                heiterate = hebegin;
                do {
                    heiterate = halfedges[heiterate].Next;
                } while (heiterate != hebegin);


                return true;
            }

            if (Face1 == Face2)
                return false;  //we don't remove non-spike edges with the same faces to not disconnect a chain

            hebegin = faces[Face2].halfedge;
            heiterate = hebegin;
            do {
                heiterate = halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            faces[Face1].halfedge = halfedges[heindex].Next;
            faces[Face2].Valid = false;

            //faces[Face2].halfedge=halfedges[halfedges[heindex].Twin].Next;

            halfedges[heindex].Valid = halfedges[halfedges[heindex].Twin].Valid = false;

            halfedges[halfedges[heindex].Next].Prev = halfedges[halfedges[heindex].Twin].Prev;
            halfedges[halfedges[halfedges[heindex].Twin].Prev].Next = halfedges[heindex].Next;

            halfedges[halfedges[halfedges[heindex].Twin].Next].Prev = halfedges[heindex].Prev;
            halfedges[halfedges[heindex].Prev].Next = halfedges[halfedges[heindex].Twin].Next;

            vertices[halfedges[heindex].vertex].halfedge = halfedges[halfedges[heindex].Twin].Next;
            vertices[halfedges[halfedges[heindex].Next].vertex].halfedge = halfedges[heindex].Next;

            //all other floating halfedges should renounce this one
            for (int i = 0; i < halfedges.size(); i++)
                if (halfedges[i].AdjFace == Face2)
                    halfedges[i].AdjFace = Face1;

            //faces[Face1].NumVertices+=faces[Face2].NumVertices-2;

            //DebugLog<<"Official number of vertices: "<<faces[Face1].NumVertices;

            hebegin = faces[Face1].halfedge;
            heiterate = hebegin;
            int currVertex = 0;
            do {
                //faces[Face1].Vertices[currVertex++]=halfedges[heiterate].vertex;
                heiterate = halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            return true;


        }

        void UnifyEdges(int heindex) {
            //if (halfedges[heindex].Twin<0)
            //  return;
            //adjusting source
            vertices[halfedges[heindex].vertex].Valid = false;
            halfedges[heindex].vertex = halfedges[halfedges[heindex].Prev].vertex;
            if (halfedges[heindex].prescribedAngle < 0.0)
                halfedges[heindex].prescribedAngle = halfedges[halfedges[heindex].Prev].prescribedAngle;
            vertices[halfedges[heindex].vertex].halfedge = heindex;

            faces[halfedges[heindex].face].halfedge = halfedges[heindex].Next;
            //faces[halfedges[heindex].face].NumVertices--;



            //adjusting halfedges
            halfedges[halfedges[heindex].Prev].Valid = false;
            halfedges[heindex].Prev = halfedges[halfedges[heindex].Prev].Prev;
            halfedges[halfedges[heindex].Prev].Next = heindex;

            //adjusting twin, if exists
            if (halfedges[heindex].Twin >= 0) {
                if (halfedges[halfedges[heindex].Twin].prescribedAngle < 0.0)
                    halfedges[halfedges[heindex].Twin].prescribedAngle = halfedges[halfedges[halfedges[heindex].Twin].Next].prescribedAngle;
                halfedges[halfedges[halfedges[heindex].Twin].Next].Valid = false;
                halfedges[halfedges[heindex].Twin].Next = halfedges[halfedges[halfedges[heindex].Twin].Next].Next;
                halfedges[halfedges[halfedges[heindex].Twin].Next].Prev = halfedges[heindex].Twin;
                faces[halfedges[halfedges[heindex].Twin].face].halfedge = halfedges[halfedges[heindex].Twin].Next;
                //faces[halfedges[halfedges[heindex].Twin].face].NumVertices--;
            }
        }


        void aggregage_dcel(DCEL<VertexData, HalfedgeData, EdgeData, FaceData>)
        {
            //TODO: aggregate new elements and reindex where necessary
        }

        // Created a Double-Connected Edge-List (a.k.a. "halfedge structure") from the usual
        // libhedra mesh representation. This data structure is very convenient for mesh editing
        // and traversing, and the data structure is again only Eigen vectors and matrices.

        //input:
        //  D           #F by 1 - face degrees
        //  F           #F by max(D) - vertex indices in face
        //  EV          #E by 2 - edge vertex indices
        //  EF          #E by 2 - edge face indices (EF(i,0) is left face, EF(i,1)=-1 if boundary
        //  EFi         #E by 2 - position of edge in face by EF
        //  innerEdges  vector of inner edges into EV

        // Output:
        // the number of halfedges can be determined by H=|HV/HE/HF|. It is 2*[Inner edges]+[Boundary Edges]
        // VH   #V by 1 - Vertex to outgoing halfedge (into HE)
        // EH   #E by 2 - edge to halfedge, where EH(i,0) halfedge is positively oriented, and EH(i,1)=-1 when boundary.
        // FH   #F by max(D) - face to (correctly oriented) halfedge s.t. the origin vertex of FH(i,j) is F(i,j)
        // HV   #H by 1 - origin vertex of the halfedge
        // HE   #H by 1 - edge carrying this halfedge. It does not say which direction.
        // HF   #F by 1 - face containing halfedge
        // nextH, prevH, twinH - #H by 1 DCEL traversing operations. halfedge[i].twin=-1 for boundary edges.

        /*DCEL(const Eigen::VectorXi& D,
             const Eigen::MatrixXi& F,
             const Eigen::MatrixXi& EV,
             const Eigen::MatrixXi& EF,
             const Eigen::MatrixXi& EFi,
             const Eigen::VectorXi& innerEdges,
             Eigen::VectorXi& VH,
             Eigen::MatrixXi& EH,
             Eigen::MatrixXi& FH,
             Eigen::VectorXi& HV,
             Eigen::VectorXi& HE,
             Eigen::VectorXi& HF,
             Eigen::VectorXi& nextH,
             Eigen::VectorXi& prevH,
             Eigen::VectorXi& twinH)
        {
            //doing a local halfedge structure for polygonal meshes
            EH=Eigen::MatrixXi::Constant(EV.rows(),2,-1);
            int numH=0;

            for (int i=0;i<EF.rows();i++){
                if (EF(i,0)!=-1)
                    EH(i,0)=numH++;
                if (EF(i,1)!=-1)
                    EH(i,1)=numH++;
            }


            //halfedges to edge
            HE.conservativeResize(numH);
            for (int i=0;i<EH.rows();i++){
                if (EH(i,0)!=-1)
                    HE(EH(i,0))=i;
                if (EH(i,1)!=-1)
                    HE(EH(i,1))=i;
            }

            //halfedge to vertex and vice versa
            HV.conservativeResize(numH);
            VH.conservativeResize(EV.maxCoeff()+1);
            for (int i=0;i<EV.rows();i++){
                if (EH(i,0)!=-1){
                    HV(EH(i,0))=EV(i,0);
                    VH(EV(i,0))=EH(i,0);
                }
                if (EH(i,1)!=-1){
                    HV(EH(i,1))=EV(i,1);
                    VH(EV(i,1))=EH(i,1);
                }
            }

            //halfedge to twin
            twinH=Eigen::VectorXi::Constant(numH, -1);
            for (int i=0;i<EH.rows();i++)
                if ((EH(i,0)!=-1)&&(EH(i,1)!=-1)){
                    twinH(EH(i,0))=EH(i,1);
                    twinH(EH(i,1))=EH(i,0);
                }

            //faces to halfedges and vice versa
            FH.resize(F.rows(), F.cols());
            HF.resize(numH);
            for (int i=0;i<EF.rows();i++){
                if (EF(i,0)!=-1){
                    FH(EF(i,0),EFi(i,0))=EH(i,0);
                    HF(EH(i,0))=EF(i,0);
                }
                if (EF(i,1)!=-1){
                    FH(EF(i,1),EFi(i,1))=EH(i,1);
                    HF(EH(i,1))=EF(i,1);
                }
            }

            //halfedge to next and prev
            nextH.conservativeResize(HE.rows());
            prevH.conservativeResize(HE.rows());
            for (int i=0;i<D.rows();i++){
                for (int j=0;j<D(i);j++){
                    nextH(FH(i,j))=FH(i,(j+1)%D(i));
                    prevH(FH(i,(j+1)%D(i)))=FH(i,j);
                }
            }

        }*/


    };

}


#endif
