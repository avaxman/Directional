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
#include <deque>


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
                        std::cout << "Valid Vertex " << i << " points to non-valid halfedge "
                                  << vertices[i].halfedge << std::endl;
                    return false;
                }


                if (halfedges[vertices[i].halfedge].vertex != i) {
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
                        std::cout << "Valid Halfedge " << i << "points to next non-valid value -1" << std::endl;
                    return false;
                }

                if (halfedges[i].prev == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to prev non-valid value -1" << std::endl;
                    return false;
                }


                if (halfedges[i].vertex == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Origin non-valid value -1" << std::endl;
                    return false;
                }

                if (halfedges[i].face == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to face non-valid value -1" << std::endl;
                    return false;
                }

                if (halfedges[halfedges[i].next].prev != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "next is " << halfedges[i].next
                                  << " which doesn't point back as prev" << std::endl;
                    return false;
                }


                if (halfedges[halfedges[i].prev].next != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "prev is " << halfedges[i].prev
                                  << " which doesn't point back as next" << std::endl;
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

                if (!edges[halfedges[i].edge].valid) {
                    if (verbose)
                        std::cout << "The edge of halfedges " << i << ", edge " << halfedges[i].edge
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

            //checking edges
            for (int i=0;i<edges.size();i++){
                if (!edges[i].valid)
                    continue;

                if (halfedges[edges[i].halfedge].edge!=i){
                    std::cout<<"Edge "<<i<<" points to halfedge "<<edges[i].halfedge<<" but not back."<<std::endl;
                    return false;
                }

                if (!halfedges[edges[i].halfedge].valid){
                    std::cout<<"Edge "<<i<<" points to halfedge "<<edges[i].halfedge<<" which is not valid."<<std::endl;
                    return false;
                }
                if ((halfedges[edges[i].halfedge].twin!=-1)&&(halfedges[halfedges[edges[i].halfedge].twin].edge!=i)){
                    std::cout<<"Edge "<<i<<" with halfedge "<<edges[i].halfedge<<" points to twin halfedge "<<halfedges[edges[i].halfedge].twin<<" which does not share the same edge."<<std::endl;
                    return false;
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
                    if (!halfedges[i].valid)
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
                            std::cout << "Face " << halfedges[i].face << " is a pure boundary face!" << std::endl;
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
                if (!halfedges[i].valid)
                    continue;

                Halfedge NewHalfedge = halfedges[i];
                NewHalfedge.ID = newHalfedges.size();
                newHalfedges.push_back(NewHalfedge);
                transHalfedges[i] = NewHalfedge.ID;
            }

            halfedges = newHalfedges;
            for (int i = 0; i < faces.size(); i++)
                faces[i].halfedge = transHalfedges[faces[i].halfedge];

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
            std::vector<Edge> newEdges;
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
                        TwinFinder(0, halfedges[halfedges[i].next].vertex, halfedges[i].vertex));
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
            if (halfedges[heindex].twin < 0)
                return true;  //there is no joining of boundary faces

            int Face1 = halfedges[heindex].face;
            int Face2 = halfedges[halfedges[heindex].twin].face;

            /*int hebegin = faces[Face1].halfedge;
            int heiterate = hebegin;
            do {
                heiterate = halfedges[heiterate].next;
            } while (heiterate != hebegin);

            hebegin = faces[Face1].halfedge;
            heiterate = hebegin;
            do {
                heiterate = halfedges[heiterate].next;
            } while (heiterate != hebegin);*/

            //check if spike edge
            if ((halfedges[heindex].prev == halfedges[heindex].twin) ||
                (halfedges[heindex].next == halfedges[heindex].twin)) {


                int CloseEdge = heindex;
                if (halfedges[heindex].prev == halfedges[heindex].twin)
                    CloseEdge = halfedges[heindex].twin;

                halfedges[CloseEdge].valid = halfedges[halfedges[CloseEdge].twin].valid = edges[halfedges[CloseEdge].edge] = false;

                vertices[halfedges[CloseEdge].vertex].halfedge = halfedges[halfedges[CloseEdge].twin].next;
                faces[Face1].halfedge = halfedges[CloseEdge].prev;

                halfedges[halfedges[CloseEdge].prev].next = halfedges[halfedges[CloseEdge].twin].next;
                halfedges[halfedges[halfedges[CloseEdge].twin].next].prev = halfedges[CloseEdge].prev;

                vertices[halfedges[halfedges[CloseEdge].twin].vertex].valid = false;
                //faces[Face1].Numvertices-=2;  //although one vertex should appear twice


                int hebegin = faces[Face1].halfedge;
                int heiterate = hebegin;
                do {

                    heiterate = halfedges[heiterate].next;
                } while (heiterate != hebegin);

                hebegin = faces[Face1].halfedge;
                heiterate = hebegin;
                do {
                    heiterate = halfedges[heiterate].next;
                } while (heiterate != hebegin);


                return true;
            }

            if (Face1 == Face2)
                return false;  //we don't remove non-spike edges with the same faces to not disconnect a chain

            /*int hebegin = faces[Face2].halfedge;
            int heiterate = hebegin;
            do {
                heiterate = halfedges[heiterate].next;
            } while (heiterate != hebegin);*/

            faces[Face1].halfedge = halfedges[heindex].next;
            faces[Face2].valid = false;

            //faces[Face2].halfedge=halfedges[halfedges[heindex].twin].next;

            halfedges[heindex].valid = halfedges[halfedges[heindex].twin].valid = false;

            halfedges[halfedges[heindex].next].prev = halfedges[halfedges[heindex].twin].prev;
            halfedges[halfedges[halfedges[heindex].twin].prev].next = halfedges[heindex].next;

            halfedges[halfedges[halfedges[heindex].twin].next].prev = halfedges[heindex].prev;
            halfedges[halfedges[heindex].prev].next = halfedges[halfedges[heindex].twin].next;

            vertices[halfedges[heindex].vertex].halfedge = halfedges[halfedges[heindex].twin].next;
            vertices[halfedges[halfedges[heindex].next].vertex].halfedge = halfedges[heindex].next;

            //all other floating halfedges should renounce this one
            for (int i = 0; i < halfedges.size(); i++)
                if (halfedges[i].face == Face2)
                    halfedges[i].face = Face1;

            //faces[Face1].NumVertices+=faces[Face2].NumVertices-2;

            //DebugLog<<"Official number of vertices: "<<faces[Face1].NumVertices;

            /*hebegin = faces[Face1].halfedge;
            heiterate = hebegin;
            int currVertex = 0;
            do {
                //faces[Face1].Vertices[currVertex++]=halfedges[heiterate].vertex;
                heiterate = halfedges[heiterate].next;
            } while (heiterate != hebegin);*/

            return true;


        }

        void unify_edges(int heindex) {
            //if (halfedges[heindex].twin<0)
            //  return;
            //adjusting source

            //std::cout<<"Unifying halfedge "<<halfedges[heindex].prev<<" into halfedge "<<heindex<<" killing edge "<<halfedges[halfedges[heindex].prev].edge<<std::endl;
            vertices[halfedges[heindex].vertex].valid = false;
            halfedges[heindex].vertex = halfedges[halfedges[heindex].prev].vertex;
            if (halfedges[heindex].data.prescribedAngle < 0.0)
                halfedges[heindex].data.prescribedAngle = halfedges[halfedges[heindex].prev].data.prescribedAngle;
            vertices[halfedges[heindex].vertex].halfedge = heindex;

            faces[halfedges[heindex].face].halfedge = halfedges[heindex].next;
            //faces[halfedges[heindex].face].NumVertices--;



            //adjusting halfedges
            halfedges[halfedges[heindex].prev].valid = false;
            edges[halfedges[halfedges[heindex].prev].edge].valid = false;
            halfedges[heindex].prev = halfedges[halfedges[heindex].prev].prev;
            halfedges[halfedges[heindex].prev].next = heindex;

            //adjusting twin, if exists
            if (halfedges[heindex].twin >= 0) {
                //if (halfedges[halfedges[heindex].twin].data.prescribedAngle < 0.0)
                //    halfedges[halfedges[heindex].twin].data.prescribedAngle = halfedges[halfedges[halfedges[heindex].twin].next].data.prescribedAngle;
                //std::cout<<"Unifying halfedge "<<halfedges[halfedges[heindex].twin].next<<" into halfedge "<<halfedges[heindex].twin<<" killing edge "<<halfedges[halfedges[halfedges[heindex].twin].next].edge<<std::endl;
                halfedges[halfedges[halfedges[heindex].twin].next].valid = false;
                edges[halfedges[halfedges[halfedges[heindex].twin].next].edge].valid = false;
                halfedges[halfedges[heindex].twin].next = halfedges[halfedges[halfedges[heindex].twin].next].next;
                halfedges[halfedges[halfedges[heindex].twin].next].prev = halfedges[heindex].twin;
                faces[halfedges[halfedges[heindex].twin].face].halfedge = halfedges[halfedges[heindex].twin].next;
                //faces[halfedges[halfedges[heindex].twin].face].NumVertices--;
            }
        }


        void aggregate_dcel(const DCEL<VertexData, HalfedgeData, EdgeData, FaceData>& aggDcel)
        {
            //TODO: aggregate new elements and reindex where necessary
            int currVOffset = vertices.size();
            int currHEOffset = halfedges.size();
            int currEOffset = edges.size();
            int currFOffset = faces.size();
            for (int i=0;i<aggDcel.vertices.size();i++){
                vertices.push_back(aggDcel.vertices[i]);
                vertices[vertices.size()-1].ID += currVOffset;
                vertices[vertices.size()-1].halfedge += currHEOffset;
            }

            for (int i=0;i<aggDcel.halfedges.size();i++){
                halfedges.push_back(aggDcel.halfedges[i]);
                halfedges[halfedges.size()-1].ID += currHEOffset;
                halfedges[halfedges.size()-1].vertex += currVOffset;
                halfedges[halfedges.size()-1].next += currHEOffset;
                halfedges[halfedges.size()-1].prev += currHEOffset;
                if (halfedges[halfedges.size()-1].twin!=-1)
                    halfedges[halfedges.size()-1].twin+= currHEOffset;
                halfedges[halfedges.size()-1].face += currFOffset;
                halfedges[halfedges.size()-1].edge += currEOffset;
            }

            for (int i=0;i<aggDcel.edges.size();i++){
                edges.push_back(aggDcel.edges[i]);
                edges[edges.size()-1].ID += currEOffset;
                edges[edges.size()-1].halfedge += currHEOffset;
            }
            for (int i=0;i<aggDcel.faces.size();i++){
                faces.push_back(aggDcel.faces[i]);
                faces[faces.size()-1].ID += currFOffset;
                faces[faces.size()-1].halfedge += currHEOffset;
            }
        }

        void RemoveVertex(int vindex, std::deque<int> &removeVertexQueue) {
            int hebegin = vertices[vindex].halfedge;
            int heiterate = hebegin;
            do {
                if (heiterate == -1) {  //boundary vertex
                    return;
                }
                heiterate = halfedges[halfedges[heiterate].prev].twin;
            } while (heiterate != hebegin);

            vertices[vindex].valid = false;

            int remainingFace = halfedges[hebegin].face;


            faces[remainingFace].halfedge = halfedges[hebegin].next;
            heiterate = hebegin;
            int infinityCounter = 0;
            do {

                int NextEdge = halfedges[heiterate].next;
                int PrevEdge = halfedges[halfedges[heiterate].twin].prev;

                halfedges[NextEdge].prev = PrevEdge;
                halfedges[PrevEdge].next = NextEdge;
                if (halfedges[NextEdge].face != remainingFace)
                    faces[halfedges[NextEdge].face].valid = false;

                if (halfedges[PrevEdge].face != remainingFace)
                    faces[halfedges[PrevEdge].face].valid = false;


                halfedges[PrevEdge].face = halfedges[NextEdge].face = remainingFace;
                halfedges[heiterate].valid = false;
                halfedges[halfedges[heiterate].twin].valid = false;
                heiterate = halfedges[halfedges[heiterate].prev].twin;
                infinityCounter++;
                if (infinityCounter > halfedges.size())
                    return;

            } while (heiterate != hebegin);

            //cleaning new face
            hebegin = faces[remainingFace].halfedge;
            //faces[remainingFace].Numvertices=0;
            heiterate = hebegin;
            infinityCounter = 0;
            do {
                //faces[remainingFace].Numvertices++;
                halfedges[heiterate].face = remainingFace;
                vertices[halfedges[heiterate].Origin].halfedge = heiterate;
                removeVertexQueue.push_front(halfedges[heiterate].Origin);
                infinityCounter++;
                if (infinityCounter > halfedges.size())
                    return;

                heiterate = halfedges[heiterate].next;
            } while (heiterate != hebegin);
        }



        //Initializing DCEL from faces, assuming this is a triangle mesh
        void init(const Eigen::MatrixXd& V,
                  const Eigen::MatrixXi& F){

            halfedges.resize(3*F.rows());
            vertices.resize(V.rows());
            faces.resize(F.rows());
            edges.clear();

            for (int i=0;i<F.rows();i++) {
                faces[i].ID=i;
                faces[i].halfedge=3*i;
                for (int j=0;j<3;j++){
                    halfedges[3*i+j].ID = 3*i+j;
                    halfedges[3*i+j].vertex=F(i,j);
                    vertices[halfedges[3*i+j].vertex].halfedge=3*i+j;
                    halfedges[3*i+j].next=3*i+(j+1)%3;
                    halfedges[3*i+j].prev=3*i+(j+2)%3;
                    halfedges[3*i+j].face=i;
                }
            }

            for (int i=0;i<vertices.size();i++)
                vertices[i].ID=i;

            struct ComparePairs {
                bool operator()(const std::pair<std::pair<int, int>, int>& a, const std::pair<std::pair<int, int>, int>& b) const {
                    if (a.first.first == b.first.first) {
                        return a.first.second < b.first.second;
                    } else {
                        return a.first.first < b.first.first;
                    }
                }
            };

            //finding twins
            typedef std::pair<std::pair<int, int>, int> pairPlusOne;
            std::set<pairPlusOne, ComparePairs> edgeSet;
            std::vector<int> EHList;
            for (int i=0;i<halfedges.size();i++){
                std::pair<int,int> oppEdge(halfedges[halfedges[i].next].vertex, halfedges[i].vertex);
                pairPlusOne oppEdgePlus(oppEdge, -1);
                std::set<pairPlusOne>::iterator si = edgeSet.find(oppEdgePlus);
                if (si == edgeSet.end()) {
                    edgeSet.insert(pairPlusOne(std::pair<int, int>( halfedges[i].vertex, halfedges[halfedges[i].next].vertex), i));
                    EHList.push_back(i);
                } else {  //found matching twin
                    halfedges[si->second].twin = i;
                    halfedges[i].twin = si->second;
                }
            }

            //creating edges
            for (int i=0;i<halfedges.size();i++) {
                if (halfedges[i].edge!=-1)
                    continue;

                edges.push_back(Edge());
                edges[edges.size()-1].ID=edges.size()-1;
                edges[edges.size()-1].halfedge=i;
                halfedges[i].edge=edges.size()-1;
                if (halfedges[i].twin!=-1)
                    halfedges[halfedges[i].twin].edge=edges.size()-1;

            }
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
