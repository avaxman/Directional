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
    class DCEL{
    public:

        Eigen::VectorXi VH, EH, FH;
        Eigen::VectorXi HV, HE, HF;
        Eigen::VectorXi nextH, prevH, twinH;

        Eigen::VectorXi VValid, HValid, FValid;

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

        bool check_consistency(const bool verbose, const bool checkHalfedgeRepetition, const bool CheckTwinGaps,
                       const bool checkPureBoundary) {
            for (int i = 0; i < VH.size(); i++) {
                if (!VValid(i))
                    continue;

                if (VH(i) == -1) {
                    if (verbose) std::cout << "Valid Vertex " << i << " points to non-valid value -1 " << std::endl;
                    return false;
                }

                if (!HValid[VH(i)]) {
                    if (verbose)
                        std::cout << "Valid Vertex " << i << " points to non-valid AdjHalfedge "
                                  << VH(i) << std::endl;
                    return false;
                }


                if (HV[VH(i)] != i) {
                    if (verbose)
                        std::cout << "Adjacent Halfedge " << VH(i) << " of vertex " << i
                                  << "does not point back" << std::endl;
                    return false;
                }

            }

            for (int i = 0; i < HV.size(); i++) {
                if (!HValid[i])
                    continue;


                if (nextH[i] == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Next non-valid value -1" << std::endl;
                    return false;
                }

                if (prevH[i] == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Prev non-valid value -1" << std::endl;
                    return false;
                }


                if (HV[i] == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to Origin non-valid value -1" << std::endl;
                    return false;
                }

                if (HF[i] == -1) {
                    if (verbose)
                        std::cout << "Valid Halfedge " << i << "points to AdjFace non-valid value -1" << std::endl;
                    return false;
                }

                if (prevH[nextH[i]] != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "Next is " << nextH[i]
                                  << " which doesn't point back as Prev" << std::endl;
                    return false;
                }


                if (nextH[prevH[i]] != i) {
                    if (verbose)
                        std::cout << "Halfedge " << i << "Prev is " << prevH[i]
                                  << " which doesn't point back as Next" << std::endl;
                    return false;
                }

                if (!VValid[HV[i]]) {
                    if (verbose)
                        std::cout << "The Origin of halfedges " << i << ", vertex " << HV[i]
                                  << " is not valid" << std::endl;
                    return false;
                }

                if (!FValid[HF[i]]) {
                    if (verbose)
                        std::cout << "The face of halfedges " << i << ", face " << HF[i]
                                  << " is not valid" << std::endl;
                    return false;
                }

                if (HV[nextH[i]] == HV[i]) {  //a degenerate edge{
                    if (verbose)
                        std::cout << "Halfedge " << i << " with twin" << twinH[i]
                                  << " is degenerate with vertex " << HV[i] << std::endl;
                    return false;
                }

                if (twinH[i] >= 0) {
                    if (twinH[twinH[i]] != i) {
                        if (verbose)
                            std::cout << "Halfedge " << i << "twin is " << twinH[i]
                                      << " which doesn't point back" << std::endl;
                        return false;
                    }

                    if (!HValid[twinH[i]]) {
                        if (verbose)
                            std::cout << "halfedge " << i << " is twin with invalid halfedge" << twinH[i]
                                      << std::endl;
                        return false;
                    }
                }

                if (!HValid[nextH[i]]) {
                    if (verbose)
                        std::cout << "halfedge " << i << " has next invalid halfedge" << nextH[i] << std::endl;
                    return false;
                }

                if (!HValid[prevH[i]]) {
                    if (verbose)
                        std::cout << "halfedge " << i << " has prev invalid halfedge" << prevH[i] << std::endl;
                    return false;
                }

                //if (Halfedges[i].isFunction) {  //checking that it is not left alone
                if (prevH[i] == twinH[i]) {
                    if (verbose)
                        std::cout << "Hex halfedge " << i << " has Halfedge " << prevH[i]
                                  << " and both prev and twin" << std::endl;
                    return false;
                }


                if (nextH[i] == twinH[i]) {
                    if (verbose)
                        std::cout << "Hex halfedge " << i << " has Halfedge " << nextH[i]
                                  << " and both next and twin" << std::endl;
                    return false;
                }
                //}
            }

            std::vector <std::set<int>> halfedgesinFace(FH.size());
            std::vector <std::set<int>> verticesinFace(FH.size());
            for (int i = 0; i < FH.size(); i++) {
                if (!FValid[i])
                    continue;

                //if (Faces[i].NumVertices<3)  //we never allow this
                //  return false;
                int hebegin = FH[i];
                int heiterate = hebegin;
                int NumEdges = 0;
                int actualNumVertices = 0;

                do {
                    if (verticesinFace[i].find(HV[heiterate]) != verticesinFace[i].end())
                        if (verbose)
                            std::cout << "Warning: Vertex " << HV[heiterate]
                                      << " appears more than once in face " << i << std::endl;

                    verticesinFace[i].insert(HV[heiterate]);
                    halfedgesinFace[i].insert(heiterate);
                    actualNumVertices++;
                    if (!HValid[heiterate])
                        return false;

                    if (HF[heiterate] != i) {
                        if (verbose)
                            std::cout << "Face " << i << " has halfedge " << heiterate << " that does not point back"
                                      << std::endl;
                        return false;
                    }

                    heiterate = nextH[heiterate];
                    NumEdges++;
                    if (NumEdges > HV.size()) {
                        if (verbose) std::cout << "Infinity loop!" << std::endl;
                        return false;
                    }


                } while (heiterate != hebegin);

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
            for (int i = 0; i < HV.size(); i++) {
                if (!HValid[i])
                    continue;
                int currFace = HF[i];
                if (halfedgesinFace[currFace].find(i) == halfedgesinFace[currFace].end()) {
                    if (verbose) std::cout << "Halfedge " << i << " is floating in face " << currFace << std::endl;
                    return false;
                }
            }

            //check if mesh is a manifold: every halfedge appears only once
            if (checkHalfedgeRepetition) {
                std::set <TwinFinder> HESet;
                for (int i = 0; i < HV.size(); i++) {
                    if (!HValid[i].Valid)
                        continue;
                    std::set<TwinFinder>::iterator HESetIterator = HESet.find(
                            TwinFinder(i, HV[i], HV[nextH[i]]));
                    if (HESetIterator != HESet.end()) {
                        if (verbose)
                            std::cout << "Warning: the halfedge (" << HV[i] << ","
                                      << HV[nextH[i]] << ") appears at least twice in the mesh"
                                      << std::endl;
                        if (verbose)
                            std::cout << "for instance halfedges " << i << " and " << HESetIterator->index << std::endl;
                        return false;
                        //return false;
                    } else {
                        HESet.insert(TwinFinder(i, HV[i], HV[nextH[i]]));
                        //if (verbose) std::cout<<"inserting halfedge "<<i<<" which is "<<Halfedges[i].Origin<<", "<<Halfedges[nextH(i)].Origin<<endl;
                    }
                }
            }

            if (CheckTwinGaps) {
                std::set <TwinFinder> HESet;
                //checking if there is a gap: two halfedges that share the same opposite vertices but do not have twins
                for (int i = 0; i < HV.size(); i++) {
                    if (!HValid[i])
                        continue;

                    std::set<TwinFinder>::iterator HESetIterator = HESet.find(
                            TwinFinder(i, HV[i], HV[nextH[i]]));
                    if (HESetIterator == HESet.end()) {
                        HESet.insert(TwinFinder(i, HV[i], HV[nextH[i]]));
                        continue;
                    }

                    HESetIterator = HESet.find(TwinFinder(i, HV[nextH[i]], HV[i]));
                    if (HESetIterator != HESet.end()) {

                        if (twinH(i) == -1) {
                            if (verbose)
                                std::cout << "Halfedge " << i << "has no twin although halfedge "
                                          << HESetIterator->index << " can be a twin" << std::endl;
                            return false;
                        }
                        if (twinH[HESetIterator->index] == -1) {
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
                for (int i = 0; i < HV.size(); i++) {
                    if (!HValid[i])
                        continue;

                    //Is this necessary?
                    /*if ((twinH(i) < 0) && (Halfedges[i].isFunction))
                        if (verbose)
                            std::cout << "WARNING: Halfedge " << i << " is a hex edge without twin!" << std::endl;
*/
                    if (twinH(i) > 0)
                        continue;

                    bool pureBoundary = true;
                    int hebegin = i;
                    int heiterate = hebegin;
                    do {
                        if (twinH[heiterate] > 0) {
                            pureBoundary = false;
                            break;
                        }
                        heiterate = nextH[heiterate];
                    } while (heiterate != hebegin);
                    if (pureBoundary) {
                        if (verbose)
                            std::cout << "Face " << HF[i] << "is a pure boundary face!" << std::endl;
                        return false;
                    }
                }

                //checking for latent valence 2 faces
                std::vector<int> Valences(VH.size());
                for (int i = 0; i < VH.size(); i++)
                    Valences[i] = 0;

                for (int i = 0; i < HV.size(); i++) {
                    if (HValid[i]) {
                        Valences[HV[i]]++;
                        //Valences[Halfedges[nextH(i)].Origin]++;
                        if (twinH(i) < 0)  //should account for the target as well
                            Valences[HV[nextH(i)]]++;
                    }
                }

                int countThree;
                for (int i = 0; i < FH.size(); i++) {
                    if (!FValid[i])
                        continue;
                    countThree = 0;
                    int hebegin = FH[i];
                    int heiterate = hebegin;
                    int numEdges = 0;
                    do {
                        if (Valences[HV[heiterate]] > 2)
                            countThree++;
                        heiterate = nextH[heiterate];
                        numEdges++;
                        if (numEdges > HV.size()) {
                            if (verbose) std::cout << "Infinity loop in face " << i << "!" << std::endl;
                            return false;
                        }
                    } while (heiterate != hebegin);
                    if (countThree < 3) {
                        if (verbose) std::cout << "Face " << i << " is a latent valence 2 face!" << std::endl;
                        if (verbose) std::cout << "Its vertices are " << std::endl;
                        do {
                            if (verbose)
                                std::cout << "Vertex " << HV[heiterate] << " halfedge " << heiterate
                                          << " valence " << Valences[HV[heiterate]] << std::endl;

                            if (Valences[HV[heiterate]] > 2)
                                countThree++;
                            heiterate = nextH[heiterate];
                            numEdges++;
                            if (numEdges > HV.size()) {
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

        //TODO: convert to DCEL format
        void clean_mesh() {
            //removing nonvalid vertices
            std::vector<int> TransVertices(VH.size());
            std::vector <Vertex> NewVertices;
            for (int i = 0; i < Vertices.size(); i++) {
                if (!Vertices[i].Valid)
                    continue;

                Vertex NewVertex = Vertices[i];
                NewVertex.ID = NewVertices.size();
                NewVertices.push_back(NewVertex);
                TransVertices[i] = NewVertex.ID;
            }


            Vertices = NewVertices;
            for (int i = 0; i < Halfedges.size(); i++)
                Halfedges[i].Origin = TransVertices[Halfedges[i].Origin];



            /*for (int i=0;i<Faces.size();i++){
              if (!Faces[i].Valid)
                continue;
              for (int j=0;j<Faces[i].NumVertices;j++){
                DebugLog<<"i is "<<i<<", j is "<<j<<", Faces[i].Vertices[j] is "<<Faces[i].Vertices[j]<<endl;
                Faces[i].Vertices[j]=TransVertices[Faces[i].Vertices[j]];
              }
            }*/



            //removing nonvalid faces
            std::vector <Face> NewFaces;
            std::vector<int> TransFaces(Faces.size());
            for (int i = 0; i < Faces.size(); i++) {
                if (!Faces[i].Valid)
                    continue;

                Face NewFace = Faces[i];
                NewFace.ID = NewFaces.size();
                NewFaces.push_back(NewFace);
                TransFaces[i] = NewFace.ID;
            }
            Faces = NewFaces;
            for (int i = 0; i < Halfedges.size(); i++)
                Halfedges[i].AdjFace = TransFaces[Halfedges[i].AdjFace];



            //removing nonvalid halfedges
            std::vector <Halfedge> NewHalfedges;
            std::vector<int> TransHalfedges(Halfedges.size());
            for (int i = 0; i < Halfedges.size(); i++) {
                if (!Halfedges[i].Valid)
                    continue;

                Halfedge NewHalfedge = Halfedges[i];
                NewHalfedge.ID = NewHalfedges.size();
                NewHalfedges.push_back(NewHalfedge);
                TransHalfedges[i] = NewHalfedge.ID;
            }


            Halfedges = NewHalfedges;
            for (int i = 0; i < Faces.size(); i++)
                Faces[i].AdjHalfedge = TransHalfedges[Faces[i].AdjHalfedge];


            for (int i = 0; i < Vertices.size(); i++)
                VH(i) = TransHalfedges[VH(i)];


            for (int i = 0; i < Halfedges.size(); i++) {
                if (twinH(i) != -1)
                    twinH(i) = TransHalfedges[twinH(i)];
                nextH(i) = TransHalfedges[nextH(i)];
                prevH(i) = TransHalfedges[prevH(i)];
            }

        }

        void ComputeTwins() {
            //twinning up edges
            std::set <TwinFinder> Twinning;
            for (int i = 0; i < Halfedges.size(); i++) {
                if (Halfedges[i].Twin >= 0)
                    continue;

                std::set<TwinFinder>::iterator Twinit = Twinning.find(
                        TwinFinder(0, Halfedges[Halfedges[i].Next].Origin, Halfedges[i].Origin));
                if (Twinit != Twinning.end()) {
                    Halfedges[Twinit->index].Twin = i;
                    Halfedges[i].Twin = Twinit->index;
                    Twinning.erase(*Twinit);
                } else {
                    Twinning.insert(TwinFinder(i, Halfedges[i].Origin, Halfedges[Halfedges[i].Next].Origin));
                }
            }
        }

        bool JoinFace(const int heindex) {
            if (Halfedges[heindex].Twin < 0)
                return true;  //there is no joining of boundary faces

            int Face1 = Halfedges[heindex].AdjFace;
            int Face2 = Halfedges[Halfedges[heindex].Twin].AdjFace;

            int hebegin = Faces[Face1].AdjHalfedge;
            int heiterate = hebegin;
            do {
                heiterate = Halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            hebegin = Faces[Face1].AdjHalfedge;
            heiterate = hebegin;
            do {
                heiterate = Halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            //check if spike edge
            if ((Halfedges[heindex].Prev == Halfedges[heindex].Twin) ||
                (Halfedges[heindex].Next == Halfedges[heindex].Twin)) {


                int CloseEdge = heindex;
                if (Halfedges[heindex].Prev == Halfedges[heindex].Twin)
                    CloseEdge = Halfedges[heindex].Twin;

                Halfedges[CloseEdge].Valid = Halfedges[Halfedges[CloseEdge].Twin].Valid = false;

                Vertices[Halfedges[CloseEdge].Origin].AdjHalfedge = Halfedges[Halfedges[CloseEdge].Twin].Next;
                Faces[Face1].AdjHalfedge = Halfedges[CloseEdge].Prev;

                Halfedges[Halfedges[CloseEdge].Prev].Next = Halfedges[Halfedges[CloseEdge].Twin].Next;
                Halfedges[Halfedges[Halfedges[CloseEdge].Twin].Next].Prev = Halfedges[CloseEdge].Prev;

                Vertices[Halfedges[Halfedges[CloseEdge].Twin].Origin].Valid = false;
                //Faces[Face1].NumVertices-=2;  //although one vertex should appear twice


                int hebegin = Faces[Face1].AdjHalfedge;
                int heiterate = hebegin;
                do {

                    heiterate = Halfedges[heiterate].Next;
                } while (heiterate != hebegin);

                hebegin = Faces[Face1].AdjHalfedge;
                heiterate = hebegin;
                do {
                    heiterate = Halfedges[heiterate].Next;
                } while (heiterate != hebegin);


                return true;
            }

            if (Face1 == Face2)
                return false;  //we don't remove non-spike edges with the same faces to not disconnect a chain

            hebegin = Faces[Face2].AdjHalfedge;
            heiterate = hebegin;
            do {
                heiterate = Halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            Faces[Face1].AdjHalfedge = Halfedges[heindex].Next;
            Faces[Face2].Valid = false;

            //Faces[Face2].AdjHalfedge=Halfedges[Halfedges[heindex].Twin].Next;

            Halfedges[heindex].Valid = Halfedges[Halfedges[heindex].Twin].Valid = false;

            Halfedges[Halfedges[heindex].Next].Prev = Halfedges[Halfedges[heindex].Twin].Prev;
            Halfedges[Halfedges[Halfedges[heindex].Twin].Prev].Next = Halfedges[heindex].Next;

            Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev = Halfedges[heindex].Prev;
            Halfedges[Halfedges[heindex].Prev].Next = Halfedges[Halfedges[heindex].Twin].Next;

            Vertices[Halfedges[heindex].Origin].AdjHalfedge = Halfedges[Halfedges[heindex].Twin].Next;
            Vertices[Halfedges[Halfedges[heindex].Next].Origin].AdjHalfedge = Halfedges[heindex].Next;

            //all other floating halfedges should renounce this one
            for (int i = 0; i < Halfedges.size(); i++)
                if (Halfedges[i].AdjFace == Face2)
                    Halfedges[i].AdjFace = Face1;

            //Faces[Face1].NumVertices+=Faces[Face2].NumVertices-2;

            //DebugLog<<"Official number of vertices: "<<Faces[Face1].NumVertices;

            hebegin = Faces[Face1].AdjHalfedge;
            heiterate = hebegin;
            int currVertex = 0;
            do {
                //Faces[Face1].Vertices[currVertex++]=Halfedges[heiterate].Origin;
                heiterate = Halfedges[heiterate].Next;
            } while (heiterate != hebegin);

            return true;


        }

        void UnifyEdges(int heindex) {
            //if (Halfedges[heindex].Twin<0)
            //  return;
            //adjusting source
            Vertices[Halfedges[heindex].Origin].Valid = false;
            Halfedges[heindex].Origin = Halfedges[Halfedges[heindex].Prev].Origin;
            if (Halfedges[heindex].prescribedAngle < 0.0)
                Halfedges[heindex].prescribedAngle = Halfedges[Halfedges[heindex].Prev].prescribedAngle;
            Vertices[Halfedges[heindex].Origin].AdjHalfedge = heindex;

            Faces[Halfedges[heindex].AdjFace].AdjHalfedge = Halfedges[heindex].Next;
            //Faces[Halfedges[heindex].AdjFace].NumVertices--;



            //adjusting halfedges
            Halfedges[Halfedges[heindex].Prev].Valid = false;
            Halfedges[heindex].Prev = Halfedges[Halfedges[heindex].Prev].Prev;
            Halfedges[Halfedges[heindex].Prev].Next = heindex;

            //adjusting twin, if exists
            if (Halfedges[heindex].Twin >= 0) {
                if (Halfedges[Halfedges[heindex].Twin].prescribedAngle < 0.0)
                    Halfedges[Halfedges[heindex].Twin].prescribedAngle = Halfedges[Halfedges[Halfedges[heindex].Twin].Next].prescribedAngle;
                Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Valid = false;
                Halfedges[Halfedges[heindex].Twin].Next = Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Next;
                Halfedges[Halfedges[Halfedges[heindex].Twin].Next].Prev = Halfedges[heindex].Twin;
                Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].AdjHalfedge = Halfedges[Halfedges[heindex].Twin].Next;
                //Faces[Halfedges[Halfedges[heindex].Twin].AdjFace].NumVertices--;
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
        // nextH, prevH, twinH - #H by 1 DCEL traversing operations. twinH(i)=-1 for boundary edges.

        DCEL(const Eigen::VectorXi& D,
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

        }


    };

}


#endif
