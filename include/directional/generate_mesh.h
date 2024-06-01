// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef FUNCTION_MESH_CLASS_HEADER_FILE
#define FUNCTION_MESH_CLASS_HEADER_FILE

#include <set>
#include <math.h>
#include <vector>
#include <queue>
#include <algorithm>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <directional/DCEL.h>
#include <directional/GMP_definitions.h>
#include <directional/NFunctionMesher.h>

namespace directional{

    //arranging a line set on a triangle
    //triangle is represented by a 3x2 matrix of (CCW) coordinates
    //lines are Nx4 matrices of (origin, direction).
    //line data is an integer associated with data on the line that gets inherited to the halfedges
    //output is the DCEL of the result
    //Outer face is deleted in post-process
    void NFunctionMesher::arrange_on_triangle(const std::vector<EVector2>& triangle,
                             const std::vector<std::pair<EVector2, EVector2>>& lines,
                             const std::vector<int>& lineData,
                             std::vector<EVector2>& V,
                             FunctionDCEL& triDcel) {

        V = triangle;

        std::vector<int> inData;  //the lines that are inside
        std::vector<std::pair<EVector2, EVector2>> inSegments;  //parameters of the line segments inside the triangle

        for (int i = 0; i < lines.size(); i++) {
            //checking for intersections with the triangle
            ENumber inParam, outParam;  //the parameters of intersection
            bool intVertex, intEdge, intFace;
            line_triangle_intersection(lines[i], triangle, intEdge, intFace, inParam, outParam);
            //checking cases of intersection
            if (!intEdge && !intFace)
                continue;   //no (non-measure-zero) intersection

            inData.push_back(lineData[i]);
            EVector2 segSource = lines[i].first + lines[i].second * inParam;
            EVector2 segTarget = lines[i].first + lines[i].second * outParam;
            inSegments.push_back(std::pair<EVector2, EVector2>(segSource, segTarget));
        }

        //pushing in triangle segments
        for (int i = 0; i < 3; i++) {
            inData.push_back(-1);  //no data
            inSegments.push_back(std::pair<EVector2, EVector2>(triangle[i], triangle[(i + 1) % 3]));
        }

        segment_arrangement(inSegments, inData, V, triDcel);
    }



    void NFunctionMesher::segment_arrangement(const std::vector<std::pair<EVector2, EVector2>>& segments,
                             const std::vector<int>& data,
                             std::vector<EVector2>& V,
                             FunctionDCEL& triDcel) {

        //First creating a graph of segment intersection

        //Creating arrangement vertices
        std::vector<EVector2> arrVertices;
        std::vector<std::vector<int>> VS;  //list of participating segments per vertex
        std::vector<std::set<std::pair<ENumber, int>>> SV;  //set of coordinates of intersection per segment
        for (int i=0;i<segments.size();i++) {
            for (int j = i + 1; j < segments.size(); j++) {
                std::vector<std::pair<ENumber, ENumber>> result = segment_segment_intersection(segments[i], segments[j]);

                if (result.empty())  //no intersection
                    continue;  //that means the segments intersect away from the triangle.

                for (int r=0;r<result.size();r++){
                    arrVertices.push_back(segments[i].first * (1 - result[i].first) + segments[i].second * result[i].first);
                    VS[arrVertices.size() - 1].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    SV[i].insert(std::pair<ENumber, int>(result[i].first, arrVertices.size() - 1));
                    SV[j].insert(std::pair<ENumber, int>(result[i].second, arrVertices.size() - 1));

                }
                /*if (result.size) {  //pointwise intersection
                    // TODO: figure out what happens if more than two lines at the same spot
                    arrVertices.push_back(segments[i].first * (1 - t1) + segments[i].second * t1);
                    VS[arrVertices.size() - 1].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    SV[i].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 1));
                    SV[j].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 1));
                }

                //Should make this aware of the double and put both into the edge data
                if (result==2) {  //subsegment; now entering two vertices, and letting the edges be entered later
                    arrVertices.push_back(segments[i].first * (ENumber(1) - t1) + segments[i].second * t1);
                    arrVertices.push_back(segments[j].first * (ENumber(1) - t2) + segments[j].second * t2);
                    VS[arrVertices.size() - 2].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    VS[arrVertices.size() - 2].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    SV[i].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 2));
                    SV[j].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 1));
                    SV[i].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 2));
                    SV[j].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 1));
                }*/
            }
        }

        //Creating the arrangement edges
        std::vector<std::pair<int, int>> arrEdges;
        std::vector<std::vector<int>> edgeData;
        for (int i=0;i<SV.size();i++){
            for (std::set<std::pair<ENumber, int>>::iterator si = SV[i].begin(); si!=SV[i].end();si++){
                std::set<std::pair<ENumber, int>>::iterator nextsi = si;
                nextsi++;
                if (nextsi!=SV[i].end()) {
                    arrEdges.push_back(std::pair<int, int>(si->second, nextsi->second));
                    std::vector<int> newEdgeData(1); newEdgeData[0]=data[i];
                    edgeData.push_back(newEdgeData);
                }
            }
        }

        //unifying vertices with the same coordinates (necessary because some segments may intersect at the same point and segment overlaps
        auto VertexCompare = [](const std::pair<EVector2, int>& a, const std::pair<EVector2, int>& b) { return a.first < b.first; };
        std::set<std::pair<EVector2, int>, decltype(VertexCompare)> uniqueVertices(VertexCompare);
        std::vector<int> uniqueVertexMap(arrVertices.size());
        std::vector<EVector2> uniqueArrVertices;
        int uniqueCounter=0;
        for (int i=0;i<arrVertices.size();i++){
            std::pair<EVector2, int> searchElement(arrVertices[i],-1);
            std::set<std::pair<EVector2, int>, decltype(VertexCompare)>::iterator si = uniqueVertices.find(searchElement);
            if (si==uniqueVertices.end()){
                uniqueVertexMap[i]=uniqueCounter;
                std::pair<EVector2, int> newElement = std::pair<EVector2, int>(arrVertices[i], uniqueCounter++);
                uniqueVertices.insert(newElement);
                uniqueArrVertices.push_back(arrVertices[i]);
            } else {
                uniqueVertexMap[i] = si->second;
            }
        }

        arrVertices = uniqueArrVertices;
        V = arrVertices;
        for (int i=0;i<arrEdges.size();i++)
            arrEdges[i]=std::pair<int, int>(uniqueVertexMap[arrEdges[i].first], uniqueVertexMap[arrEdges[i].second]);

        //unifying edges with the same vertices (aggregating data)
        Eigen::VectorXi isDeadEdge=Eigen::VectorXi::Constant(arrEdges.size(),0);
        for (int i=0;i<arrEdges.size();i++) {
            for (int j = i + 1; j < arrEdges.size(); j++) {
                if (((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second)) ||
                    ((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second))) {
                    isDeadEdge(j) = 1;
                    edgeData[i].insert(edgeData[i].end(), edgeData[j].begin(), edgeData[j].end());
                }
            }
        }
        //cleaning dead edges
        std::vector<std::pair<int, int>> newArrEdges;
        std::vector<std::vector<int>> newEdgeData;
        for (int i=0;i<arrEdges.size();i++){
            if (isDeadEdge[i])
                continue;
            newArrEdges.push_back(arrEdges[i]);
            newEdgeData.push_back(edgeData[i]);
        }
        arrEdges=newArrEdges;
        edgeData=newEdgeData;

        //Generating the DCEL
        triDcel.vertices.resize(arrVertices.size());
        triDcel.edges.resize(arrEdges.size());
        triDcel.halfedges.resize(2*arrEdges.size());

        for (int i=0;i<arrVertices.size();i++){
            triDcel.vertices[i].ID = i;
        }

        for (int i=0;i<arrEdges.size();i++) {
            triDcel.edges[i].ID = i;
            //triDcel.edges[i].data = edgeData[i];   //TODO: halfedges (and edges) don't get any data!
            triDcel.edges[i].halfedge=2*i;

            triDcel.halfedges[2*i].ID=2*i;
            triDcel.halfedges[2*i+1].ID=2*i+1;
            triDcel.halfedges[2*i].vertex=arrEdges[i].first;
            triDcel.halfedges[2*i+1].vertex=arrEdges[i].second;
            triDcel.vertices[arrEdges[i].first].halfedge=2*i;
            triDcel.vertices[arrEdges[i].second].halfedge=2*i+1;
            triDcel.halfedges[2*i].edge=triDcel.halfedges[2*i+1].edge = i;
            triDcel.halfedges[2*i].twin=2*i+1;
            triDcel.halfedges[2*i+1].twin=2*i;

        }

        //Orienting segments around each vertex by CCW order
        for (int i=0;i<arrVertices.size();i++) {
            std::vector<std::pair<int,bool>> adjArrEdges;  //second is direction
            for (int j=0;j<arrEdges.size();j++) {
                if (arrEdges[j].first==i)
                    adjArrEdges.push_back(std::pair<int,bool>(j,true));
                if (arrEdges[j].second==i)
                    adjArrEdges.push_back(std::pair<int,bool>(j,false));
            }//not very efficient but probably not that bad

            std::set<std::pair<ENumber, int>> CCWSegments;
            //using this slope function: https://math.stackexchange.com/questions/1450498/rational-ordering-of-vectors
            for (int j=0;j<adjArrEdges.size();j++) {
                EVector2 edgeVec = arrVertices[arrEdges[adjArrEdges[j].first].second] - arrVertices[arrEdges[adjArrEdges[j].first].first];
                edgeVec = (adjArrEdges[j].second ? edgeVec : -edgeVec);
                ENumber slopeFunc = slope_function(edgeVec);
                CCWSegments.insert(std::pair<ENumber, int>(slopeFunc, j));
            }
            int currHE = -1;
            for (std::set<std::pair<ENumber, int>>::iterator si = CCWSegments.begin(); si!=CCWSegments.end();si++) {
                bool outgoing = adjArrEdges[si->second].second;
                int outCurrHE = (outgoing ? triDcel.edges[adjArrEdges[si->second].second].halfedge  : triDcel.halfedges[triDcel.edges[adjArrEdges[si->second].second].halfedge].twin);
                std::set<std::pair<ENumber, int>>::iterator nextsi = si; nextsi++;
                if (nextsi==CCWSegments.end())
                    nextsi = CCWSegments.begin();

                outgoing = adjArrEdges[nextsi->second].second;
                int outNextHE = (outgoing ? triDcel.edges[adjArrEdges[nextsi->second].second].halfedge  : triDcel.halfedges[triDcel.edges[adjArrEdges[nextsi->second].second].halfedge].twin);
                triDcel.halfedges[triDcel.halfedges[outCurrHE].twin].next=outNextHE;
                triDcel.halfedges[outNextHE].prev = triDcel.halfedges[outCurrHE].twin;
            }
        }

        //generating faces (at this stage, there is also an outer face
        int currFace=0;
        for (int i=0;i<triDcel.halfedges.size();i++) {
            if (triDcel.halfedges[i].face != -1)
                continue;  //already been assigned

            FunctionDCEL::Face newFace;
            newFace.ID = currFace++;
            int beginHE = i;
            newFace.halfedge=beginHE;
            int currHE = beginHE;
            int counter=0;
            do {
                triDcel.halfedges[currHE].face = newFace.ID;
                currHE=triDcel.halfedges[currHE].next;
                counter++;
                assert ("something wrong with the face" && counter<1000);
            }while (currHE!=beginHE);
            triDcel.faces.push_back(newFace);
        }
        int numFaces=currFace;

        //EXPENSIVE! comment out after being sure
        triDcel.check_consistency(true, true, true, false);

        //Removing the outer face and deleting all associated halfedges
        //identifying it by the only polygon with negative signed area (expensive?)
        int outerFace=-1;
        for (int f = 0;f<numFaces;f++){
            std::vector<EVector2> faceVectors;
            int beginHE = triDcel.faces[f].halfedge;
            int currHE = beginHE;
            do{
                faceVectors.push_back(V[triDcel.halfedges[triDcel.halfedges[currHE].next].vertex] - V[triDcel.halfedges[currHE].vertex]);
                currHE= triDcel.halfedges[currHE].next;
            }while(currHE!=beginHE);
            ENumber sfa = signed_face_area(faceVectors);
            if (sfa<ENumber(0)){
                outerFace=f;
                break;
            }
        }
        assert("Didn't find outer face!" && outerFace!=-1);

        //invalidating outer face
        triDcel.faces[outerFace].valid=false;
        for (int i=0;i< triDcel.halfedges.size();i++) {
            if (triDcel.halfedges[i].face != outerFace)
                continue;

            triDcel.halfedges[i].valid=false;
            triDcel.halfedges[triDcel.halfedges[i].twin].twin = -1;
            triDcel.edges[triDcel.halfedges[i].edge].halfedge=triDcel.halfedges[i].twin;
            triDcel.vertices[triDcel.halfedges[i].vertex].halfedge = triDcel.halfedges[triDcel.halfedges[i].twin].next;
            //dcel.VH( dcel.HV( dcel.nextH(i)))= dcel.twinH(i);
            //dcel.VH( dcel.HV(i))= dcel.nextH( dcel.twinH(i));
        }

        //removing dead edges
        triDcel.clean_mesh();
        //EXPENSIVE! comment out after being sure
        triDcel.check_consistency(true, true, true, false);
    }


    void NFunctionMesher::generate_mesh(const unsigned long resolution = 1e7) {

        using namespace std;
        using namespace Eigen;

        //int numNFunction = mesher.exactNFunctions[0].size();

        //DebugLog.open("Debugging.txt");

        //resolution is set to 10e-6 of bounding box of mesh
        /*vector<RowVector3d> coordList;
        for (int i = 0; i < origMesh.V.rows(); i++)
            coordList.push_back(origMesh.Vertices[i].Coordinates);*/

        //Bbox_3 boundBox = ::CGAL::bbox_3  ( coordList.begin(), coordList.end());

        /*double minRange = 3276700.0;
         for (int i=0;i<2;i++)
         minRange=std::min(minRange, boundBox.max(i)-boundBox.min(i));*/

        //pow(10,ceil(10/log10(minRange)));
        //cout<<"Resolution: "<<Resolution<<endl;

        //Looping over all triangles, and arranging parametric lines for each
        for (int findex = 0; findex < origMesh.F.rows(); findex++) {

            //building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face

            //int ebegin = origMesh.Faces[findex].AdjHalfedge;
            //int eiterate = ebegin;
            //vector<Point2D> TriPoints2D;

            //basis for triangle
            /*do{
             Point2D Location((Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis1,(Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis2);
             TriPoints2D.push_back(Location);
             eiterate=Halfedges[eiterate].Next;
             }while (eiterate!=ebegin);
             Triangle2D CurrTri(TriPoints2D[0], TriPoints2D[1], TriPoints2D[2]);*/

            vector<ENumber> triExactNFunction = exactNFunction[findex];

            //vector<vector<ENumber> > funcValues(3);

            //DebugLog<<"Working on triangle "<<findex<<"\n";
            vector<ENumber> minFuncs(mfiData.N);
            vector<ENumber> maxFuncs(mfiData.N);
            for (int k = 0; k < mfiData.N; k++) {
                minFuncs[k] = mpq_class(327600);
                maxFuncs[k] = mpq_class(-327600);
            }

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < mfiData.N; j++) {
                    if (triExactNFunction[3 * i + j] > maxFuncs[j])
                        maxFuncs[j] = triExactNFunction[3 * i + j];
                    if (triExactNFunction[3 * i + j] < minFuncs[j])
                        minFuncs[j] = triExactNFunction[3 * i + j];
                }
            }

            ////////////////////////building the one-triangle arrangement
            //ebegin = origMesh.Faces[findex].AdjHalfedge;
            //eiterate = ebegin;
            //vector<RowVector2ed> ETriPoints2D;
            //vector<Point2D> TriPoints;
            //vector<RowVector2ed> ETriPoints3D;
            //vector<NFunctionMesher::EdgeData> EdgeDatas;
            std::vector<EVector2> ETriPoints2D(3);
            std::vector<EVector3> ETriPoints3D(3);
            ETriPoints2D[0][0] = 0;
            ETriPoints2D[0][1] = 0;
            ETriPoints2D[1][0] = 1;
            ETriPoints2D[1][1] = 0;
            ETriPoints2D[2][0] = 0;
            ETriPoints2D[2][1] = 1;

            for (int i = 0; i < 3; i++) {
                RowVector3d position = origMesh.V.row(origMesh.F(findex, i));
                //ENumber cx=ENumber((int)(Location.x()*Resolution),Resolution);
                //ENumber cy=ENumber((int)(Location.y()*Resolution),Resolution);
                ENumber x = ENumber((signed long) round((long double) (position(0)) * resolution), resolution);
                ENumber y = ENumber((signed long) round((long double) (position(1)) * resolution), resolution);
                ENumber z = ENumber((signed long) round((long double) (position(2)) * resolution), resolution);

                EVector3 xyz;
                xyz[0] = x;
                xyz[1] = y;
                xyz[2] = z;
                ETriPoints3D.push_back(xyz);
            }

            /*ETriPoints2D.push_back(RowVector2ed(0,0));
            ETriPoints2D.push_back(RowVector2ed(1,0));
            ETriPoints2D.push_back(RowVector2ed(0,1));*/
            //do {
            //cout<<"Halfedges[eiterate].Origin: "<<Halfedges[eiterate].Origin<<endl;
            //Point2D Location((Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis1,(Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis2);
            //cout<<"Location: "<<Location<<endl;*/


            /*if (abs(x.to_double() - Position.x()) > 10e-7) {
                cout << "x.to_double(): " << x.to_double() << endl;
                cout << "Position.x(): " << Position.x() << endl;
            }*/
            //ETriPoints.push_back(EPoint2D(cx,cy));
            //TriPoints.push_back(Location);

            /*int DomEdge;

            if ((origMesh.Halfedges[eiterate].Twin < 0) || (origMesh.Halfedges[eiterate].Twin > eiterate))
                DomEdge = eiterate;
            else
                DomEdge = origMesh.Halfedges[eiterate].Twin;
            NFunctionMesher::EdgeData ed;
            ed.OrigHalfedge = DomEdge;
            ed.isBoundary = (origMesh.Halfedges[eiterate].Twin < 0);
            EdgeDatas.push_back(ed);
            eiterate = origMesh.Halfedges[eiterate].Next;*/
            //} while (ebegin != eiterate);

            /*for (int i = 0; i < 3; i++) {
                X_monotone_curve_2 c = ESegment2D(ETriPoints2D[i], ETriPoints2D[(i + 1) % 3]);
                Halfedge_handle he = CGAL::insert_non_intersecting_curve(TriangleArr, c);
                he->set_data(EdgeDatas[i]);
                if (EdgeDatas[i].isBoundary)
                    he->source()->data() = he->target()->data() = 0;
                else
                    he->source()->data() = he->target()->data() = 1;

                he->twin()->set_data(EdgeDatas[i]);
            }

            for (Face_iterator fi = TriangleArr.faces_begin(); fi != TriangleArr.faces_end(); fi++) {
                if (fi->is_unbounded())
                    fi->data() = 0;
                else
                    fi->data() = 1;
            }*/

            //Generating the parametric lines for the canonical triangle
            vector<std::pair<EVector2, EVector2>> paramLines;
            vector<int> lineData;
            //VectorXi lineData;
            //vector<ELine2> paramLines;
            vector<EVector2> isoDirections(mfiData.N);
            //int jumps = (numNFunction%2==0 ? 2 : 1);
            for (int funcIter = 0; funcIter < mfiData.N/*/jumps*/; funcIter++) {

                vector<EInt> isoValues;
                //cout<<"isoValues: "<<endl;
                EInt q, r;
                div_mod(minFuncs[funcIter].get_num_mpz_t(), minFuncs[funcIter].get_den_mpz_t(), q, r);
                EInt minIsoValue = q + (r < 0 ? -1 : 0);
                div_mod(maxFuncs[funcIter].get_num_mpz_t(), maxFuncs[funcIter].get_den_mpz_t(), q, r);
                EInt maxIsoValue = q + (r < 0 ? 0 : -1);
                for (EInt isoValue = minIsoValue - 2; isoValue <= maxIsoValue + 2; isoValue++) {
                    //cout<<"isoValue: "<<isoValue<<endl;
                    isoValues.push_back(isoValue);
                }

                //computing gradient of function in plane
                EVector2 e01 = ETriPoints2D[1] - ETriPoints2D[0];
                EVector2 e12 = ETriPoints2D[2] - ETriPoints2D[1];
                EVector2 e20 = ETriPoints2D[0] - ETriPoints2D[2];

                //a and b values of lines
                EVector2 gradVector = triExactNFunction[2 * mfiData.N + funcIter] * EVector2({-e01[1], e01[0]}) +
                                      triExactNFunction[0 * mfiData.N + funcIter] * EVector2({-e12[1], e12[0]}) +
                                      triExactNFunction[1 * mfiData.N + funcIter] * EVector2({-e20[1], e20[0]});

                isoDirections[funcIter] = gradVector;

                //Number avgFuncValue = (funcValues[0](funcIter)+funcValues[1](funcIter)+funcValues[2](funcIter))/3.0;
                //TODO: find c = z1*u+z2 of ax+by+c(u) ad then use it to generate all values between floor and ceil.

                //pinv of [a 1;b 1;c 1] is [           2*a - b - c,           2*b - a - c,           2*c - b - a]
                //[ b^2 - a*b + c^2 - a*c, a^2 - b*a + c^2 - b*c, a^2 - c*a + b^2 - c*b]/(2*a^2 - 2*a*b - 2*a*c + 2*b^2 - 2*b*c + 2*c^2)

                ENumber a = triExactNFunction[0 * mfiData.N + funcIter];
                ENumber b = triExactNFunction[1 * mfiData.N + funcIter];
                ENumber c = triExactNFunction[2 * mfiData.N + funcIter];
                if ((a == b) && (b == c))
                    continue;  //that means a degenerate function on the triangle

                //cout<<"a,b,c: "<<a.to_double()<<","<<b.to_double()<<","<<c.to_double()<<endl;

                ENumber rhs[3];
                rhs[0] = -gradVector[0] * ETriPoints2D[0][0] - gradVector[1] * ETriPoints2D[0][1];
                rhs[1] = -gradVector[0] * ETriPoints2D[1][0] - gradVector[1] * ETriPoints2D[1][1];
                rhs[2] = -gradVector[0] * ETriPoints2D[2][0] - gradVector[1] * ETriPoints2D[2][1];

                ENumber invM[2][3];
                invM[0][0] = 2 * a - b - c;
                invM[0][1] = 2 * b - a - c;
                invM[0][2] = 2 * c - b - a;
                invM[1][0] = b * b - a * b + c * c - a * c;
                invM[1][1] = a * a - b * a + c * c - b * c;
                invM[1][2] = a * a - c * a + b * b - c * b;
                for (int row = 0; row < 2; row++)
                    for (int col = 0; col < 3; col++)
                        invM[row][col] /= (ENumber(2) * (a * a - a * b - a * c + b * b - b * c + c * c));

                //cout<<(ENumber(2)*(a*a - a*b - a*c + b*b- b*c + c*c)).to_double()<<endl;

                ENumber x[2];
                x[0] = invM[0][0] * rhs[0] + invM[0][1] * rhs[1] + invM[0][2] * rhs[2];
                x[1] = invM[1][0] * rhs[0] + invM[1][1] * rhs[1] + invM[1][2] * rhs[2];


                //RowVectorXd x = lhs.colPivHouseholderQr().solve(rhs).transpose();

                //sanity check
                ENumber error[3];
                error[0] = x[0] * a + x[1] - rhs[0];
                error[1] = x[0] * b + x[1] - rhs[1];
                error[2] = x[0] * c + x[1] - rhs[2];


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
                for (int isoIndex = 0; isoIndex < isoValues.size(); isoIndex++) {
                    //ENumber isoVec[2];
                    //isoVec[0]=isoValues[isoIndex];
                    //isoVec[1]= ENumber(1);
                    ENumber currc = isoValues[isoIndex] * x[0] + x[1];
                    // ENumber a=ENumber((int)(gradVector[0]*Resolution),Resolution);
                    //ENumber b=ENumber((int)(gradVector[1]*Resolution),Resolution);
                    //ENumber c=ENumber((int)(currc(0)*Resolution),Resolution);
                    EVector2 lineVector;
                    lineVector[0] = -gradVector[1];
                    lineVector[1] = gradVector[0];
                    EVector2 linePoint;
                    if (gradVector[1] != ENumber(0)) {
                        linePoint[0] = ENumber(0);
                        linePoint[1] = -currc / gradVector[1];
                    } else {
                        linePoint[1] = ENumber(0);
                        linePoint[0] = -currc / gradVector[0];
                    }

                    paramLines.push_back(std::pair<EVector2, EVector2>(linePoint, lineVector));
                    lineData.push_back(funcIter);
                    //cout<<"paramLine: "<<gradVector[0]<<","<<gradVector[1]<<","<<currc<<endl;
                }
            }

            FunctionDCEL localArrDcel;
            vector<EVector2> localV;
            arrange_on_triangle(ETriPoints2D, paramLines, lineData, localV, localArrDcel);

            //vector<EVector3> ELocalV3D(localV.size());
           // MatrixXd localV3D(localV.size(), 3);
            //converting the vertices to 3D
            for (int i = 0; i < localV.size(); i++) {
                //checking if this is an original vertex
                bool isOrigTriangle = false;
                for (int j = 0; j < 3; j++) {
                    if (localV[i] == ETriPoints2D[j]) {
                        localArrDcel.vertices[i].data.eCoords = ETriPoints3D[j];
                        isOrigTriangle = true;
                    }
                }

                if (isOrigTriangle)
                    continue;    //this is probably not needed but covered by barycentric coordinates but w/e

                //finding out barycentric coordinates
                ENumber baryValues[3];
                ENumber sum = 0;
                for (int j = 0; j < 3; j++) {
                    //ETriangle2D t(vi->point(), ETriPoints2D[(i + 1) % 3], ETriPoints2D[(i + 2) % 3]);
                    vector<EVector2> inTri(3);
                    inTri[0] = localV[i];
                    inTri[1] = ETriPoints2D[(j + 1) % 3];
                    inTri[2] = ETriPoints2D[(j + 2) % 3];

                    baryValues[j] = triangle_area(inTri);
                    sum += baryValues[j];
                }
                for (int j = 0; j < 3; j++)
                    baryValues[j] /= sum;

                localArrDcel.vertices[i].data.eCoords = EVector3({0, 0, 0});
                for (int j = 0; j < 3; j++)
                    localArrDcel.vertices[i].data.eCoords = localArrDcel.vertices[i].data.eCoords + ETriPoints3D[i] * baryValues[i];

                localArrDcel.vertices[i].data.coords << localArrDcel.vertices[i].data.eCoords[0].get_d(),
                        localArrDcel.vertices[i].data.eCoords[1].get_d(),
                        localArrDcel.vertices[i].data.eCoords[2].get_d();

            }

            //aggregating to the general DCEL
            genDcel.aggregage_dcel(localArrDcel);
        }
    }
    
} //namespace directional


#endif
