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
#include <gmp.h>
#include <gmpxx.h>
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
                             const Eigen::VectorXi& lineData,
                             std::vector<EVector2>& V,
                             DCEL<bool,  void, std::vector<int>, void>& triDcel) {

        V = triangle;

        std::vector<int> inData;  //the lines that are inside
        std::vector<std::pair<EVector2, EVector2>> inSegments;  //parameters of the line segments inside the triangle

        for (int i = 0; i < lines.size(); i++) {
            //checking for intersections with the triangle
            ENumber inParam, outParam;  //the parameters of intersection
            bool intVertex, intEdge, intFace;
            line_triangle_intersection(lines[i], triangle, intEdge, intFace, inParam, outParam);
            //checking cases of intersection
            if ((intEdge < 0) && (intFace < 0))
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
                             DCEL<bool,  void, std::vector<int>, void>& triDcel) {

        //First creating a graph of segment intersection

        //Creating arrangement vertices
        std::vector<EVector2> arrVertices;
        std::vector<std::vector<int>> VS;  //list of participating segments per vertex
        std::vector<std::set<std::pair<ENumber, int>>> SV;  //set of coordinates of intersection per segment
        for (int i=0;i<segments.size();i++) {
            for (int j = i + 1; j < segments.size(); j++) {
                ENumber t1, t2;
                int result = segment_segment_intersection(segments[i].first, segments[i].second - segments[i].first,
                                                  segments[j].first, segments[j].second - segments[j].first, t1, t2));

                if (!result)  //no intersection
                    continue;  //that means the lines intersect away from the triangle.

                if (result==1) {  //pointwise intersection
                    // TODO: figure out what happens if more than two lines at the same spot
                    arrVertices.push_back(segments[i].first * (1 - t1) + segments[i].second * t1);
                    VS[arrVertices.size() - 1].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    SV[i].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 1));
                    SV[i].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 1));
                }

                //Should make this aware of the double and put both into the edge data
                if (result==2) {  //subsegment; now entering two vertices, and letting the edges be entered later
                    arrVertices.push_back(segments[i].first * (ENumber(1) - t1) + segments[i].second * t1);
                    arrVertices.push_back(segments[j].first * (ENumber(1) - t2) + segments[j].second * t2);
                    VS[arrVertices.size() - 2].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    VS[arrVertices.size() - 2].push_back(i);
                    VS[arrVertices.size() - 1].push_back(j);
                    SV[i].insert(std::pair<ENumber, int>(t1, arrVertices.size() - 1));
                    SV[i].insert(std::pair<ENumber, int>(t2, arrVertices.size() - 1));
                }
            }
        }

        //Creating the arrangement edges
        std::vector<std::pair<int, int>> arrEdges;
        std::vector<int> edgeData;
        for (int i=0;i<SV.size();i++){
            for (std::set<std::pair<ENumber, int>>::iterator si = SV[i].begin(); si!=SV[i].end();si++){
                std::set<std::pair<ENumber, int>>::iterator nextsi = si;
                nextsi++;
                if (nextsi!=SV[i].end()) {
                    arrEdges.push_back(std::pair<int, int>(si->second, nextsi->second));
                    edgeData.push_back(data[i]);
                }
            }
        }

        //unifying vertices with the same coordinates
        std::set<std::pair<EVector2, int>, vertexFinder> uniqueVertices;
        std::vector<int> uniqueVertexMap(arrVertices.size());
        std::vector<EVector2> uniqueArrVertices;
        int uniqueCounter=0;
        for (int i=0;i<arrVertices.size();i++){
            std::set<std::pair<EVector2, int>, vertexFinder>::iterator si = uniqueVertices.find(arrVertices[i]);
            if (si==uniqueVertices.end()){
                uniqueVertexMap[i]=uniqueCounter;
                uniqueVertices.insert(std::pair<EVector2, int>(arrVertices[i], uniqueCounter++));
                uniqueArrVertices.push_back(arrVertices[i]);
            } else {
                uniqueVertexMap[i] = si->second;
            }
        }

        arrVertices = uniqueArrVertices;
        V = arrVertices;
        for (int i=0;i<arrEdges.size();i++)
            arrEdges[i]=std::pair<int, int>(uniqueVertexMap[arrEdges[i].first], uniqueVertexMap[arrEdges[i].second]);

        //unifying edges with the same vertices (in case of segment-segment intersections; using the higher data
        Eigen::VectorXi isDeadEdge=Eigen::VectorXi::Constant(arrEdges.size(),0);
        for (int i=0;i<arrEdges.size();i++) {
            for (int j = i + 1; j < arrEdges.size(); j++) {
                if (((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second)) ||
                    ((arrEdges[i].first == arrEdges[j].first) && (arrEdges[i].second == arrEdges[j].second))) {
                    isDeadEdge(j) = 1;
                    edgeData[i] = (edgeData[i] > edgeData[j] ? edgeData[i] : edgeData[j]);
                }
            }
        }
        std::vector<std::pair<int, int>> newArrEdges;
        std::vector<int> newEdgeData;
        for (int i=0;i<arrEdges.size();i++){
            if (isDeadEdge[i])
                continue;
            newArrEdges.push_back(arrEdges[i]);
            newEdgeData.push_back(edgeData[i]);
        }
        arrEdges=newArrEdges;
        edgeData=newEdgeData;

        //TODO:
        //generating the halfedge structure
        //Everything will be twinned at this point
        triDcel.vertices.resize(arrVertices.size());
        triDcel.edges.resize(arrEdges.size());
        triDcel.halfedges.resize(2*arrEdges.size());

        for (int i=0;i<arrVertices.size();i++){
            triDcel.vertices[i].ID = i;
        }

        for (int i=0;i<arrEdges.size();i++) {
            triDcel.edges[i].ID = i;
            triDcel.edges[i].data = edgeData[i];
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

        //orienting segments around each vertex by CCW order
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

            DCEL<bool,  void, std::vector<int>, void>::Face newFace;
            newFace.ID = currFace++;
            int beginHE = i;
            newFace.halfedge=beginHE;
            int currHE = beginHE;
            int counter=0;
            do {
                triDcel.halfedges[currHE]=newFace.ID;
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
            if (sfa<0.0){
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
            triDcel.vertices[triDcel.halfedges[i].vertex].halfedge = triDcel.halfedges[triDcel.halfedges[i].twin].next
            //dcel.VH( dcel.HV( dcel.nextH(i)))= dcel.twinH(i);
            //dcel.VH( dcel.HV(i))= dcel.nextH( dcel.twinH(i));
        }

        //removing dead edges
        triDcel.clean_mesh();
        //EXPENSIVE! comment out after being sure
        triDcel.check_consistency(true, true, true, false);
    }


    void NFunctionMesher::generate_mesh() {

        using namespace std;
        using namespace Eigen;
        //using namespace ::CGAL;

        mesher.clear();

        int numNFunction = mesher.exactNFunctions[0].size();

        //DebugLog.open("Debugging.txt");

        //resolution is set to 10e-6 of bounding box of mesh
        vector<RowVector3d> coordList;
        for (int i = 0; i < origMesh.Vertices.size(); i++)
            coordList.push_back(origMesh.Vertices[i].Coordinates);

        //Bbox_3 boundBox = ::CGAL::bbox_3  ( coordList.begin(), coordList.end());

        /*double minRange = 3276700.0;
         for (int i=0;i<2;i++)
         minRange=std::min(minRange, boundBox.max(i)-boundBox.min(i));*/

        unsigned long Resolution = 1e7; //pow(10,ceil(10/log10(minRange)));
        //cout<<"Resolution: "<<Resolution<<endl;

        for (int findex = 0; findex < origMesh.Faces.size(); findex++) {

            //building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face

            int ebegin = origMesh.Faces[findex].AdjHalfedge;
            int eiterate = ebegin;
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
            vector<ENumber> minFuncs(numNFunction);
            vector<ENumber> maxFuncs(numNFunction);
            for (int k = 0; k < numNFunction; k++) {
                minFuncs[k] = mpq_class(327600);
                maxFuncs[k] = mpq_class(-327600);
            }

            //Arr_2 ParamArr,TriangleArr, FullArr;
            ebegin = origMesh.Faces[findex].AdjHalfedge;
            eiterate = ebegin;
            int currVertex = 0;
            do {
                for (int i = 0; i < numNFunction; i++) {
                    if (origMesh.Halfedges[eiterate].exactNFunction[i] > maxFuncs[i])
                        maxFuncs[i] = origMesh.Halfedges[eiterate].exactNFunction[i];
                    if (origMesh.Halfedges[eiterate].exactNFunction[i] < minFuncs[i])
                        minFuncs[i] = origMesh.Halfedges[eiterate].exactNFunction[i];
                }
                funcValues[currVertex++] = origMesh.Halfedges[eiterate].exactNFunction;
                eiterate = origMesh.Halfedges[eiterate].Next;
            } while (eiterate != ebegin);

            ////////////////////////building the one-triangle arrangement
            ebegin = origMesh.Faces[findex].AdjHalfedge;
            eiterate = ebegin;
            //vector<RowVector2ed> ETriPoints2D;
            //vector<Point2D> TriPoints;
            //vector<RowVector2ed> ETriPoints3D;
            vector<NFunctionMesher::EdgeData> EdgeDatas;
            std::vector<EVector2> ETriPoints2D(3);
            std::vector<EVector3> ETriPoints3D(3);
            ETriPoints2D[0][0] = 0;
            ETriPoints2D[0][1] = 0;
            ETriPoints2D[1][0] = 1;
            ETriPoints2D[1][1] = 0;
            ETriPoints2D[2][0] = 0;
            ETriPoints2D[2][1] = 1;

            /*ETriPoints2D.push_back(RowVector2ed(0,0));
            ETriPoints2D.push_back(RowVector2ed(1,0));
            ETriPoints2D.push_back(RowVector2ed(0,1));*/
            do {
                //cout<<"Halfedges[eiterate].Origin: "<<Halfedges[eiterate].Origin<<endl;
                //Point2D Location((Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis1,(Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis2);
                //cout<<"Location: "<<Location<<endl;*/

                RowVector3d Position = origMesh.Vertices[origMesh.Halfedges[eiterate].Origin].Coordinates;
                //ENumber cx=ENumber((int)(Location.x()*Resolution),Resolution);
                //ENumber cy=ENumber((int)(Location.y()*Resolution),Resolution);
                ENumber x = ENumber((signed long) round((long double) (Position.x()) * Resolution), Resolution);
                ENumber y = ENumber((signed long) round((long double) (Position.y()) * Resolution), Resolution);
                ENumber z = ENumber((signed long) round((long double) (Position.z()) * Resolution), Resolution);
                /*if (abs(x.to_double() - Position.x()) > 10e-7) {
                    cout << "x.to_double(): " << x.to_double() << endl;
                    cout << "Position.x(): " << Position.x() << endl;
                }*/
                //ETriPoints.push_back(EPoint2D(cx,cy));
                //TriPoints.push_back(Location);
                EVector xyz(3);
                xyz[0] = x;
                xyz[1] = y;
                xyz[2] = z;
                ETriPoints3D.push_back(xyz);
                int DomEdge;

                if ((origMesh.Halfedges[eiterate].Twin < 0) || (origMesh.Halfedges[eiterate].Twin > eiterate))
                    DomEdge = eiterate;
                else
                    DomEdge = origMesh.Halfedges[eiterate].Twin;
                NFunctionMesher::EdgeData ed;
                ed.OrigHalfedge = DomEdge;
                ed.isBoundary = (origMesh.Halfedges[eiterate].Twin < 0);
                EdgeDatas.push_back(ed);
                eiterate = origMesh.Halfedges[eiterate].Next;
            } while (ebegin != eiterate);

            for (int i = 0; i < 3; i++) {
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
            }
            * /

            //creating the primal arrangement of lines
            vector<ELine2>
            paramLines;
            vector<EDirection2> isoDirections(numNFunction);
            //int jumps = (numNFunction%2==0 ? 2 : 1);
            for (int funcIter = 0; funcIter < numNFunction/*/jumps*/; funcIter++) {

                vector<EInt> isoValues;
                //cout<<"isoValues: "<<endl;
                EInt q, r;
                CGAL::div_mod(minFuncs[funcIter].numerator(), minFuncs[funcIter].denominator(), q, r);
                EInt minIsoValue = q + (r < 0 ? -1 : 0);
                CGAL::div_mod(maxFuncs[funcIter].numerator(), maxFuncs[funcIter].denominator(), q, r);
                EInt maxIsoValue = q + (r < 0 ? 0 : -1);
                for (EInt isoValue = minIsoValue - 2; isoValue <= maxIsoValue + 2; isoValue++) {
                    //cout<<"isoValue: "<<isoValue<<endl;
                    isoValues.push_back(isoValue);
                }

                //computing gradient of function in plane
                EVector2D e01 = ETriPoints2D[1] - ETriPoints2D[0];
                EVector2D e12 = ETriPoints2D[2] - ETriPoints2D[1];
                EVector2D e20 = ETriPoints2D[0] - ETriPoints2D[2];

                //a and b values of lines
                EVector2D gradVector = funcValues[2][funcIter] * EVector2D(-e01.y(), e01.x()) +
                                       funcValues[0][funcIter] * EVector2D(-e12.y(), e12.x()) +
                                       funcValues[1][funcIter] * EVector2D(-e20.y(), e20.x());

                isoDirections[funcIter] = EDirection2D(gradVector);

                //Number avgFuncValue = (funcValues[0](funcIter)+funcValues[1](funcIter)+funcValues[2](funcIter))/3.0;
                //TODO: find c = z1*u+z2 of ax+by+c(u) ad then use it to generate all values between floor and ceil.

                //pinv of [a 1;b 1;c 1] is [           2*a - b - c,           2*b - a - c,           2*c - b - a]
                //[ b^2 - a*b + c^2 - a*c, a^2 - b*a + c^2 - b*c, a^2 - c*a + b^2 - c*b]/(2*a^2 - 2*a*b - 2*a*c + 2*b^2 - 2*b*c + 2*c^2)

                ENumber a = funcValues[0][funcIter];
                ENumber b = funcValues[1][funcIter];
                ENumber c = funcValues[2][funcIter];
                if ((a == b) && (b == c))
                    continue;  //that means a degenerate function on the triangle

                //cout<<"a,b,c: "<<a.to_double()<<","<<b.to_double()<<","<<c.to_double()<<endl;

                ENumber rhs[3];
                rhs[0] = -gradVector[0] * ETriPoints2D[0].x() - gradVector[1] * ETriPoints2D[0].y();
                rhs[1] = -gradVector[0] * ETriPoints2D[1].x() - gradVector[1] * ETriPoints2D[1].y();
                rhs[2] = -gradVector[0] * ETriPoints2D[2].x() - gradVector[1] * ETriPoints2D[2].y();

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
                    paramLines.push_back(ELine2D(gradVector[0], gradVector[1], currc));
                    //cout<<"paramLine: "<<gradVector[0]<<","<<gradVector[1]<<","<<currc<<endl;
                }
            }

            //cout<<"paramLines.size() :"<<paramLines.size()<<endl;
            CGAL::insert(ParamArr, paramLines.begin(), paramLines.end());

            //giving edge data to curve arrangement
            Arr_2::Edge_iterator eit;
            Arr_2::Originating_curve_iterator ocit;
            for (eit = ParamArr.edges_begin(); eit != ParamArr.edges_end(); ++eit) {
                for (ocit = ParamArr.originating_curves_begin(eit);
                     ocit != ParamArr.originating_curves_end(eit); ++ocit) {
                    EDirection2D thisDirection = EDirection2D(ocit->supporting_line().a(), ocit->supporting_line().b());
                    //cout<<"thisDirection: "<<thisDirection<<endl;
                    for (int paramIter = 0; paramIter < numNFunction/*/jumps*/; paramIter++) {
                        //cout<<"isoDirections[paramIter]: "<<isoDirections[paramIter]<<endl;
                        if ((thisDirection == isoDirections[paramIter]) ||
                            (thisDirection == -isoDirections[paramIter])) {
                            eit->data().funcNum = paramIter;
                            eit->twin()->data().funcNum = paramIter;
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
            overlay(TriangleArr, ParamArr, FullArr, ot);

            /*cout<<"FullArr edges: "<<endl;
             for (eit = FullArr.edges_begin(); eit != FullArr.edges_end(); ++eit)
             cout<<"FullArr eit->data().funcNum: "<<eit->data().funcNum<<endl;*/


            for (Face_iterator fi = FullArr.faces_begin(); fi != FullArr.faces_end(); fi++) {
                if (!fi->data())
                    continue;  //not participating

                Ccb_halfedge_circulator hebegin = fi->outer_ccb();
                Ccb_halfedge_circulator heiterate = hebegin;
                do {

                    if (heiterate->source()->data() < 0) {  //new vertex
                        Vertex NewVertex;
                        NewVertex.ID = funcMesh.Vertices.size();
                        NewVertex.isFunction = (heiterate->source()->data() == -2);
                        funcMesh.Vertices.push_back(NewVertex);
                        heiterate->source()->data() = NewVertex.ID;
                    }

                    if (heiterate->data().ID < 0) {  //new halfedge
                        Halfedge NewHalfedge;
                        NewHalfedge.ID = funcMesh.Halfedges.size();
                        NewHalfedge.isFunction = (heiterate->data().ID == -2);
                        NewHalfedge.Origin = heiterate->source()->data();
                        NewHalfedge.OrigHalfedge = heiterate->data().OrigHalfedge;
                        NewHalfedge.OrigNFunctionIndex = heiterate->data().funcNum;
                        //cout<<"NewHalfedge.OrigParamFunc :"<<NewHalfedge.OrigParamFunc<<endl;
                        funcMesh.Vertices[heiterate->source()->data()].AdjHalfedge = NewHalfedge.ID;
                        funcMesh.Halfedges.push_back(NewHalfedge);
                        heiterate->data().ID = NewHalfedge.ID;
                    }
                    heiterate++;
                } while (heiterate != hebegin);

                //now assigning nexts and prevs
                do {
                    funcMesh.Halfedges[heiterate->data().ID].Next = heiterate->next()->data().ID;
                    funcMesh.Halfedges[heiterate->data().ID].Prev = heiterate->prev()->data().ID;
                    funcMesh.Halfedges[heiterate->data().ID].Twin = heiterate->twin()->data().ID;
                    if (heiterate->twin()->data().ID >= 0)
                        funcMesh.Halfedges[heiterate->twin()->data().ID].Twin = heiterate->data().ID;

                    heiterate++;
                } while (heiterate != hebegin);
            }

            //constructing the actual vertices
            for (Vertex_iterator vi = FullArr.vertices_begin(); vi != FullArr.vertices_end(); vi++) {
                if (vi->data() < 0)
                    continue;

                //finding out barycentric coordinates
                ENumber BaryValues[3];
                ENumber Sum = 0;
                for (int i = 0; i < 3; i++) {
                    ETriangle2D t(vi->point(), ETriPoints2D[(i + 1) % 3], ETriPoints2D[(i + 2) % 3]);
                    BaryValues[i] = t.area();
                    Sum += BaryValues[i];
                }
                for (int i = 0; i < 3; i++)
                    BaryValues[i] /= Sum;

                EPoint3D ENewPosition(0, 0, 0);
                for (int i = 0; i < 3; i++)
                    ENewPosition = ENewPosition + (ETriPoints3D[i] - CGAL::ORIGIN) * BaryValues[i];

                Point3D NewPosition(to_double(ENewPosition.x()), to_double(ENewPosition.y()),
                                    to_double(ENewPosition.z()));
                funcMesh.Vertices[vi->data()].Coordinates = NewPosition;
                funcMesh.Vertices[vi->data()].ECoordinates = ENewPosition;

                //DebugLog<<"Creating Vertex "<<vi->data()<<" with 2D coordinates ("<<vi->point().x()<<","<<vi->point().y()<<") "<<" and 3D Coordinates ("<<std::setprecision(10) <<NewPosition.x()<<","<<NewPosition.y()<<","<<NewPosition.z()<<")\n";
            }

            for (Face_iterator fi = FullArr.faces_begin(); fi != FullArr.faces_end(); fi++) {
                if (!fi->data())
                    continue;

                int FaceSize = 0;
                Ccb_halfedge_circulator hebegin = fi->outer_ccb();
                Ccb_halfedge_circulator heiterate = hebegin;
                do {
                    FaceSize++;
                    heiterate++;
                }
                while (heiterate != hebegin);
                int CurrPlace = 0;

                Face NewFace;
                NewFace.ID = funcMesh.Faces.size();
                //NewFace.NumVertices=FaceSize;
                NewFace.AdjHalfedge = hebegin->data().ID;

                do {
                    //NewFace.Vertices[CurrPlace++]=heiterate->source()->data();
                    funcMesh.Halfedges[heiterate->data().ID].AdjFace = NewFace.ID;
                    heiterate++;
                } while (heiterate != hebegin);
                funcMesh.Faces.push_back(NewFace);
            }

        }

        //devising angles from differences in functions
        //int ratio = (numNFunction%2==0 ? 1 : 2);
        /*for (int hi=0;hi<funcMesh.Halfedges.size();hi++){
          //cout<<"funcMesh.Halfedges[hi].OrigParamFunc: "<<funcMesh.Halfedges[hi].OrigParamFunc<<endl;
          //cout<<"funcMesh.Halfedges[Halfedges[hi].Prev].OrigParamFunc: "<<funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigParamFunc<<endl;
          if ((funcMesh.Halfedges[hi].OrigNFunctionIndex==-1)||(funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigNFunctionIndex==-1))
            funcMesh.Halfedges[hi].prescribedAngle=-1.0;  //one of the edges is a triangle edge, and it will be devised later.
          else{
            //int func1 =(ratio*(funcMesh.Halfedges[hi].OrigParamFunc)) % (numNFunction/(3-ratio));
            //int func2 =(ratio*(funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigParamFunc)) % (numNFunction/(3-ratio));

            double funcOrient1 = funcOrientations(funcMesh.Halfedges[hi].OrigNFunctionIndex);
            double funcOrient2 = funcOrientations(funcMesh.Halfedges[funcMesh.Halfedges[hi].Prev].OrigNFunctionIndex);
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
        }*/
    }
    
} //namespace directional


#endif
