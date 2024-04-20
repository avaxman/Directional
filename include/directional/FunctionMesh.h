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
#include <directional/GMP_definitions.h>
#include <directional/generated_mesh_simplification.h>

namespace directional{

        //arranging a line set on a triangle
        //triangle is represented by a 3x2 matrix of (CCW) coordinates
        //lines are Nx4 matrices of (origin, direction).
        //line data is an integer associated with data on the line that gets inherited to the halfedges
        void arrange_on_triangle(const std::vector<EVector2>& triangle,
                                 const std::vector<std::pair<EVector2, EVector2>>& lines,
                                 const Eigen::VectorXi& lineData,
                                 std::vector<EVector2>& V,
                                 Eigen::MatrixXi& F,
                                 Eigen::MatrixXi& HV,
                                 Eigen::MatrixXi& VH,
                                 Eigen::MatrixXi& FH,
                                 Eigen::MatrixXi& HF,
                                 Eigen::MatrixXi& nextH,
                                 Eigen::MatrixXi& prevH,
                                 Eigen::MatrixXi& twinH,
                                 Eigen::VectorXi& dataH) {
        }



    void GenerateMesh(NFunctionMesher& origMesh, NFunctionMesher& funcMesh){

            using namespace std;
            using namespace Eigen;
            //using namespace ::CGAL;


            funcMesh.Vertices.clear();
            funcMesh.Halfedges.clear();
            funcMesh.Faces.clear();

            int numNFunction=origMesh.Halfedges[0].exactNFunction.size();

            //DebugLog.open("Debugging.txt");

            //resolution is set to 10e-6 of bounding box of mesh
            vector<RowVector3d> coordList;
            for (int i=0;i<origMesh.Vertices.size();i++)
                coordList.push_back(origMesh.Vertices[i].Coordinates);

            //Bbox_3 boundBox = ::CGAL::bbox_3  ( coordList.begin(), coordList.end());

            /*double minRange = 3276700.0;
             for (int i=0;i<2;i++)
             minRange=std::min(minRange, boundBox.max(i)-boundBox.min(i));*/

            unsigned long Resolution=1e7; //pow(10,ceil(10/log10(minRange)));
            //cout<<"Resolution: "<<Resolution<<endl;

            for (int findex=0;findex<origMesh.Faces.size();findex++){

                //building small face overlays of one triangle and a few roughly surrounding hexes to retrieve the structure in the face

                int ebegin=origMesh.Faces[findex].AdjHalfedge;
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
                vector<ENumber> minFuncs(numNFunction);
                vector<ENumber> maxFuncs(numNFunction);
                for (int k=0;k<numNFunction;k++){
                    minFuncs[k] = mpq_class(327600);
                    maxFuncs[k] = mpq_class(-327600);
                }

                //Arr_2 ParamArr,TriangleArr, FullArr;
                ebegin=origMesh.Faces[findex].AdjHalfedge;
                eiterate=ebegin;
                int currVertex=0;
                do{
                    for(int i=0;i<numNFunction;i++){
                        if (origMesh.Halfedges[eiterate].exactNFunction[i]>maxFuncs[i]) maxFuncs[i]=origMesh.Halfedges[eiterate].exactNFunction[i];
                        if (origMesh.Halfedges[eiterate].exactNFunction[i]<minFuncs[i]) minFuncs[i]=origMesh.Halfedges[eiterate].exactNFunction[i];
                    }
                    funcValues[currVertex++]=origMesh.Halfedges[eiterate].exactNFunction;
                    eiterate=origMesh.Halfedges[eiterate].Next;
                }while (eiterate!=ebegin);

                ////////////////////////building the one-triangle arrangement
                ebegin=origMesh.Faces[findex].AdjHalfedge;
                eiterate=ebegin;
                //vector<RowVector2ed> ETriPoints2D;
                //vector<Point2D> TriPoints;
                //vector<RowVector2ed> ETriPoints3D;
                vector<NFunctionMesher::EdgeData> EdgeDatas;
                std::vector<EVector(2)> ETriPoints2D(3);
                std::vector<EVector(3)> ETriPoints3D(3);
                ETriPoints2D[0][0]=0; ETriPoints2D[0][1]=0;
                ETriPoints2D[1][0]=1; ETriPoints2D[1][1]=0;
                ETriPoints2D[2][0]=0; ETriPoints2D[2][1]=1;

                /*ETriPoints2D.push_back(RowVector2ed(0,0));
                ETriPoints2D.push_back(RowVector2ed(1,0));
                ETriPoints2D.push_back(RowVector2ed(0,1));*/
                do{
                    //cout<<"Halfedges[eiterate].Origin: "<<Halfedges[eiterate].Origin<<endl;
                    //Point2D Location((Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis1,(Vertices[Halfedges[eiterate].Origin].Coordinates-Faces[findex].Centroid)*Faces[findex].Basis2);
                    //cout<<"Location: "<<Location<<endl;*/

                    RowVector3d Position=origMesh.Vertices[origMesh.Halfedges[eiterate].Origin].Coordinates;
                    //ENumber cx=ENumber((int)(Location.x()*Resolution),Resolution);
                    //ENumber cy=ENumber((int)(Location.y()*Resolution),Resolution);
                    ENumber x=ENumber((signed long)round((long double)(Position.x())*Resolution),Resolution);
                    ENumber y=ENumber((signed long)round((long double)(Position.y())*Resolution),Resolution);
                    ENumber z=ENumber((signed long)round((long double)(Position.z())*Resolution),Resolution);
                    /*if (abs(x.to_double() - Position.x()) > 10e-7) {
                        cout << "x.to_double(): " << x.to_double() << endl;
                        cout << "Position.x(): " << Position.x() << endl;
                    }*/
                //ETriPoints.push_back(EPoint2D(cx,cy));
                //TriPoints.push_back(Location);
                EVector xyz(3); xyz[0]=x; xyz[1]=y; xyz[2]=z;
                ETriPoints3D.push_back(xyz);
                int DomEdge;

                if ((origMesh.Halfedges[eiterate].Twin<0)||(origMesh.Halfedges[eiterate].Twin>eiterate))
                    DomEdge=eiterate;
                else
                    DomEdge=origMesh.Halfedges[eiterate].Twin;
                NFunctionMesher::EdgeData ed; ed.OrigHalfedge=DomEdge;
                ed.isBoundary=(origMesh.Halfedges[eiterate].Twin<0);
                EdgeDatas.push_back(ed);
                eiterate=origMesh.Halfedges[eiterate].Next;
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
            }*/

             //creating the primal arrangement of lines
             vector<ELine2> paramLines;
            vector<EDirection2> isoDirections(numNFunction);
            //int jumps = (numNFunction%2==0 ? 2 : 1);
            for (int funcIter=0;funcIter<numNFunction/*/jumps*/;funcIter++){

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
                    for (int paramIter = 0;paramIter<numNFunction/*/jumps*/;paramIter++){
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
                        NewHalfedge.OrigNFunctionIndex=heiterate->data().funcNum;
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


        void init(const TriMesh& origMesh,
                  const Eigen::MatrixXd& cutV,
                  const Eigen::MatrixXi& cutF,
                  const Eigen::VectorXd& vertexNFunction,
                  const int N,
                  const Eigen::SparseMatrix<double>& vertexToCornerMat,
                  const Eigen::SparseMatrix<int>& exactVertexToCornerMat,
                  const Eigen::VectorXi& integerVars,
                  const unsigned long resolution=1e7){

            using namespace std;
            using namespace Eigen;
            Vertices.resize(origMesh.V.rows());
            Halfedges.resize(origMesh.HE.rows());
            Faces.resize(origMesh.F.rows());

            //int NFull=(N%2==0 ? N/2: N);

            for (int i=0;i<origMesh.V.rows();i++){
                Vertices[i].Coordinates=origMesh.V.row(i);
                Vertices[i].AdjHalfedge=origMesh.VH(i);
                Vertices[i].ID=i;
            }

            for (int i=0;i<origMesh.HE.rows();i++){
                Halfedges[i].ID=i;
                Halfedges[i].Origin=origMesh.HV(i);
                Halfedges[i].Next=origMesh.nextH(i);
                Halfedges[i].Prev=origMesh.prevH(i);
                Halfedges[i].Twin=origMesh.twinH(i);
                Halfedges[i].AdjFace=origMesh.HF(i);
            }


            for (int i=0;i<origMesh.FH.rows();i++)
                for (int j=0;j<origMesh.FH.cols();j++)


                    for (int i=0;i<origMesh.F.rows();i++){
                        Faces[i].ID=i;
                        Faces[i].AdjHalfedge=origMesh.FH(i);
                    }

            //computing exact rational corner values by quantizing the free variables d and then manually performing the sparse matrix multiplication
            vector<ENumber> exactVertexNFunction(vertexNFunction.size());
            for (int i=0;i<vertexNFunction.size();i++){
                exactVertexNFunction[i]=ENumber((signed long)round((long double)(vertexNFunction(i)*resolution)),(unsigned long)resolution);
                /*if (abs(exactVertexNFunction[i].to_double() - vertexNFunction(i))>10e-8) {
                    cout << "exactVertexNFunction[i].to_double(): " << exactVertexNFunction[i].to_double() << endl;
                    cout << "vertexNFunction(i): " << vertexNFunction(i) << endl;
                    cout << "(long double)(vertexNFunction(i)*resolution): " << (long double)(vertexNFunction(i) * resolution) << endl;
                }*/
            }

            for (int i=0;i<integerVars.size();i++){
                exactVertexNFunction[integerVars(i)]=ENumber((long)round(vertexNFunction(integerVars(i))));
                //cout<<"rounding diff of integer var "<<integerVars(i)<<" is "<<exactVertexNFunction[integerVars(i)].to_double()-vertexNFunction(integerVars(i))<<endl;
            }

            VectorXd cutNFunctionVec = vertexToCornerMat*vertexNFunction;
            vector<ENumber> exactCutNFunctionVec;
            exactSparseMult(exactVertexToCornerMat, exactVertexNFunction,exactCutNFunctionVec);

            //sanity check - comparing exact to double
            double maxError2 = -32767000.0;
            for (int i=0;i<exactCutNFunctionVec.size();i++){
                double fromExact = exactCutNFunctionVec[i].to_double();
                if (abs(fromExact-cutNFunctionVec[i])>maxError2)
                    maxError2 =abs(fromExact-cutNFunctionVec[i]);
            }

            //cout<<"double from exact in halfedges maxError2: "<<maxError2<<endl;

            for (int i=0;i<FH.rows();i++)
                for (int j=0;j<FH.cols();j++){
                    Halfedges[FH(i,j)].exactNFunction.resize(N);
                    Halfedges[FH(i,j)].NFunction = cutNFunctionVec.segment(N*cutF(i,j), N).transpose();
                    for (int k=0;k<N;k++)
                        Halfedges[FH(i,j)].exactNFunction[k] = exactCutNFunctionVec[N*cutF(i,j)+k];
                }

            //sanity check
            double maxError = -32767000.0;
            for (int i=0;i<Halfedges.size();i++){
                for (int j=0;j<N;j++){
                    double fromExact = Halfedges[i].exactNFunction[j].to_double();
                    //cout<<"fromExact: "<<fromExact<<endl;
                    //cout<<"Halfedges[i].NFunction[j]: "<<Halfedges[i].NFunction[j]<<endl;
                    if (abs(fromExact-Halfedges[i].NFunction[j])>maxError)
                        maxError =abs(fromExact-Halfedges[i].NFunction[j]);
                }

            }
            //cout<<"double from exact in halfedges maxError: "<<maxError<<endl;
        }


        //corner angles is per vertex in each F
        void toHedra(Eigen::MatrixXd& generatedV, Eigen::VectorXi& generatedD, Eigen::MatrixXi& generatedF){
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

            /*generatedFfuncNum.resize(Faces.size(),generatedD.maxCoeff());
            cornerAngles=Eigen::MatrixXd::Constant(Faces.size(),generatedD.maxCoeff(),-1.0);
            //prescribedAnglesInt.resize(Faces.size(),generatedD.maxCoeff());
            for (int i=0;i<Faces.size();i++){
              int hebegin = Faces[i].AdjHalfedge;
              int vCount=0;
              int heiterate=hebegin;
              do{
                generatedFfuncNum(i,vCount)=Halfedges[heiterate].OrigNFunctionIndex;
                cornerAngles(i,vCount++)=Halfedges[heiterate].prescribedAngle;
                //prescribedAnglesInt(i,vCount++)=Halfedges[heiterate].prescribedAngleDiff;
                //cout<<"Halfedges[heiterate].prescribedAngleDiff: "<<Halfedges[heiterate].prescribedAngleDiff<<endl;
                heiterate=Halfedges[heiterate].Next;
              }while (heiterate!=hebegin);
            }*/


        }

        NFunctionMesher(){}
        ~NFunctionMesher(){}
    };

} //namespace directional


#endif
