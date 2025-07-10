// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_EXACT_GEOMETRIC_DEFINITIONS_H
#define DIRECTIONAL_EXACT_GEOMETRIC_DEFINITIONS_H

#include <set>
#include <math.h>
#include <vector>
#include <queue>
#include <cassert>
#include <algorithm>
#include <utility>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#ifdef USE_GMP_ENABLED
#include <directional/ENumber_GMP.h>
#else
#include <directional/ENumber_internal.h>
#endif


//This header file concentrates geometric operations on vectors, segments, lines, and arrangement in exact rational numbers.
namespace directional{
    
    template <size_t Size>
    class EVector {
    public:
        EVector() : data(Size) {}
        
        // Other methods can be added as needed
        const ENumber& operator[](const size_t index) const {
            return data[index];
        }
        
        ENumber& operator[](const size_t index) {
            return data[index];
        }
        
        EVector<Size> operator+(const EVector<Size>& ev) const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=data[i]+ev.data[i];
            
            return newVec;
        }
        EVector<Size> operator-() const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=-data[i];
            
            return newVec;
        }
        EVector<Size> operator-(const EVector<Size>& ev) const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=data[i]-ev.data[i];
            
            return newVec;
        }
        
        EVector<Size> operator*(const ENumber s) const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=data[i]*s;
            
            return newVec;
        }
        
        bool operator==(const EVector<Size>& ev) const{
            bool equal=true;
            for (int i=0;i<Size;i++)
                equal = equal & (ev.data[i]==data[i]);
            return equal;
        }
        
        //for the sake of sorting
        bool operator<(const EVector<Size>& ev) const{
            for (int i=0;i<Size;i++)
                if (data[i]!=ev.data[i])
                    return data[i]<ev.data[i];
            return false;
        }
        
        EVector(const std::initializer_list<ENumber>& args){
            data.insert(data.end(), args.begin(), args.end());
        }
        
        EVector operator=(const EVector<Size>& evec){
            data=evec.data;
            return *this;
        }
        
        Eigen::RowVectorXd to_double() const{
            Eigen::RowVectorXd doubleVec(Size);
            for (int i=0;i<Size;i++)
                doubleVec(i)=data[i].to_double();
            return doubleVec;
        }
        
        ENumber cross(const EVector<Size>& evec) const{
            static_assert("This method only works for Size==2" && Size==2);
            return data[0]*evec.data[1] - data[1]*evec.data[0];
        }
        
        ENumber max_abs() const{
            ENumber maxAbs(-1);
            for (int i=0;i<Size;i++)
                if (data[i].abs()>maxAbs) maxAbs = data[i].abs();
            return maxAbs;
        }
        
        ENumber operator*(const EVector<Size>& evec) const{
            ENumber dotProd(0);
            for (int i=0;i<Size;i++)
                dotProd+=data[i]*evec.data[i];
            return dotProd;
        }
        
        template<size_t _Size>
        friend std::ostream& operator<<(std::ostream& os, const EVector<_Size>& evec);
        
        /*void canonicalize(){
         for (int i=0;i<Size;i++)
         data[i].canonicalize();
         }*/
        
        //protected:
        std::vector<ENumber> data;
    };
    
    template<size_t Size>
    EVector<Size> operator*(ENumber scalar, const EVector<Size>& vec) {
        return vec * scalar; // Leverage the previous operator*
    }
    
    template<size_t Size>
    std::ostream& operator<<(std::ostream& os, const EVector<Size>& evec) {
        os << "(";
        for (int i=0;i<Size-1;i++)
            os<< evec[i].to_double()<<",";
        os<<evec[Size-1].to_double()<<")";
        return os;
    }
    
    
    
    typedef EVector<2> EVector2;
    typedef EVector<3> EVector3;
    
    struct Segment2 {
    public:
        EVector2 source, target;
        Segment2(const EVector2& _source, const EVector2& _target){
            source=_source; target=_target;
        }
        
        Segment2 operator=(const Segment2& seg2){
            source = seg2.source;
            target = seg2.target;
            return *this;
        }
        
        Segment2(){}
        
        friend std::ostream& operator<<(std::ostream& os, const Segment2& seg);
    };
    
    std::ostream& operator<<(std::ostream& os, const Segment2& seg) {
        os << "Segment2(" << seg.source << "->" << seg.target << ")";
        return os;
    }
    
    struct Line2{
    public:
        EVector2 point, direction;
        Line2(const EVector2& _point, const EVector2& _direction){
            point=_point; direction=_direction;
        }
        
        friend std::ostream& operator<<(std::ostream& os, const Line2& seg);
        
        ENumber point_param(const EVector2& p) const {  //if the point is not on the line, this is the parameter of the orthogonally-projected point
            return (direction*(p-point))/(direction*direction);
        }
    };
    
    std::ostream& operator<<(std::ostream& os, const Line2& line) {
        os << "Line2(" << line.point << " + " << line.direction << ")";
        return os;
    }
    
    
    struct LinePencil{
        int numLines;
        EVector2 direction;   //the mutual direction along the line
        EVector2 p0, pVec;  //p0 is the origin of the first line. pVec is the vector between the origins of the lines (p0(I+1)-p0(I) = pvec)
        
        inline Line2 line(const int lineNum) const{
            return Line2(p0 + pVec*ENumber(lineNum), direction);
        }
    };
    
    
    ENumber squaredDistance(const EVector3& v1, const EVector3& v2){
        ENumber sd(0);
        //assert("vectors are of different dimensions!" && v1.size()==v2.size());
        for (int i=0;i<3;i++)
            sd+=(v1[i]-v2[i])*(v1[i]-v2[i]);  //maybe it's not efficient
        
        return sd;
    }
    
    //produces y = M*x
    void exactSparseMult(const Eigen::SparseMatrix<int> M, const std::vector<ENumber>& x,std::vector<ENumber>& y){
        y.resize(M.rows());
        
        for (int i=0;i<y.size();i++)
            y[i]=ENumber(0);
        
        for (int k=0; k<M.outerSize(); ++k)
            for (Eigen::SparseMatrix<int>::InnerIterator it(M,k); it; ++it){
                /*if (it.row()==219){
                 std::cout<<"it.value(): "<<it.value()<<std::endl;
                 std::cout<<"x[it.col()]: "<<x[it.col()].to_double()<<std::endl;
                 std::cout<<"ENumber((long)it.value())*x[it.col()]: "<<(ENumber((long)it.value())*x[it.col()]).to_double()<<std::endl;
                 std::cout<<"y[it.row()]: "<<y[it.row()].to_double()<<std::endl;
                 }*/
                y[it.row()]+=ENumber((long)it.value())*x[it.col()];
                
            }
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
    
    //This assumes components is already resized to the correct |v|
    //not very efficient but probably not terrible
    int connectedComponents(const std::vector<std::pair<int,int>>& matches, std::vector<int>& components)
    {
        for (int i=0;i<components.size();i++)
            components[i]=-1;
        
        std::vector<std::vector<int>> VV(components.size());
        for (int i=0;i<matches.size();i++){
            VV[matches[i].first].push_back(matches[i].second);
            VV[matches[i].second].push_back(matches[i].first);
        }
        
        std::deque<int> nextVertexQueue;
        for (int i=0;i<components.size();i++)
            nextVertexQueue.push_front(i);
        
        int numComponents=0;
        while (!nextVertexQueue.empty()) {
            int nextVertex=nextVertexQueue.front();
            nextVertexQueue.pop_front();
            if (components[nextVertex]==-1) {  //first components
                //std::cout << "New component " << numComponents << " seed vertex " << nextVertex << std::endl;
                components[nextVertex] = numComponents++;
            }
            
            //Otherwise, doing DFS on edges
            for (int i=0;i<VV[nextVertex].size();i++){
                if (components[VV[nextVertex][i]]==-1){
                    components[VV[nextVertex][i]]=components[nextVertex];
                    //std::cout<<"adding vertex "<<VV[nextVertex][i]<<" to component "<<components[nextVertex]<<std::endl;
                    nextVertexQueue.push_front(VV[nextVertex][i]);
                } else assert(components[VV[nextVertex][i]]==components[nextVertex]);
            }
        }
        return numComponents;
    }
    
    int line_line_intersection(const Line2& line1,
                               const Line2& line2,
                               ENumber& t1,
                               ENumber& t2){
        //std::cout<<"Computing intersection between line: "<<line1<<" and line "<<line2<<std::endl;
        ENumber v1v2 = line1.direction.cross(line2.direction);
        if (v1v2==ENumber(0)){
            EVector2 pointVec = line1.point-line2.point;
            //std::cout<<"lines are parallel"<<std::endl;
            return  (pointVec.cross(line1.direction)==ENumber(0) ? 2 : 0);
        }
        
        //std::cout<<"exact line2.point[0]-line1.point[0]: "<<(line2.point[0]-line1.point[0]).to_double()<<std::endl;
        //std::cout<<"double line2.point[0]-line1.point[0]: "<<(line2.point[0].to_double()-line1.point[0].to_double())<<std::endl;
        
        //std::cout<<"exact (line2.point[0]-line1.point[0])*(line2.direction[1]): "<<((line2.point[0]-line1.point[0])*(line2.direction[1])).to_double()<<std::endl;
        //std::cout<<"double (line2.point[0]-line1.point[0])*(line2.direction[1]):  "<<(line2.point[0].to_double()-line1.point[0].to_double())*(line2.direction[1].to_double())<<std::endl;
        EVector2 p12 = line2.point-line1.point;
        t1 = p12.cross(line2.direction)/v1v2;
        t2 = p12.cross(line1.direction)/v1v2;
        //std::cout<<"t1, t2: "<<t1.to_double()<<","<<t2.to_double()<<std::endl;
        //std::cout<<"line1.point+t1*line1.direction: "<<line1.point+t1*line1.direction<<std::endl;
        //std::cout<<"line2.point+t2*line2.direction: "<<line2.point+t2*line2.direction<<std::endl;
        //EVector2 diff = (line1.point+t1*line1.direction) - (line2.point+t2*line2.direction);
        //diff.canonicalize();
        //assert("line_line_intersection is wrong!" && diff == EVector2());
        //std::cout<<"lines intersect at "<<(line1.point+t1*line1.direction)<<std::endl;
        //std::cout<<"parameters: t1:"<<t1.get_d()<<", t2: "<<t2.get_d()<<std::endl;
        return 1;
        
    }
    
    
    //returns a generator for the grid of intersections, parameterized by p00 + pVec1*isoValue1 + pVec2*isoValue2,
    //txp00 is the t(1 or 2) of the p00 point in each respective line
    //dtx (1 or 2) is the dt going along each pVecx (1 or 2)
    //for iso1 and iso2 in the respective line pencil ranges
    //result = 2 is only acceptable if |lp2| = 1, not handling parallel full line pencils (shouldn't be unless the parameterization is degenerate).
    inline int linepencil_intersection(const LinePencil& lp1,
                                       const LinePencil& lp2,
                                       Eigen::Matrix<ENumber, 2, 1>& t00,
                                       Eigen::Matrix<ENumber, 2, 2>& I2dt,    //t of lines is I2dt*[I1;I2]+t00
                                       EInt& iso1Overlap){   //overlap line in lp1 matching that of lp2. ignored unless result==2
        
        
        //std::cout<<"Line pencil intersection between "<<std::endl;
        //std::cout<<"lp1 numLines: "<<lp1.numLines<<" p0: "<<lp1.p0<<" direction: "<<lp1.direction<<" pVec: "<<lp1.pVec<<std::endl;
        //std::cout<<"lp2 numLines: "<<lp2.numLines<<" p0: "<<lp2.p0<<" direction: "<<lp2.direction<<" pVec: "<<lp2.pVec<<std::endl;
        
        ENumber v1v2 = lp1.direction.cross(lp2.direction);
        if (v1v2==ENumber(0)){  //the two pencils are parallel; looking which line overlaps lp2
            assert("linepencil_intersection: parallel line pencils must have |lp2|=1; parameterization must be degenerate! " && lp2.numLines==1);
            EVector2 p012vec = lp2.p0-lp1.p0;
            ENumber pVecv1 = lp1.pVec.cross(lp1.direction);
            assert("pVecv1 shouldn't be zero! " && pVecv1!=ENumber(0));
            ENumber overlapIso = p012vec.cross(lp1.direction) / pVecv1;
            if (overlapIso.den()==1){
                iso1Overlap = overlapIso.num().convert();
                return ((iso1Overlap >= 0 && iso1Overlap <= lp1.numLines) ? 2 : 0);  //only an overlap if integer and within range of lp1
            } else return 0;  //not an integer; no overlap for certain
        }
        
        //line pencils are next not overlapping; finding grid parameters
        
        //p00 is the intersection of both min isovalues
        EVector2 p012 = lp2.p0-lp1.p0;
        t00<<p012.cross(lp2.direction) / v1v2, p012.cross(lp1.direction) / v1v2;
        //EVector2 diff = (lp1.p0 + t00(0) * lp1.direction) - (lp2.p0 + t00(1) * lp2.direction);
        //assert("line_pencil original point intersection is wrong!" && diff == EVector2());
        //computing the dts
        I2dt<<-lp1.pVec.cross(lp2.direction), lp2.pVec.cross(lp2.direction),
        -lp1.pVec.cross(lp1.direction), lp2.pVec.cross(lp1.direction);
        I2dt.array()/=v1v2;
        //tests (expensive):
        //Eigen::Matrix<ENumber, 2,1> ITest; ITest<<ENumber(3),ENumber(2);
        //Eigen::Matrix<ENumber, 2,1> gridPointtValues = t00 + I2dt*ITest;
        //EVector2 pLine1 = (lp1.p0 + lp1.pVec*(ITest(0))+ gridPointtValues(0) * lp1.direction);
        //EVector2 pLine2 = (lp2.p0 + lp2.pVec*(ITest(1))+ gridPointtValues(1) * lp2.direction);
        //diff = pLine1 - pLine2;
        //assert("line_pencil dts computation is wrong!" && diff == EVector2());
        
        //std::cout<<"t00: "<<t00(0).to_double()<<", "<<t00(1).to_double()<<std::endl;
        //std::cout<<"I2dt: "<<I2dt(0,0).to_double()<<", "<<I2dt(0,1).to_double()<<", "<<I2dt(1,0).to_double()<<", "<<I2dt(1,1).to_double()<<std::endl;
        //std::cout<<"iso1Overlap "<<iso1Overlap<<std::endl;
        return 1;
        
    }
    
    std::vector<std::pair<ENumber, ENumber>> segment_segment_intersection(const Segment2& seg1,
                                                                          const Segment2& seg2){
        
        ENumber t1, t2;
        //std::cout<<"Computing intersection of "<<seg1<<" and "<<seg2<<std::endl;
        int result = line_line_intersection(Line2(seg1.source, seg1.target-seg1.source),
                                            Line2(seg2.source, seg2.target-seg2.source),t1, t2);
        
        if (result==0) {
            //std::cout<<"supporting lines don't intersect"<<std::endl;
            return std::vector<std::pair<ENumber, ENumber>>(); //no intersection
        }
        
        if (result==1) {  //a single intersection at most; should check t1 and t2
            //std::cout<<"single intersection"<<std::endl;
            if ((t1>=ENumber(0))&&(t1<=ENumber(1))&&(t2>=ENumber(0))&&(t2<=ENumber(1))){
                std::vector<std::pair<ENumber, ENumber>> point(1);
                point[0]=std::pair<ENumber, ENumber>(t1,t2);
                //std::cout<<"Intersecting at parameters "<<point[0].first.get_d()<<","<<point[0].second.get_d()<<std::endl;
                return point;
            } else{
                //std::cout<<"Intersecting out of parameter bounds"<<std::endl;
                return std::vector<std::pair<ENumber, ENumber>>(); //no intersection
            }
        }
        
        if (result==2){  //lines overlap; should check the segments overlap and then return both overlap points (order not important)
            //std::cout<<"Supporting lines overlap"<<std::endl;
            EVector2 vec = seg1.target-seg1.source;
            int axis = (vec[0]!=ENumber(0) ? 0 : 1);
            Segment2 sortSeg1, sortSeg2;
            if (seg1.source[axis]<seg1.target[axis]) sortSeg1=seg1; else sortSeg1=Segment2(seg1.target,seg1.source);
            if (seg2.source[axis]<seg2.target[axis]) sortSeg2=seg2; else sortSeg2=Segment2(seg2.target,seg2.source);
            EVector2 startPoint = (sortSeg1.source[axis]>sortSeg2.source[axis] ? sortSeg1.source : sortSeg2.source);
            EVector2 endPoint = (sortSeg1.target[axis]>sortSeg2.target[axis] ? sortSeg2.target : sortSeg1.target);
            
            //TODO: what is the policy with segment-vertex tangency?
            if (startPoint[axis]<endPoint[axis]) { //there is a (non-zero) intersection
                ENumber startAtSeg1 = (startPoint[axis]-seg1.source[axis])/(seg1.target[axis]-seg1.source[axis]);
                ENumber startAtSeg2 = (startPoint[axis]-seg2.source[axis])/(seg2.target[axis]-seg2.source[axis]);
                ENumber endAtSeg1 = (endPoint[axis]-seg1.source[axis])/(seg1.target[axis]-seg1.source[axis]);
                ENumber endAtSeg2 = (endPoint[axis]-seg2.source[axis])/(seg2.target[axis]-seg2.source[axis]);
                std::vector<std::pair<ENumber, ENumber>> points(2);
                points[0] = std::pair<ENumber, ENumber>(startAtSeg1, startAtSeg2);
                points[1] = std::pair<ENumber, ENumber>(endAtSeg1, endAtSeg2);
                //std::cout<<"Intersecting at parameters "<<points[0].first.to_double()<<","<<points[0].second.to_double()<<" and "<<points[1].first.to_double()<<","<<points[1].second.to_double()<<std::endl;
                return points;
            } else{
                //std::cout<<"parameters don't overlap"<<std::endl;
                return std::vector<std::pair<ENumber, ENumber>>(); //no intersection
            }
            
        }
        
        
    }
    
    
    
    std::vector<ENumber> line_segment_intersection(const Line2& line,
                                                   const Segment2& segment){
        
        //std::cout<<"Computing intersection between line :"<<line<<" and segment "<<segment<<std::endl;
        Line2 segLine(segment.source, segment.target-segment.source);
        ENumber t1, t2;
        int intersectType=line_line_intersection(line, segLine, t1, t2);
        if (intersectType==0) {  //no intersection
            //std::cout<<"No intersection (lines parallel)"<<std::endl;
            return std::vector<ENumber>();
        } else if (intersectType==2) { //the entire segment is contained in the line
            std::vector<ENumber> result(2); result[0]=ENumber(0); result[1]=ENumber(1);   //This is wrong! should have been the line parameters
            std::cout<<"Entire segment is contained in line"<<std::endl;
            return result;
        } else { //(intersectType==1)
            if ((t2>=ENumber(0)) && (t2<=ENumber(1))){
                std::vector<ENumber> result(1); result[0]=t1;
                //std::cout<<"Intersecting at line parameter t1:"<<t1.get_d()<<std::endl;
                return result;
            } else {
                //std::cout<<"No intersection  (out of segment parameter)"<<std::endl;
                return std::vector<ENumber>();
            }
        }
    }
    
    void line_triangle_intersection(const Line2& line,
                                    const std::vector<EVector2> triangle,
                                    bool& intEdge,
                                    bool& intFace,
                                    ENumber& inParam,
                                    ENumber& outParam){
        
        inParam = ENumber(3276700);
        outParam = ENumber(-3276700);
        intFace=intEdge=false;
        for (int i=0;i<3;i++){
            Segment2 edgeSegment(triangle[i], triangle[(i+1)%3]);
            std::vector<ENumber> result = line_segment_intersection(line, edgeSegment);
            for (int j=0;j<result.size();j++){
                inParam = (inParam < result[j] ? inParam : result[j]);
                outParam = (outParam > result[j] ? outParam : result[j]);
            }
            if (result.size()==2){
                std::cout<<"line triangle overlap!!!"<<std::endl;
                intEdge=true;
                intFace=false;
                return;
            }
            if (result.size()==1)
                intFace=true;
        }
        if (inParam==outParam)  //intersecting the triangle only by a vertex; ignored
            intFace=intEdge=false;
        /*if (intFace){
         std::cout<<"Intersecting within the face with parameters "<<inParam.get_d()<<"->"<<outParam.get_d()<<std::endl;
         }*/
        /*else if (intEdge){
         std::cout<<"Intersecting an edge with parameters "<<inParam.get_d()<<"->"<<outParam.get_d()<<std::endl;
         } else std::cout<<"Line doesn't intersect triangle"<<std::endl;*/
        
    }
    
    void linepencil_triangle_intersection(const LinePencil& lp,
                                          const std::vector<EVector2> triangle,
                                          std::vector<bool>& intEdges,
                                          std::vector<bool>& intFaces,
                                          std::vector<ENumber>& inParams,
                                          std::vector<ENumber>& outParams,
                                          std::vector<std::vector<ENumber>>& triParams){
        
        using namespace std;
        inParams.resize(lp.numLines);
        outParams.resize(lp.numLines);
        intEdges.resize(lp.numLines);
        intFaces.resize(lp.numLines);
        for (int i=0;i<inParams.size();i++){
            inParams[i] = ENumber(3276700);
            outParams[i] = ENumber(-3276700);
            intEdges[i]=intFaces[i]=false;
        }
        
        directional::EVector2 p00, pVec1, pVec2;
        ENumber t1p00, t2p00, dt1, dt2;
        triParams.resize(3);
        for (int i=0;i<3;i++){
            LinePencil triEdgePencil;
            triEdgePencil.numLines = 1;
            triEdgePencil.p0 = triangle[i];
            triEdgePencil.direction = triangle[(i+1)%3] - triangle[i];
            triEdgePencil.pVec[0] = -triEdgePencil.direction[1];
            triEdgePencil.pVec[1] = triEdgePencil.direction[0];
            Eigen::Matrix<ENumber, 2, 1> t00;
            Eigen::Matrix<ENumber, 2, 2> I2dt;
            EInt iso1Overlap;
            int lpResult = linepencil_intersection(lp, triEdgePencil, t00, I2dt, iso1Overlap);
            if (lpResult==2){  //overlap of iso1Overlap isovalue with the segment
                std::cout<<"line triangle overlap!!!"<<std::endl;
                intEdges[iso1Overlap.convert()]=true;
                //updating only the params of this line, not those of the others who stay open
                ENumber tSource = lp.line(iso1Overlap.convert()).point_param(triangle[i]);
                ENumber tTarget = lp.line(iso1Overlap.convert()).point_param(triangle[(i+1)%3]);
                if (tSource < tTarget){
                    inParams[iso1Overlap.convert()]=tSource; outParams[iso1Overlap.convert()]=tTarget;
                } else {
                    inParams[iso1Overlap.convert()]=tTarget; outParams[iso1Overlap.convert()]=tSource;
                }
                triParams[i].push_back(ENumber(0));
                triParams[i].push_back(ENumber(1));
                //if (tsource==tTarget)  //intersecting the triangle edge only by a vertex; ignored
                //    intFace[iso1Overlap]=intEdge[iso1Overlap]=false;
            }
            if (lpResult==1){  //intersecting the triangle edge non-trivially
                //TODO: here
                //generating the intersection t pairs for all lines
                Eigen::Matrix<ENumber, 2, Eigen::Dynamic> IPairs(2, lp.numLines);
                Eigen::Matrix<ENumber, 2, Eigen::Dynamic> tPairs(2, lp.numLines);
                //cout<<"tPairs: "<<endl;
                for (int j=0;j<lp.numLines;j++){
                    //Eigen::Matrix<ENumber, 2, 1> currI; currI<<ENumber(j,0.0), ENumber(0);
                    tPairs.col(j) = t00+I2dt.col(0)*ENumber(j);//*currI;  //can be replaced by a simple extraction from the matrix
                    //cout<<tPairs(0,j).to_double()<<", "<<tPairs(1,j).to_double()<<endl;
                    //testing:
                    //Line2 l1 = lp.line(j);
                    //Line2 l2 = triEdgePencil.line(0);
                    //ENumber t1, t2;
                    //line_line_intersection(l1, l2, t1, t2);
                    //cout<<"t1: "<<t1.to_double()<<" t2: "<<t2.to_double()<<endl;
                }
                ENumber tTriSource = triEdgePencil.line(0).point_param(triangle[i]);
                ENumber tTriTarget = triEdgePencil.line(0).point_param(triangle[(i+1)%3]);
                //cout<<"Intersecting edge " <<i<<endl;
                //cout<<"tTriSource: "<<tTriSource.to_double()<<endl;
                //cout<<"tTriTarget: "<<tTriTarget.to_double()<<endl;
                for (int j=0;j<lp.numLines;j++){
                    if ((tPairs(1,j)>=tTriSource)&&(tPairs(1,j)<=tTriTarget)){
                        triParams[i].push_back(tPairs(1,j));
                        if (tPairs(0,j)>outParams[j]) outParams[j]=tPairs(0,j);
                        if (tPairs(0,j)<inParams[j]) inParams[j]=tPairs(0,j);
                        
                    }
                }
            }
        }
        
        //for (int i=0;i<lp.numLines;i++)
        //    assert("didn't find an entry and exit points" && inParams[i]!=ENumber(3276700) && outParams[i]!=ENumber(-3276700));
        //Setting everything that has nontrivial inParams->outParams and not intEdge to intFaces = true
        for (int i=0;i<lp.numLines;i++)
            if ((!intEdges[i])&&(inParams[i]!=outParams[i])&&(inParams[i]!=ENumber(3276700))&&(outParams[i]!=ENumber(-3276700)))
                intFaces[i]=true;
    }
    
    //according to this: https://math.stackexchange.com/questions/1450498/rational-ordering-of-vectors
    ENumber slope_function(const EVector2& vec){
        //predicates might be expensive, so precomputing
        bool x0 = vec[0]>ENumber(0);
        bool y0 = vec[1]>ENumber(0);
        //bool xy = (y0 && vec[1]>vec[0])||(!y0 && vec[1]<=vec[0]);
        bool xy = vec[1].abs() > vec[0].abs();
        
        //std::cout<<"vec: "<<vec<<std::endl;
        
        if (xy){
            if (y0) return ENumber(1) - vec[0]/vec[1]; // case 1
            else return ENumber(5) - vec[0]/vec[1];  //case 3
        } else {
            if (x0){
                if (y0) return vec[1]/vec[0] - ENumber(1);  //case 0
                else return vec[1]/vec[0] + ENumber(7); //case 4
            }else{
                return vec[1]/vec[0] + ENumber(3); //case 2
            }
        }
    }
    
    
    double slope_function_double(const Eigen::RowVector2d& vec){
        //predicates might be expensive, so precomputing
        bool x0 = vec[0]>0.0;
        bool y0 = vec[1]>0.0;
        //bool xy = (y0 && vec[1]>vec[0])||(!y0 && vec[1]<=vec[0]);
        bool xy = std::abs(vec[1]) > std::abs(vec[0]);
        
        //std::cout<<"vec: "<<vec<<std::endl;
        
        if (xy){
            if (y0) return 1.0 - vec[0]/vec[1]; // case 1
            else return 5.0 - vec[0]/vec[1];  //case 3
        } else {
            if (x0){
                if (y0) return vec[1]/vec[0] - 1.0;  //case 0
                else return vec[1]/vec[0] + 7.0; //case 4
            }else{
                return vec[1]/vec[0] + 3.0; //case 2
            }
        }
    }
    
    
    
    
    double signed_face_area(const std::vector<EVector2>& faceVectors)
    {
        Eigen::RowVector2d currVertex=Eigen::RowVector2d::Zero(); //currVertex[0]=currVertex[1]=ENumber(0);
        double sfa=0.0;
        for (int i=0;i<faceVectors.size();i++){
            Eigen::RowVector2d nextVector = faceVectors[i].to_double();
            Eigen::RowVector2d nextVertex = currVertex+nextVector;
            sfa = sfa + currVertex[0]*nextVertex[1]-currVertex[1]*nextVertex[0];
            //std::cout<<"curr sfa: "<<sfa.get_d()<<std::endl;
            currVertex=nextVertex;
        }
        return sfa;
    }
    
    ENumber triangle_area(const std::vector<EVector2>& tri){
        EVector2 e12 = tri[1]-tri[0];
        EVector2 e13 = tri[2]-tri[0];
        return (e12[0]*e13[1]-e13[0]*e12[1])/ENumber(2);
    }
    
    void div_mod(const EInt a, const EInt b, EInt& q, EInt& r){
        //mpz_tdiv_qr(q.get_mpz_t(),r.get_mpz_t(),a,b);
        q = a / b;
        r = a - b * q;
    }
}







#endif //DIRECTIONAL_GMP_DEFINITIONS_H
