//
// Created by Amir Vaxman on 20.04.24.
//

#ifndef DIRECTIONAL_GMP_DEFINITIONS_H
#define DIRECTIONAL_GMP_DEFINITIONS_H

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

namespace directional{
    struct MpqAllocator {
        typedef mpq_t value_type;
        mpq_t* allocate(size_t n) const {
            return static_cast<mpq_t*>(malloc(n * sizeof(mpq_t)));
        }
        void deallocate(mpq_t* p, size_t) const {
            free(p);
        }
    };

    typedef mpz_class EInt;
    typedef mpq_class ENumber;
    //typedef Eigen::Matrix<mpq_t, Eigen::Dynamic, Eigen::Dynamic> MatrixXed;

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
                newVec.data[i]=data[i]+ev.data[i];

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
                if (data[i]>=ev.data[i])
                    return false;
            return true;
        }

        EVector(const std::initializer_list<ENumber>& args){
            data.insert(data.end(), args.begin(), args.end());
        }

        EVector operator=(const EVector<Size>& evec){
            data=evec.data;
            return *this;
        }

        template<size_t _Size>
        friend std::ostream& operator<<(std::ostream& os, const EVector<_Size>& evec);

    protected:
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
            os<< evec[i].get_d()<<",";
        os<<evec[Size-1].get_d()<<")";
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
    };

    std::ostream& operator<<(std::ostream& os, const Line2& line) {
        os << "Line2(" << line.point << " + " << line.direction << ")";
        return os;
    }

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
            for (Eigen::SparseMatrix<int>::InnerIterator it(M,k); it; ++it)
                y[it.row()]+=ENumber((long)it.value())*x[it.col()];
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

        std::queue<int> nextVertexQueue;
        for (int i=0;i<components.size();i++)
            nextVertexQueue.push(i);

        int numComponents=0;
        while (!nextVertexQueue.empty()) {
            int nextVertex=nextVertexQueue.front();
            nextVertexQueue.pop();
            if (components[nextVertex]==-1)  //first components
                components[nextVertex]=numComponents++;

            //Otherwise, doing DFS on edges
            for (int i=0;i<VV[nextVertex].size();i++){
                if (components[VV[nextVertex][i]]==-1){
                    components[VV[nextVertex][i]]=components[nextVertex];
                    nextVertexQueue.push(VV[nextVertex][i]);
                } else assert(components[VV[nextVertex][i]]==components[nextVertex]);
            }
        }
        return numComponents-1;
    }

    int line_line_intersection(const Line2& line1,
                               const Line2& line2,
                               ENumber& t1,
                               ENumber& t2){
        std::cout<<"Computing intersection between line: "<<line1<<" and line "<<line2<<std::endl;
        ENumber denom = line1.direction[0]*line2.direction[1]-line1.direction[1]*line2.direction[0];
        if (denom==0){
            EVector2 pointVec = line1.point-line2.point;
            return  (pointVec[0]*line1.direction[1]-pointVec[1]*line1.direction[0]==0 ? 2 : 0);
        }

        t1 = ((line2.point[0]-line1.point[0])*(line2.direction[1])-(line2.point[1]-line1.point[1])*(line2.direction[0]))/denom;
        t2 = ((line2.point[0]-line1.point[0])*(line1.direction[1])-(line2.point[1]-line1.point[1])*(line1.direction[0]))/denom;
        assert("line_line_intersection is wrong!" && line1.point+t1*line1.direction == line2.point+t2*line2.direction);
        return 1;

    }

    std::vector<std::pair<ENumber, ENumber>> segment_segment_intersection(const Segment2& seg1,
                                                             const Segment2& seg2){

        ENumber t1, t2;
        int result = line_line_intersection(Line2(seg1.source, seg1.target-seg1.source),
                                            Line2(seg2.source, seg2.target-seg2.source),t1, t2);

        if (result==0)
            return std::vector<std::pair<ENumber, ENumber>>(); //no intersection

        if (result==1) {  //a single intersection at most; should check t1 and t2
            if ((t1>=ENumber(0))&&(t1<=ENumber(1))&&(t2>=ENumber(0))&&(t2<=ENumber(1))){
                std::vector<std::pair<ENumber, ENumber>> point(1);
                point[0]=std::pair<ENumber, ENumber>(t1,t2);
                return point;
            }
        }

        if (result==2){  //lines overlap; should check the segments overlap and then return both overlap points (order not important)
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
                return points;
            }

        }


    }



    std::vector<ENumber> line_segment_intersection(const Line2& line,
                                                   const Segment2& segment){

        std::cout<<"Computing intersection between line :"<<line<<" and segment "<<segment<<std::endl;
        Line2 segLine(segment.source, segment.target-segment.source);
        ENumber t1, t2;
        int intersectType=line_line_intersection(line, segLine, t1, t2);
        if (intersectType==0)   //no intersection
            return std::vector<ENumber>();
        if (intersectType==2) { //the entire segment is contained in the line
            std::vector<ENumber> result(2); result[0]=ENumber(0); result[1]=ENumber(1);
            return result;
        }
        if (intersectType==1){
            std::vector<ENumber> result(1); result[0]=t2;
            std::cout<<"t2: "<<t2.get_d()<<std::endl;
            return result;
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
                intEdge=true;
                intFace=false;
                return;
            }
            if (result.size()==1)
                intFace=true;
        }

    }

    //according to this: https://math.stackexchange.com/questions/1450498/rational-ordering-of-vectors
    ENumber slope_function(const EVector2& vec){
        //predicates might be expensive, so precomputing
        bool x0 = vec[0]>0;
        bool y0 = vec[1]>0;
        //bool xy = (y0 && vec[1]>vec[0])||(!y0 && vec[1]<=vec[0]);
        bool xy = abs(vec[1])>abs(vec[0]);

        std::cout<<"vec: "<<vec<<std::endl;

        if (xy){
            if (y0) return (vec[1]-vec[0])/(vec[1]); // case 1
            else return (vec[1]-vec[0])/(vec[1])+ENumber(4);  //case 3
        } else {
            if (x0){
                if (y0) return (vec[1]-vec[0])/(vec[0]);  //case 0
                else return (vec[1]-vec[0])/(vec[0])+ENumber(8); //case 4
            }else{
                return (vec[1]-vec[0])/(vec[0])+ENumber(4); //case 2
            }
        }
    }


    ENumber signed_face_area(const std::vector<EVector2>& faceVectors)
    {
        EVector2 currVertex; currVertex[0]=currVertex[1]=ENumber(0);
        ENumber sfa(0);
        for (int i=0;i<faceVectors.size();i++){
            EVector2 nextVector = faceVectors[(i!=faceVectors.size()-1) ? i : 0];
            EVector2 nextVertex = currVertex+nextVector;
            sfa = sfa + currVertex[0]*nextVertex[1]-currVertex[1]*nextVertex[0];
            currVertex=nextVertex;
        }
        return sfa;
    }

    ENumber triangle_area(const std::vector<EVector2>& tri){
       EVector e12 = tri[1]-tri[0];
       EVector e13 = tri[2]-tri[0];
       return (e12[0]*e13[1]-e13[0]*e12[1])/ENumber(2);
    }

    void div_mod(const mpz_ptr a, mpz_ptr b, EInt& q, EInt& r){
        mpz_tdiv_qr(q.get_mpz_t(),r.get_mpz_t(),a,b);
    }
}

#endif //DIRECTIONAL_GMP_DEFINITIONS_H
