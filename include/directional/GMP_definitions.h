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
        }
        EVector<Size> operator-() const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=-data[i];
        }
        EVector<Size> operator-(const EVector<Size>& ev) const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=data[i]+ev.data[i];
        }

        EVector<Size> operator*(const ENumber s) const{
            EVector<Size> newVec;
            for (int i=0;i<Size;i++)
                newVec.data[i]=data[i]*s;
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

    protected:
        std::vector<ENumber> data;
    };

    template<size_t Size>
    EVector<Size> operator*(ENumber scalar, const EVector<Size>& vec) {
        return vec * scalar; // Leverage the previous operator*
    }

    typedef EVector<2> EVector2;
    typedef EVector<3> EVector3;

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

    int line_line_intersection(const std::pair<EVector2, EVector2>& line1,
                               const std::pair<EVector2, EVector2>& line2,
                               ENumber t1,
                               ENumber t2){
        ENumber denom = line1.second[0]*line2.second[1]-line1.second[0]*line2.second[0];
        if (denom==0)
            return (line1.first==line2.first ? 2 : 0);

        t1 = ((line2.first[0]-line1.first[0])*(line2.second[1])-(line2.first[1]-line1.first[1])*(line2.second[0]))/denom;
        t2 = ((line2.first[0]-line1.first[0])*(line1.second[1])-(line2.first[1]-line1.first[1])*(line1.second[0]))/denom;
        assert("line_line_intersection is wrong!" && line1.first+t1*line1.second == line2.first+t2*line2.second);

    }

    std::vector<std::pair<ENumber, ENumber>> segment_segment_intersection(const std::pair<EVector2, EVector2>& seg1,
                               const std::pair<EVector2, EVector2>& seg2){

        ENumber t1, t2;
        int result = line_line_intersection(std::pair<EVector2, EVector2>(seg1.first, seg1.second-seg1.first),
                                            std::pair<EVector2, EVector2>(seg2.first, seg2.second-seg2.first),t1, t2);

        if (result==0)
            return std::vector<std::pair<ENumber, ENumber>>(); //no intersection

        if (result==1) {  //a single intersection at most; should check t1 and t2
            if ((t1>=ENumber(0))&&(t1<=ENumber(1))&&(t2>=ENumber(0))&&(t2<=ENumber(1))){
                std::vector<std::pair<ENumber, ENumber>> result(1);
                result[0]=std::pair<ENumber, ENumber>(t1,t2);
                return result;
            }
        }

        if (result==2){  //lines overlap; should check the segments overlap and then return both overlap points (order not important)
            EVector2 vec = seg1.second-seg1.first;
            int axis = (vec[0]!=ENumber(0) ? 0 : 1);
            std::pair<EVector2, EVector2> sortSeg1, sortSeg2;
            if (seg1.first[axis]<seg1.second[axis]) sortSeg1=seg1; else sortSeg1=std::pair<EVector2, EVector2>(seg1.second,seg1.first);
            if (seg2.first[axis]<seg2.second[axis]) sortSeg2=seg2; else sortSeg2=std::pair<EVector2, EVector2>(seg2.second,seg2.first);
            EVector2 startPoint = (sortSeg1.first[axis]>sortSeg2.first[axis] ? sortSeg1.first : sortSeg2.first);
            EVector2 endPoint = (sortSeg1.second[axis]>sortSeg2.second[axis] ? sortSeg2.second : sortSeg1.second);

            //TODO: what is the policy with segment-vertex tangency?
            if (startPoint[axis]<endPoint[axis]) { //there is a (non-zero) intersection
                ENumber startAtSeg1 = (startPoint[axis]-seg1.first[axis])/(seg1.second[axis]-seg1.first[axis]);
                ENumber startAtSeg2 = (startPoint[axis]-seg2.first[axis])/(seg2.second[axis]-seg2.first[axis]);
                ENumber endAtSeg1 = (endPoint[axis]-seg1.first[axis])/(seg1.second[axis]-seg1.first[axis]);
                ENumber endAtSeg2 = (endPoint[axis]-seg2.first[axis])/(seg2.second[axis]-seg2.first[axis]);
                std::vector<std::pair<ENumber, ENumber>> result(2);
                result[0] = std::pair<ENumber, ENumber>(startAtSeg1, startAtSeg2);
                result[1] = std::pair<ENumber, ENumber>(endAtSeg1, endAtSeg2);
                return result;
            }

        }


    }



    std::vector<ENumber> line_segment_intersection(const std::pair<EVector2, EVector2>& line,
                                   const std::pair<EVector2, EVector2>& segment){

        std::pair<EVector2, EVector2> segLine(segment.first, segment.second-segment.first);
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
            return result;
        }
    }

    void line_triangle_intersection(const std::pair<EVector2, EVector2>& line,
                                    const std::vector<EVector2> triangle,
                                    bool intEdge,
                                    bool intFace,
                                    ENumber inParam,
                                    ENumber outParam){

        inParam = ENumber(3276700);
        outParam = ENumber(-3276700);
        intFace=intEdge=false;
        for (int i=0;i<3;i++){
            std::pair<EVector2, EVector2> edgeSegment(triangle[0], triangle[1]);
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

    ENumber slope_function(const EVector2& vec){

    }
}

#endif //DIRECTIONAL_GMP_DEFINITIONS_H
