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
        mpq_class& operator[](size_t index) {
            return data[index];
        }

    private:
        std::vector<mpq_class> data;
    };

    typedef EVector<2> EVector2;
    typedef EVector<3> EVector3;

    ENumber squaredDistance(const EVector3& v1, const EVector3& v2){
        ENumber sd(0);
        //assert("vectors are of different dimensions!" && v1.size()==v2.size());
        for (int i=0;i<3;i++)
            ENumber+=(v1[i]-v2[i])*(v1[i]-v2[i]);  //maybe it's not efficient

        return ENumber;
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
}

#endif //DIRECTIONAL_GMP_DEFINITIONS_H
