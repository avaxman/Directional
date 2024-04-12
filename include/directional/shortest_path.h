//
// Created by Amir Vaxman on 10.04.24.
//

#ifndef DIRECTIONAL_SHORTEST_PATH_H
#define DIRECTIONAL_SHORTEST_PATH_H

#include <vector>
#include <set>
#include <queue>
#include <Eigen/Core>
#include <directional/TriMesh.h>

namespace directional{
    void inline shortest_path(const TriMesh& mesh,
                              const int initialVertex,
                              std::vector<int> targetNodes,
                              std::vector<int> shortestPath){

        using namespace std;
        using namespace Eigen;
        double infWeight = mesh.avgEdgeLength*(double)mesh.HV.size()*10000.0;
        Eigen::VectorXd distances = Eigen::VectorXd::Constant(mesh.V.rows(), infWeight);
        int currVertex = 0;
        distances(currVertex)=0.0;
        Eigen::VectorXi isVisited = Eigen::VectorXi::Zero(mesh.V.rows());
        Eigen::VectorXi targetNodeMask = Eigen::VectorXi::Zero(mesh.V.rows());;
        for (int i=0;i<targetNodes.size();i++)  //for easy access of having been visiting
            targetNodeMask(targetNodes[i]) = 1;

        if (targetNodeMask[initialVertex]==1)
            return; //initialvertex is already a target vertex


        int sumVisited = 0;
        set<pair<double, pair<int,int>> > distVertices;
        distVertices.insert(pair<double, pair<int, int>>(0.0,pair<int, int>(currVertex,-1)));
        Eigen::VectorXi predHEList = Eigen::VectorXi::Constant(mesh.V.rows(),-1);

        while ((sumVisited<=mesh.V.rows())&&(targetNodeMask[currVertex]==0)){
            //calculating distance to all neighbors
            assert("Trying an already visited vertex!" && (!isVisited[currVertex]));
            int hebegin = mesh.VH(currVertex);
            int hecurr = hebegin;
            do{
                int nextVertex = mesh.HV(mesh.nextH(hecurr));
                if (!isVisited(nextVertex)) {
                    double length = (mesh.V.row(nextVertex) - mesh.V.row(currVertex)).norm();
                    if (distances(nextVertex)>distances(currVertex)+length)
                        distances(nextVertex) = distances(currVertex)+length;

                    distVertices.insert(pair<double, pair<int, int>>(distances(nextVertex),pair<int,int>(nextVertex, hecurr)));
                }
                if (mesh.twinH(hecurr)==-1)
                    break;
                hecurr = mesh.twinH(mesh.prevH(hecurr));
            }while (hecurr!=hebegin);
            isVisited[currVertex]=1;
            sumVisited++;
            if (sumVisited>=mesh.V.rows())
                break;  //shouldn't happen! this means we haven't reached any target node!
            //going to next vertex with shortest distance
            int predHE = -1;
            do{
                currVertex = distVertices.begin()->second.first;
                predHE =distVertices.begin()->second.second;
                distVertices.erase(distVertices.begin());
            }while ((isVisited[currVertex])&&(!distVertices.empty()));
            predHEList[currVertex]=predHE;
        }

        //Backtracking the path
        assert("Haven't reached any target node!" && targetNodeMask[currVertex]==1);
        while (currVertex!=initialVertex){
            shortestPath.push_back(predHEList[currVertex]);
            currVertex = mesh.HV(predHEList[currVertex]);
        }

    }
}

#endif //DIRECTIONAL_SHORTEST_PATH_H
