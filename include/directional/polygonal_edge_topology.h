// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_POLYGONAL_EDGE_TOPOLOGY_H
#define HEDRA_POLYGONAL_EDGE_TOPOLOGY_H

#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <vector>


namespace hedra
{
    // Initialize Edges and their topological relations
    
    //input:
    //  D  eigen int vector     #F by 1 - face degrees
    //  F  eigen int matrix     #F by max(D) - vertex indices in face
  
    // Output:
    // EV   #E by 2, Stores the edge description as pair of indices to vertices
    // FE : #F by max(D), Stores the Face-Edge relation
    // EF : #E by 2: Stores the Edge-Face relation
    // EFi: #E by 2: corresponding to EF and stores the relative position of the edge in the face (e.g., if the edge is (v1,v2) and the face has (vx,vy,v2,v1,vz,va), then the value is 3)
    // EFs: #E by 2: if the edge is oriented positively or negatively in the face (e.g. in the example above we get -1)
    // InnerEdges: indices into EV of which edges are internal (not boundary)

    IGL_INLINE void polygonal_edge_topology(const Eigen::VectorXi& D,
                                            const Eigen::MatrixXi& F,
                                            Eigen::MatrixXi& EV,
                                            Eigen::MatrixXi& FE,
                                            Eigen::MatrixXi& EF,
                                            Eigen::MatrixXi& EFi,
                                            Eigen::MatrixXd& FEs,
                                            Eigen::VectorXi& InnerEdges)
    {
        // Only needs to be edge-manifold
        std::vector<std::vector<int> > ETT;
        for(int f=0;f<D.rows();++f)
            for (int i=0;i<D(f);++i)
            {
                // v1 v2 f vi
                int v1 = F(f,i);
                int v2 = F(f,(i+1)%D(f));
                if (v1 > v2) std::swap(v1,v2);
                std::vector<int> r(4);
                r[0] = v1; r[1] = v2;
                r[2] = f;  r[3] = i;
                ETT.push_back(r);
            }
        std::sort(ETT.begin(),ETT.end());
        
        // count the number of edges (assume manifoldness)
        int En = 1; // the last is always counted
        for(unsigned i=0;i<ETT.size()-1;++i)
            if (!((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1])))
                ++En;
        
        EV = Eigen::MatrixXi::Constant((int)(En),2,-1);
        FE = Eigen::MatrixXi::Constant((int)(F.rows()),(int)(F.cols()),-1);
        EF = Eigen::MatrixXi::Constant((int)(En),2,-1);
        En = 0;
        
        for(unsigned i=0;i<ETT.size();++i)
        {
            if (i == ETT.size()-1 ||
                !((ETT[i][0] == ETT[i+1][0]) && (ETT[i][1] == ETT[i+1][1]))
                )
            {
                // Border edge
                std::vector<int>& r1 = ETT[i];
                EV(En,0)     = r1[0];
                EV(En,1)     = r1[1];
                EF(En,0)    = r1[2];
                FE(r1[2],r1[3]) = En;
            }
            else
            {
                std::vector<int>& r1 = ETT[i];
                std::vector<int>& r2 = ETT[i+1];
                EV(En,0)     = r1[0];
                EV(En,1)     = r1[1];
                EF(En,0)    = r1[2];
                EF(En,1)    = r2[2];
                FE(r1[2],r1[3]) = En;
                FE(r2[2],r2[3]) = En;
                ++i; // skip the next one
            }
            ++En;
        }
        
        // Sort the relation EF, accordingly to EV
        // the first one is the face on the left of the edge
        
        for(unsigned i=0; i<EF.rows(); ++i)
        {
            int fid = EF(i,0);
            bool flip = true;
            // search for edge EV.row(i)
            for (unsigned j=0; j<D(fid); ++j)
            {
                if ((F(fid,j) == EV(i,0)) && (F(fid,(j+1)%D(fid)) == EV(i,1)))
                    flip = false;
            }
            
            if (flip)
            {
                int tmp = EF(i,0);
                EF(i,0) = EF(i,1);
                EF(i,1) = tmp;
            }
        }
        
        
        std::vector<int> InnerEdgesVec;
        EFi=Eigen::MatrixXi::Constant(EF.rows(), 2,-1);
        FEs=Eigen::MatrixXd::Zero(FE.rows(),FE.cols());
        for (int i=0;i<EF.rows();i++)
            for (int k=0;k<2;k++){
                if (EF(i,k)==-1)
                    continue;
                    
                for (int j=0;j<D(EF(i,k));j++)
                    if (FE(EF(i,k),j)==i)
                        EFi(i,k)=j;
            }
        
        for (int i=0;i<EF.rows();i++){
            if (EFi(i,0)!=-1) FEs(EF(i,0),EFi(i,0))=1.0;
            if (EFi(i,1)!=-1) FEs(EF(i,1),EFi(i,1))=-1.0;
            if ((EF(i,0)!=-1)&&(EF(i,1)!=-1))
                InnerEdgesVec.push_back(i);
        }
        
        InnerEdges.resize(InnerEdgesVec.size());
        for (int i=0;i<InnerEdgesVec.size();i++)
            InnerEdges(i)=InnerEdgesVec[i];
        

    }
}


#endif
