// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_BARYCENTRIC_SUBDIVISION_H
#define DIRECTIONAL_BARYCENTRIC_SUBDIVISION_H


#include <Eigen/Core>
#include <vector>

namespace directional
{

//Subdivides a mesh trivially in the middle
inline void barycentric_subdivision(const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    Eigen::MatrixXd& VFine,
                                    Eigen::MatrixXi& FFine)
{
    VFine.resize(V.rows()+F.rows(),3);
    VFine.block(0,0,V.rows(),3);
    for (int i=0;i<FFine.rows();i++)
        VFine.row(V.rows()+i) = (V.row(F(i,0))+V.row(F(i,1))+V.row(F(i,2)))/3.0;
    
    FFine.resize(3*F.rows(),3);
    for (int i=0;i<F.rows();i++){
        FFine.row(3*i)<<F(i,0), F(i,1), V.rows()+i;
        FFine.row(3*i+1)<<F(i,1), F(i,2), V.rows()+i;
        FFine.row(3*i+2)<<F(i,2), F(i,0), V.rows()+i;
    }
    
}
}




#endif


