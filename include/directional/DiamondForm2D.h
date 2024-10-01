// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_DIAMOND_FORM_2D_H
#define DIRECTIONAL_DIAMOND_FORM_2D_H

#include <eigen/sparse>
#include <directional/VolumeForm2D.h>
#include <directional/TriMesh.h>

namespace directional{

    //an interface to a single element in the class
    template<typename NumberType>
    class DiamondForm2D : public VolumeForm2D<NumberType>{
    public:

        int N;
        const TriMesh* mesh;
        Eigen::Vector<NumberType, Eigen::Dynamic> diamondValues;

        DiamondForm2D(){}
        ~DiamondForm2D(){}


        void init(const TriMesh* _mesh, Eigen::Vector<NumberType, Eigen::Dynamic>& _diamondValues, const int _N=1){
            N=_N;
            mesh=_mesh;
            diamondValues=_diamondValues;
        }

        Eigen::SparseMatrix<NumberType> mass_matrix(){
            assert("Currently only defined for N=1 " && N==1);
            Eigen::VectorXd diamondAreas(mesh->EF.rows());
            for (int i=0;i<mesh->EF.rows();i++){
                double diamondArea=0.0;
                if (mesh->EF(i,0)!=-1)
                    diamondArea+=mesh->faceAreas(mesh->EF(i,0));
                if (mesh->EF(i,1)!=-1)
                    diamondArea+=mesh->faceAreas(mesh->EF(i,1));

                diamondAreas(i)=diamondArea/3.0;
            }

        }
    };

}


#endif
