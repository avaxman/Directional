// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PIECEWISE_LINEAR_FUNCTION_H
#define DIRECTIONAL_PIECEWISE_LINEAR_FUNCTION_H

#include <eigen/sparse>
#include <directional/ScalarFunction2D.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/gradient_matrix.h>

namespace directional{

    template<typename NumberType>
    class PLFunction2D: public ScalarFunction2D<NumberType>{
    public:

        int N;  //the order of the field
        const TriMesh* mesh;
        Eigen::Vector<NumberType, Eigen::Dynamic> vertexValues;

        //nodeValues are values on vertices
        PLFunction2D(){}
        ~PLFunction2D(){}

        void init(const TriMesh& _mesh, Eigen::Vector<NumberType, Eigen::Dynamic>& _vertexValues, const int _N=1){
            N=_N;
            mesh=&_mesh;
            vertexValues=_vertexValues;
        }

        //Only good for triangle meshes right now
        NumberType  value(const int faceIndex,
                         const Eigen::VectorXd& baryCoords){

            NumberType currValue = 0;
            for (int i=0;i<baryCoords.size();i++)
                currValue += baryCoords[i]*vertexValues(mesh->F(faceIndex,i));

            return currValue;
        }

        //TODO: save the matrix for future use
        Eigen::SparseMatrix<NumberType> gradient_matrix(){
            return directional::gradient_matrix<NumberType>(mesh, 1, N);
        }

        void gradient(directional::CartesianField& gradField){
            assert("PLFunction::gradient(): gradField is no of the correct type! " && gradField.N==N && gradField.tb->discTangType()==discTangTypeEnum::FACE_SPACES);
            gradField.fieldType=fieldTypeEnum::RAW_FIELD;
            Eigen::SparseMatrix<NumberType> G = gradient_matrix();
            Eigen::Vector<NumberType, Eigen::Dynamic> fieldVector = G*vertexValues;
            Eigen::Matrix<NumberType, Eigen::Dynamic, 3> extField(fieldVector.size()/(3*N), 3*N);
            for (int i=0;i<extField.rows();i++)
                extField.row(i)<<fieldVector.segment(3*N*i,3*N).transpose();
            gradField.set_extrinsic_field(extField);
        }

        //the original unlumped mass matrix of inner product of hat functions
        Eigen::SparseMatrix<NumberType> mass_matrix(){
            assert("Currently only implemented for N=1" && N==1);
            Eigen::SparseMatrix<NumberType> M(vertexValues.size(), vertexValues.size());
            std::vector<Eigen::Triplet<NumberType>> MTris;
            for (int i=0;i<mesh->F.rows();i++){
                for (int j=0;j<3;j++)
                    for (int k=0;k<3;k++)
                        MTris.push_back(Eigen::Triplet<NumberType>(mesh->F(i,j),mesh->F(i,k), mesh->faceAreas(i)*(j==k ? 1.0/6.0 : 1.0/12.0)));
            }
            M.setFromTriplets(MTris.begin(), MTris.end());
            return M;
        }

        //uses the lumped matrix which is not a true inverse
        Eigen::SparseMatrix<NumberType> inv_mass_matrix(){
            assert("Currently only implemented for N=1" && N==1);
            Eigen::SparseMatrix<NumberType> M(vertexValues.size(), vertexValues.size());
            std::vector<Eigen::Triplet<NumberType>> MTris;
            for (int i=0;i<mesh->F.rows();i++){
                //adding the 1/3 of each face's area to the vertex
                for (int j=0;j<3;j++)
                    MTris.push_back(Eigen::Triplet<NumberType>(mesh->F(i,j),mesh->F(i,j), 3.0/mesh->faceAreas(i)));
            }
            M.setFromTriplets(MTris.begin(), MTris.end());
            return M;
        }
    };

}


#endif
