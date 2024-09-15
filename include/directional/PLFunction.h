// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PIECEWISE_LINEAR_FUNCTION_H
#define DIRECTIONAL_PIECEWISE_LINEAR_FUNCTION_H

#include <eigen/sparse>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>

namespace directional{

    //an interface to a single element in the class
    template<typename NumberType>
    class PLFunction{
    public:

        const TriMesh* mesh;
        Eigen::Vector<NumberType, Eigen::Dynamic> vertexValues;

        PLFunction(){}
        ~PLFunction(){}

        void init(const TriMesh& _mesh, const Eigen::Vector<NumberType, Eigen::Dynamic> _vertexValues){
            mesh = &_mesh;
            vertexValues = _vertexValues;
        }

        //Only good for triangle meshes right now
        NumberType value(const int faceIndex,
                         const Eigen::VectorXd& baryCoords){

            NumberType currValue = 0;
            for (int i=0;i<baryCoords.size();i++)
                currValue += baryCoords[i]*vertexValues(mesh->F(faceIndex,i));

            return currValue;
        }

        //TODO: save the matrix for future use
        Eigen::SparseMatrix<NumberType> gradient_matrix(){
            Eigen::SparseMatrix<NumberType> G(3*mesh->F.rows(), vertexValues.size());
            std::vector<Eigen::Triplet<NumberType>> GTris;
            for (int i=0;i<mesh->F.rows();i++){
                Eigen::Matrix3d ep; ep<<mesh->faceNormals.row(i).cross(mesh->V.row(mesh->F(i, 2)) - mesh->V.row(mesh->F(i, 1))),
                mesh->faceNormals.row(i).cross(mesh->V.row(mesh->F(i, 0)) - mesh->V.row(mesh->F(i, 2))),
                mesh->faceNormals.row(i).cross(mesh->V.row(mesh->F(i, 1)) - mesh->V.row(mesh->F(i, 0)));
                double faceArea = mesh->faceAreas(i);
                for (int j=0;j<3;j++)
                    for (int k=0;k<3;k++)
                        GTris.push_back(Eigen::Triplet<NumberType>(3*i+k, mesh->F(i,j), ep(j,k)/(2.0*faceArea)));

            }
            G.set_from_triplets(GTris.begin(), GTris.end());
        }

        void gradient(directional::CartesianField& gradField){
            assert("PLFunction::gradient(): gradField is no of the correct type! " && gradField.N==1 && gradField.tv->discTangType()==discTangTypeEnum::FACE_SPACES);
            gradField.fieldType=fieldTypeEnum::RAW_FIELD;
            Eigen::SparseMatrix<NumberType> G = gradient_matrix();
            Eigen::Vector<NumberType, Eigen::Dynamic> fieldVector = G*vertexValues;
            Eigen::Matrix<NumberType, Eigen::Dynamic, 3> extField(fieldVector.size()/3, 3);
            for (int i=0;i<fieldVector.size()/3;i++)
                extField.row(i)<<fieldVector.segment(3*i,3).transpose();
            gradField.set_extrinsic_field(extField);
        }

        //the original unlumped mass matrix of inner product of hat functions
        Eigen::SparseMatrix<NumberType> mass_matrix(){
            Eigen::SparseMatrix<NumberType> M(vertexValues.size(), vertexValues.size());
            std::vector<Eigen::Triplet<NumberType>> MTris;
            for (int i=0;i<mesh->F.rows();i++){
                for (int j=0;j<3;j++)
                    for (int k=0;k<3;k++)
                        MTris.push_back(Eigen::Triplet<NumberType>(mesh->F(i,j),mesh->F(i,k), mesh->faceAreas(i)*(j==k ? 1.0/6.0 : 1.0/12.0)));
            }
            M.set_from_triplets(MTris.begin(), MTris.end());
            return M;
        }
        Eigen::SparseMatrix<NumberType> lumped_mass_matrix(){
            Eigen::SparseMatrix<NumberType> M(vertexValues.size(), vertexValues.size());
            std::vector<Eigen::Triplet<NumberType>> MTris;
            for (int i=0;i<mesh->F.rows();i++){
                //adding the 1/3 of each face's area to the vertex
                for (int j=0;j<3;j++)
                    MTris.push_back(Eigen::Triplet<NumberType>(mesh->F(i,j),mesh->F(i,j), mesh->faceAreas(i)/3.0));
            }
            M.set_from_triplets(MTris.begin(), MTris.end());
            return M;
        }
    };

}


#endif
