#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/CochainComplex.h>
#include <directional/gradient_matrices.h>
#include <directional/curl_matrices.h>
#include <directional/mass_matrices.h>

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
Eigen::VectorXd origField, gradFieldVec, rotCogradFieldVec, harmFieldVec, vertexFunction, diamondForm;
directional::CartesianField gradField, rotCogradField, harmField;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/1146164.off",mesh);
    ftb.init(mesh);

    //TODO: create the actual field

    Eigen::SparseMatrix<double> G = directional::conf_gradient_matrix_2D<double>(&mesh);
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(&mesh);
    Eigen::SparseMatrix<double> Mx = directional::face_vectors_mass_matrix_2D<double>(&mesh);
    Eigen::SparseMatrix<double> Mc = directional::edge_diamond_mass_matrix_2D<double>(&mesh);

    directional::hodge_decomposition(G, C, Mx, Mc, origField, vertexFunction, gradFieldVec, rotCogradFieldVec, diamondForm, harmFieldVec);

    gradField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    gradField.set_extrinsic_field(gradFieldVec);
    //std::cout<<"exactField.extField: "<<exactField.extField<<std::endl;
    rotCogradField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    rotCogradField.set_extrinsic_field(rotCogradFieldVec);
    harmField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    harmField.set_extrinsic_field(harmFieldVec);

    //triangle mesh setup
    viewer.init();
    viewer.set_mesh(mesh);
    viewer.set_vertex_data(vertexFunction, vertexFunction.minCoeff(), vertexFunction.maxCoeff(),"Vertex Function", 0);
    viewer.set_field(gradField,"Gradient Field", 0, 1, 2.0);
    viewer.set_edge_data(diamondForm, diamondForm.minCoeff(), diamondForm.maxCoeff(),"Diamond Function", 0);
    viewer.set_field(rotCogradField,"Rotated Cogradient Field", 0, 1, 20.0);
    viewer.set_field(harmField,"Harmonic Field", 0, 2, 20.0);
    viewer.launch();
}
