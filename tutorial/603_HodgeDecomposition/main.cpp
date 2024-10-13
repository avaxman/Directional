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
#include <directional/extrinsic_intrinsic_matrices.h>

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
Eigen::VectorXd gradFieldVec, rotCogradFieldVec, harmFieldVec, vertexFunction, diamondForm;
directional::CartesianField origField, gradField, rotCogradField, harmField;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/1146164.off",mesh);
    ftb.init(mesh);

    //TODO: create the actual field

    //Must use intrinsic since otherwise the harmonic field will have spurious normal components
    Eigen::SparseMatrix<double> G = directional::conf_gradient_matrix_2D<double>(&mesh, true);
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(&mesh, true);
    Eigen::SparseMatrix<double> Mx = directional::face_vectors_mass_matrix_2D<double>(&mesh, true);
    Eigen::SparseMatrix<double> Mc = directional::edge_diamond_mass_matrix_2D<double>(&mesh, true);
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(&mesh);

    Eigen::MatrixXd harmBasis;
    int bettiNumber;
    directional::cohomology_basis(G, C, Mx,  harmBasis, bettiNumber);
    std::cout<<"Euler characteristic: "<<mesh.V.rows()-mesh.EV.rows()+mesh.F.rows()<<std::endl;
    std::cout<<"Betti number: "<<bettiNumber<<std::endl;
    std::cout<<"harmBasis: "<<harmBasis.block(0,0,20,harmBasis.cols())<<std::endl;

    origField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    origField.set_extrinsic_field(IE*harmBasis.col(1));

    //directional::hodge_decomposition(G, C, Mx, Mc, origField, vertexFunction, gradFieldVec, rotCogradFieldVec, diamondForm, harmFieldVec);

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
    viewer.set_field(origField,"Original Field", 0, 0, 2.0);
    /*viewer.set_vertex_data(vertexFunction, vertexFunction.minCoeff(), vertexFunction.maxCoeff(),"Vertex Function", 0);
    viewer.set_field(gradField,"Gradient Field", 0, 1, 2.0);
    viewer.set_edge_data(diamondForm, diamondForm.minCoeff(), diamondForm.maxCoeff(),"Diamond Function", 0);
    viewer.set_field(rotCogradField,"Rotated Cogradient Field", 0, 2, 20.0);
    viewer.set_field(harmField,"Harmonic Field", 0, 3, 20.0);*/
    viewer.launch();
}
