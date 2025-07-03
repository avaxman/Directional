#include <iostream>
#include <Eigen/Core>
#include <directional/readOBJ.h>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/gradient_matrices.h>
#include <directional/mass_matrices.h>
#include <directional/curl_matrices.h>
#include <directional/div_matrices.h>

directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField exactField, coexactField;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOBJ(TUTORIAL_DATA_PATH "/60745__sf.obj",mesh);
    ftb.init(mesh);
    
    Eigen::VectorXd confVec(mesh.dcel.vertices.size()), nonConfVec(mesh.dcel.edges.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        confVec[i] = 50.0*sin(mesh.V(i,2)/20.0);
    for (int i=0;i<mesh.dcel.edges.size();i++)
        nonConfVec[i] = (mesh.midEdges.row(i)-mesh.midEdges.row(17160)).norm();
    
    Eigen::SparseMatrix<double> Gv = directional::conf_gradient_matrix_2D<double>(mesh);
    Eigen::SparseMatrix<double> Ge = directional::non_conf_gradient_matrix_2D<double>(mesh);
    Eigen::SparseMatrix<double> J =  directional::face_vector_rotation_matrix_2D<double>(mesh);
    
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(mesh);
    Eigen::SparseMatrix<double> D = directional::div_matrix_2D<double>(mesh);
    
    //demonstrating the exact sequences
    std::cout<<"max abs curl of exact field (should be numerically zero): "<<(C*Gv*confVec).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"max abs divergence of coexact field (should be numerically zero): "<<(D*J*Ge*nonConfVec).cwiseAbs().maxCoeff()<<std::endl;
    
    Eigen::MatrixXd gradRawField(mesh.F.rows(), 3), rotCogradRawField(mesh.F.rows(), 3);
    Eigen::VectorXd gradFieldVec = Gv*confVec;
    Eigen::VectorXd rotCogradFieldVec = J*Ge*nonConfVec;
    for (int i=0;i<mesh.F.rows();i++){
        gradRawField.row(i)<<gradFieldVec.segment(3*i,3).transpose();
        rotCogradRawField.row(i)<<rotCogradFieldVec.segment(3*i,3).transpose();
    }
    
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_surface_edge_data(nonConfVec, "Non-Conforming Function");
    viewer.set_surface_vertex_data(confVec ,"Conforming Function");
    viewer.set_raw_field(mesh.barycenters, gradRawField, "Gradient field", 0, 1.0*mesh.avgEdgeLength);
    viewer.set_raw_field(mesh.barycenters, rotCogradRawField,  "Rot. Cogradient field", 1, 1.0*mesh.avgEdgeLength);
    viewer.toggle_raw_field(false, 1);
    viewer.launch();
}
