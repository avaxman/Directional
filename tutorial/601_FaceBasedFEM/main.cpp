#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/gradient_matrix.h>
#include <directional/mass_matrices.h>

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField exactField, coexactField;
directional::DirectionalViewer viewer;



int main()
{

  directional::readOFF(TUTORIAL_DATA_PATH "/Input_wheel.off",mesh);
  ftb.init(mesh);

  Eigen::VectorXd confVec(mesh.dcel.vertices.size()), nonConfVec(mesh.dcel.edges.size());
  for (int i=0;i<mesh.dcel.vertices.size();i++)
      confVec[i] = 10.0*sin(mesh.V(i,0))*cos(mesh.V(i,1));
  for (int i=0;i<mesh.dcel.edges.size();i++){
      Eigen::RowVector3d midEdgePoint = 0.5*(mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1)));
      nonConfVec[i] = 10.0*sin(midEdgePoint(1))*cos(midEdgePoint(2));
  }

  Eigen::SparseMatrix<double> Gv = directional::conf_gradient_matrix<double>(&mesh, 1,1);
  Eigen::SparseMatrix<double> Ge = directional::non_conf_gradient_matrix<double>(&mesh, 1,1);
  Eigen::SparseMatrix<double> J =  directional::face_vector_rotation_matrix<double>(&mesh, 3, 1);

  Eigen::SparseMatrix<double> C = directional::curl_matrix<double>(&mesh, 3, 1, 1);
  Eigen::SparseMatrix<double> D = directional::div_matrix<double>(&mesh, 3, 1, 1);

  //demonstrating the exact sequences
  std::cout<<"max abs curl of exact field (should be numerically zero): "<<(C*Gv*confVec).cwiseAbs().maxCoeff()<<std::endl;
  std::cout<<"max abs divergence of coexact field (should be numerically zero): "<<(D*J*Ge*nonConfVec).cwiseAbs().maxCoeff()<<std::endl;

    exactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    exactField.set_extrinsic_field(Gv*confVec);
    coexactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    coexactField.set_extrinsic_field(J*Ge*nonConfVec);

  //triangle mesh setup
  viewer.init();
  viewer.set_mesh(mesh);
  viewer.set_vertex_data(confVec, confVec.minCoeff(), confVec.maxCoeff(),"Conforming Function", 0);
  viewer.set_field(exactField,"Gradient field", 0, 0);
  viewer.set_edge_data(nonConfVec, nonConfVec.minCoeff(), nonConfVec.maxCoeff(),"Non-Conforming Function", 0);
  viewer.set_field(coexactField,"Rot. Cogradient field", 0, 1);
  viewer.launch();
}
