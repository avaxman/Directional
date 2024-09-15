#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/PLFunction.h>
#include <directional/NonConfPLFunction.h>
#include <directional/DiamondForm.h>
#include <directional/VoronoiForm.h>
#include <directional/directional_viewer.h>


directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField exactField, coexactField, harmField;
Eigen::MatrixXd harmBasis;
directional::DirectionalViewer viewer;
directional::PLFunction confScalarField;
directional::NonConfPLFunction nonConfScalarField;

int main()
{

  directional::readOFF(TUTORIAL_DATA_PATH "/Input_wheel.off",mesh);
  ftb.init(mesh);

  Eigen::VectorXd confVec(mesh.dcel.vertices.size()), nonConfVec(mesh.dcel.edges.size());
  for (int i=0;i<mesh.dcel.vertices.size();i++)
      confVec[i] = 10.0*sin(mesh.V(i,0))*cos(mesh.V(i,1));
  for (int i=0;i<mesh.dcel.edges.size();i++){
      Eigen::RowVector3d midEdgePoint = 0.5*(mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1)));
      confVec[i] = 10.0*sin(midEdgePoint(i,1))*cos(midEdgePoint(i,2));
  }
  confScalarField.from_vector(confVec);
  nonConfScalarField.from_vector(nonConfVec);

  confScalarField.gradient(exactField);
  nonConfScalarField.gradient(coexactField);
  nonConfScalarField.dualize();

  //demonstrating the exact sequences
  std::cout<<"max abs curl of exact field: "<<exactField.curl().cwiseAbs().maxCoeff()<<std::endl;
  std::cout<<"max abs divergence of coexact field: "<<exactField.div().cwiseAbs().maxCoeff()<<std::endl;

  //triangle mesh setup
  viewer.init();
  viewer.set_mesh(mesh);
  viewer.set_vertex_data(confVec, confVec.minCoeff(), confVec.maxCoeff(),"Conforming Function", 0);
  viewer.set_field(confScalarField,"Gradient field", 0, 0);
  viewer.set_edge_data(nonConfVec, nonConfVec.minCoeff(), nonConfVec.maxCoeff(),"Non-Conforming Function", 0);
  viewer.set_field(coexactField,"Rot. Cogradient field", 0, 1);
  viewer.launch();
}
