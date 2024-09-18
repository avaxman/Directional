#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/PLFunction.h>
#include <directional/CochainComplex.h.h>
#include <directional/DiamondForm.h>
#include <directional/VoronoiForm.h>
#include <directional/directional_viewer.h>
#include <directional/ccochain_complex.h>


directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField exactField, coexactField;
directional::DirectionalViewer viewer;
directional::PLFunction<double> PLScalarFunc;
directional::CochainComplex<double> PLComplex;
directional::CochainComplex<double> NonConfPLComplex;

int main()
{

  directional::readOFF(TUTORIAL_DATA_PATH "/Input_wheel.off",mesh);
  ftb.init(mesh);

  Eigen::VectorXd confVec(mesh.dcel.vertices.size()), nonConfVec(mesh.dcel.edges.size());
  for (int i=0;i<mesh.dcel.vertices.size();i++)
      confVec[i] = 10.0*sin(mesh.V(i,0))*cos(mesh.V(i,1));
  for (int i=0;i<mesh.dcel.edges.size();i++){
      Eigen::RowVector3d midEdgePoint = 0.5*(mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1)));
      nonConfVec[i] = 10.0*sin(midEdgePoint(i,1))*cos(midEdgePoint(i,2));
  }

  PLScalarFunc.init(mesh, confVec);
  exactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
  coexactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
  PLComplex.differentials.resize(2);
  PLComplex.metrics.resize(3);
  PLComplex.differentials[0]=PLScalarFunc.gradient_matrix();
  PLComplex.differentials[1]=exactField.curl_matrix();
  PLComplex.metrics[0] = PLScalarFunc.mass_matrix();
  PLComplex.metrics[1] = exactField.mass_matrix();
  PLComplex.metrics[2] = PLDiamondForm.mass_matrix();

  //Generating the dual complex (representing PL non-conforming functions -> PC vector field -> Voronoi forms
  NonConfPLComplex = PLComplex.dual_complex();
  Eigen::VectorXd coexactVec = NonConfPLComplex.differentials[0]*nonConfVec;  //automaticaly get rotated cogradient field of the PL non-conforming function
  coexactField.from_vector(coexactVec);

  //demonstrating the exact sequences
  std::cout<<"max abs curl of exact field (should be numerically zero): "<<exactField.curl().cwiseAbs().maxCoeff()<<std::endl;
  std::cout<<"max abs divergence of coexact field (should be numerically zero): "<<coexactField.div().cwiseAbs().maxCoeff()<<std::endl;

  //triangle mesh setup
  viewer.init();
  viewer.set_mesh(mesh);
  viewer.set_vertex_data(confVec, confVec.minCoeff(), confVec.maxCoeff(),"Conforming Function", 0);
  viewer.set_field(PLScalarFunc,"Gradient field", 0, 0);
  viewer.set_edge_data(nonConfVec, nonConfVec.minCoeff(), nonConfVec.maxCoeff(),"Non-Conforming Function", 0);
  viewer.set_field(coexactField,"Rot. Cogradient field", 0, 1);
  viewer.launch();
}
