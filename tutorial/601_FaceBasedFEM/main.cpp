#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/ScalarFunction2D.h>
#include <directional/PLFunction2D.h>
#include <directional/VolumeForm2D.h>
#include <directional/DiamondForm2D.h>
#include <directional/directional_viewer.h>
#include <directional/CochainComplex.h>


directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField exactField, coexactField;
directional::DirectionalViewer viewer;
directional::PLFunction2D<double> PLScalarFunc;
directional::DiamondForm2D<double> DiamondForm;
directional::CochainComplex<double> PLComplex;
directional::CochainComplex<double> NonConfPLComplex;

template<typename NumberType>
directional::CochainComplex<NumberType> vector_field_complex_2D(directional::ScalarFunction2D<NumberType>& scalarFunc, directional::CartesianField& field, directional::VolumeForm2D<NumberType>& form){
    directional::CochainComplex<NumberType> complex;
    complex.differentials.resize(2);
    complex.metrics.resize(3);
    complex.invMetrics.resize(3);
    complex.differentials[0] = scalarFunc.gradient_matrix();
    complex.differentials[1] = field.curl_matrix();
    complex.metrics[0] = scalarFunc.mass_matrix();
    complex.metrics[1] = field.mass_matrix();
    complex.metrics[2] = form.mass_matrix();
    complex.invMetrics[0] = scalarFunc.inv_mass_matrix();
    complex.invMetrics[1] = field.inv_mass_matrix();
    complex.invMetrics[2] = form.inv_mass_matrix();
    return complex;
}

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

  //Exact (Gradient field):
  PLScalarFunc.init(mesh, confVec);
  exactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
  PLScalarFunc.gradient(exactField);
  DiamondForm.init(mesh, Eigen::VectorXd::Zero(mesh.EV.rows()), 1);

  coexactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
  PLComplex = vector_field_complex_2D(PLScalarFunc, exactField, DiamondForm);

  //Generating the dual complex (representing PL non-conforming functions -> PC vector field -> Voronoi forms
  NonConfPLComplex = PLComplex.dual_complex();
  coexactField.set_extrinsic_field(NonConfPLComplex.invMetrics[1]*NonConfPLComplex.differentials[0]*nonConfVec); //Emulating rotated cogradient field of the PL non-conforming function

  //computing divergence
  Eigen::VectorXd divergence = PLScalarFunc.gradient_matrix().transpose()* coexactField.mass_matrix()*coexactField.flatten();

  //demonstrating the exact sequences
  std::cout<<"max abs curl of exact field (should be numerically zero): "<<exactField.curl().cwiseAbs().maxCoeff()<<std::endl;
  std::cout<<"max abs divergence of coexact field (should be numerically zero): "<<divergence.cwiseAbs().maxCoeff()<<std::endl;

  //triangle mesh setup
  viewer.init();
  viewer.set_mesh(mesh);
  viewer.set_vertex_data(confVec, confVec.minCoeff(), confVec.maxCoeff(),"Conforming Function", 0);
  viewer.set_field(exactField,"Gradient field", 0, 0);
  viewer.set_edge_data(nonConfVec, nonConfVec.minCoeff(), nonConfVec.maxCoeff(),"Non-Conforming Function", 0);
  viewer.set_field(coexactField,"Rot. Cogradient field", 0, 1);
  viewer.launch();
}
