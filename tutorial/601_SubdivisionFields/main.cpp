#include <iostream>
#include <fstream>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/local_basis.h>
#include <igl/avg_edge_length.h>
#include <igl/is_border_vertex.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/false_barycentric_subdivision.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <igl/edge_flaps.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/subdivide_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include "tutorial_shared_path.h"
#include <directional/cut_mesh_with_singularities.h>
#include <unordered_set>


using namespace std;

Eigen::VectorXi matchingCoarse, matchingFine, combedMatchingCoarse, combedMatchingFine, singVerticesCoarse, singVerticesFine, singIndicesCoarse, singIndicesFine;
Eigen::VectorXd effortCoarse, effortFine, combedEffortCoarse, combedEffortFine;
Eigen::MatrixXi FCoarse, FFine, FCutCoarse, FCutFine;
Eigen::MatrixXi EVCoarse, EFCoarse, FECoarse,EVFine, EFFine, FEFine;

Eigen::MatrixXd VCoarse, VCutCoarse, VFine, VCutFine;
Eigen::MatrixXd CMeshCoarse, CFieldCoarse, CSingsCoarse, CSeamsCoarse, CMeshFine, CFieldFine, CSingsFine, CSeamsFine;
Eigen::MatrixXd rawFieldCoarse, rawFieldFine;
Eigen::MatrixXd combedFieldCoarse, combedFieldFine;
Eigen::VectorXd curlCoarse, curlFine; // norm of curl per edge

Eigen::MatrixXd cutReducedUVCoarse, cutFullUVCoarse, cornerWholeUVCoarse;
Eigen::MatrixXd cutReducedUVFine, cutFullUVFine, cornerWholeUVFine;

Eigen::MatrixXcd powerField;


// The igl viewer
directional::DirectionalViewer viewer;

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;

//for averaging curl to faces, for visualization
Eigen::SparseMatrix<double> AE2FCoarse, AE2FFine;

int N = 4;
int targetLevel = 2;
int iter = 0;

// The set of parameters for calculating the curl-free fields
directional::polycurl_reduction_parameters params;

// Solver data (needed for precomputation)
directional::PolyCurlReductionSolverData pcrdata;

typedef enum {COARSE_FIELD, COARSE_CURL, COARSE_PARAMETERIZATION,FINE_FIELD, FINE_CURL, FINE_PARAMETERIZATION} ViewingModes;
ViewingModes viewingMode=COARSE_FIELD;


void update_viewer()
{
  viewer.toggle_mesh(viewingMode==COARSE_FIELD, 0);
  viewer.toggle_mesh(viewingMode==FINE_FIELD, 1);
  viewer.toggle_mesh(viewingMode==COARSE_PARAMETERIZATION, 2);
  viewer.toggle_mesh(viewingMode==FINE_PARAMETERIZATION, 3);
  viewer.toggle_field(viewingMode==COARSE_FIELD,0);
  viewer.toggle_field(viewingMode==FINE_FIELD,1);
  viewer.toggle_seams(viewingMode==COARSE_FIELD,0);
  viewer.toggle_seams(viewingMode==FINE_FIELD,1);
  viewer.toggle_edge_data(viewingMode==COARSE_CURL,0);
  viewer.toggle_edge_data(viewingMode==FINE_CURL,1);
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = COARSE_FIELD; break;
    case '2': viewingMode = COARSE_CURL; break;
    case '3': viewingMode = COARSE_PARAMETERIZATION; break;
    case '4': viewingMode = FINE_FIELD; break;
    case '5': viewingMode = FINE_CURL; break;
    case '6': viewingMode = FINE_PARAMETERIZATION; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/knight-coarse-param.obj", VCutCoarse, FCutCoarse, emptyMat, emptyMat, cutReducedUVCoarse, FCutCoarse);
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/knight-fine-param.obj", VCutFine, FCutFine, emptyMat, emptyMat, cutReducedUVFine, FCutFine);
      break;
  }
  update_viewer();
  return true;
}

// Parameterize a mesh from a field
void parameterize_mesh(const Eigen::MatrixXd& VWhole,
                       const Eigen::MatrixXi& FWhole,
                       const Eigen::MatrixXi& EV,
                       const Eigen::MatrixXi& EF,
                       const Eigen::MatrixXi& FE,
                       const Eigen::VectorXi& singVertices,
                       const Eigen::MatrixXd& rawField,
                       const Eigen::VectorXi& matching,
                       Eigen::MatrixXd& VCut,
                       Eigen::MatrixXi& FCut,
                       Eigen::MatrixXd& cutFullUV)
{
  Eigen::MatrixXd combedField;
  Eigen::VectorXi combedMatching;

  directional::IntegrationData intData(N);
  std::cout << "Setting up Integration" << std::endl;

  intData.verbose = false;
  intData.integralSeamless = false;
  intData.lengthRatio = 0.03;

  directional::setup_integration(VWhole, FWhole, EV, EF, FE, rawField, matching, singVertices, intData, VCut, FCut, combedField, combedMatching);

  std::cout << "Integrating" << std::endl;
  Eigen::MatrixXd cornerWholeUV;
  directional::integrate(VWhole, FWhole, FE, combedField, intData, VCut, FCut, cutFullUV, cornerWholeUV);
  std::cout << "Done!" << std::endl;

  cutFullUV = cutFullUV.block(0, 0, cutFullUV.rows(), 2).eval();
}


int main(int argc, char *argv[])
{
  auto keyAction = [](const std::string& key, const std::string& description)
  {
    std::cout << "  " << key << "      " << description << std::endl;
  };
  //keyAction("A", "Optimize 60 batches for curl reduction.");
  keyAction("1", "Show coarse field.");
  keyAction("2", "Show curl of coarse field.");
  keyAction("3", "Show coarse parameterization.");
  keyAction("4", "Show fine field.");
  keyAction("5", "Show curl of fine field.");
  keyAction("6", "Show fine parameterization.");
  
  igl::readOBJ(TUTORIAL_SHARED_PATH "/bunny1k.obj", VCoarse, FCoarse);
 
  igl::edge_topology(VCoarse, FCoarse, EVCoarse, FECoarse, EFCoarse);

  //Reducing curl from coarse mesh
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  b.resize(0);
  bc.resize(0, 3);

  cout<<"Computing initial field"<<endl;
  directional::power_field(VCoarse, FCoarse, b, bc, N, powerField);
  powerField.array()/=powerField.cwiseAbs().array();
  directional::power_to_raw(VCoarse,FCoarse,powerField, N, rawFieldCoarse);
   
  // Precompute polycurl reduction data.
  b.resize(1); b << 0;
  bc.resize(1, 6); bc << rawFieldCoarse.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b << 1;
  directional::polycurl_reduction_precompute(VCoarse, FCoarse, b, bc, blevel, rawFieldCoarse, pcrdata);
  
  printf("--Improving Curl--\n");
  //Performing 60 Batches of curl reduction
  for (int bi = 0; bi < 30; ++bi)
  {
    
    printf("\n\n **** Batch %d ****\n", iter);
    // Solve
    directional::polycurl_reduction_solve(pcrdata, params, rawFieldCoarse, iter == 0);
    ++iter;
    // Adjust the smoothness weight
    params.wSmooth *= params.redFactor_wsmooth;
  }
  
  directional::curl_matching(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, rawFieldCoarse, matchingCoarse, effortCoarse, curlCoarse, singVerticesCoarse, singIndicesCoarse);
  directional::combing(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, rawFieldCoarse, matchingCoarse, combedFieldCoarse);
  directional::curl_matching(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, combedFieldCoarse, combedMatchingCoarse, combedEffortCoarse, curlCoarse, singVerticesCoarse, singIndicesCoarse);
  double curlMax = curlCoarse.maxCoeff();
  std::cout << "Coarse optimized absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing coarse parameterization"<<endl;
  parameterize_mesh(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, singVerticesCoarse, rawFieldCoarse, matchingCoarse, VCutCoarse, FCutCoarse, cutFullUVCoarse);
  
  cout<<"Doing field subdivision"<<endl;
  //Subdividing field and mesh. Curl precision should result in very low curl on the fine mesh as well
  directional::subdivide_field(VCoarse, FCoarse, rawFieldCoarse, targetLevel,  VFine, FFine, rawFieldFine);
  
  igl::edge_topology(VFine, FFine, EVFine, FEFine, EFFine);
  directional::curl_matching(VFine, FFine, EVFine, EFFine, FEFine, rawFieldFine, matchingFine, effortFine, curlFine, singVerticesFine, singIndicesFine);
  directional::combing(VFine, FFine, EVFine, EFFine, FEFine, rawFieldFine, matchingFine, combedFieldFine);
  directional::curl_matching(VFine, FFine, EVFine, EFFine, FEFine, combedFieldFine, combedMatchingFine, combedEffortFine, curlFine, singVerticesCoarse, singIndicesCoarse);
  curlMax = curlFine.maxCoeff();
  std::cout << "Fine subdivided absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing fine parameterization"<<endl;
  parameterize_mesh(VFine, FFine, EVFine, EFFine, FEFine, singVerticesFine, rawFieldFine, matchingFine, VCutFine, FCutFine, cutFullUVFine);
  
  cout<<"Done!"<<endl;
  
  //coarse mesh
  viewer.set_mesh(VCoarse, FCoarse, Eigen::MatrixXd(), 0);
  viewer.set_field(combedFieldCoarse,Eigen::MatrixXd(),0);
  viewer.set_singularities(singVerticesCoarse,singIndicesCoarse,0);
  viewer.set_seams(EVCoarse, FECoarse, EFCoarse, combedMatchingCoarse,0);
  viewer.set_edge_data(curlCoarse, curlCoarse.minCoeff(), curlCoarse.maxCoeff(), EVCoarse, FECoarse, EFCoarse,0);
  
  //fine mesh
  viewer.set_mesh(VFine, FFine, Eigen::MatrixXd(), 1);
  viewer.set_field(combedFieldFine,Eigen::MatrixXd(),1);
  viewer.set_singularities(singVerticesFine,singIndicesFine,1);
  viewer.set_seams(EVFine, FEFine, EFFine, combedMatchingFine,1);
  viewer.set_edge_data(curlFine, curlCoarse.minCoeff(), curlCoarse.maxCoeff(), EVFine, FEFine, EFFine,1);
  
  //coarse texture mesh
  viewer.set_mesh(VCutCoarse, FCutCoarse, Eigen::MatrixXd(), 2);
  viewer.set_uv(cutFullUVCoarse,2);
  viewer.toggle_texture(true,2);
  
  //coarse texture mesh
  viewer.set_mesh(VCutFine, FCutCoarse, Eigen::MatrixXd(), 3);
  viewer.set_uv(cutFullUVFine,3);
  viewer.toggle_texture(true,3);
  
  // Update view
  update_viewer();
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  return 0;
}
