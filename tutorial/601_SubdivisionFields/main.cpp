#include <iostream>
#include <fstream>
#include <unordered_set>
#include <igl/readOBJ.h>
#include <directional/readOBJ.h>
#include <directional/writeOBJ.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/subdivide_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/cut_mesh_with_singularities.h>
#include "tutorial_shared_path.h"



using namespace std;

directional::TriMesh meshCoarse, meshCutCoarse, meshFine, meshCutFine;
directional::IntrinsicFaceTangentBundle ftbCoarse, ftbFine;
directional::CartesianField rawFieldCoarse, rawFieldFine, combedFieldCoarse, combedFieldFine, powerField;
Eigen::VectorXd curlCoarse, curlFine; // norm of curl per edge
Eigen::MatrixXd cutReducedUVCoarse, cutFullUVCoarse, cornerWholeUVCoarse;
Eigen::MatrixXd cutReducedUVFine, cutFullUVFine, cornerWholeUVFine;

// The igl viewer
directional::DirectionalViewer viewer;

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
      directional::writeOBJ(TUTORIAL_SHARED_PATH "/knight-coarse-param.obj", meshCutCoarse,  cutReducedUVCoarse, meshCutCoarse.F);
      directional::writeOBJ(TUTORIAL_SHARED_PATH "/knight-fine-param.obj", meshCutFine,cutReducedUVFine, meshCutFine.F);
      break;
  }
  update_viewer();
  return true;
}

// Parameterize a mesh from a field
void parameterize_mesh(const directional::TriMesh& wholeMesh,
                       const directional::CartesianField& rawField,
                       directional::TriMesh& cutMesh,
                       Eigen::MatrixXd& cutFullUV)
{
  directional::CartesianField combedField;

  directional::IntegrationData intData(N);
  std::cout << "Setting up Integration" << std::endl;

  intData.verbose = false;
  intData.integralSeamless = false;
  intData.lengthRatio = 0.03;

  directional::setup_integration(rawField,intData, cutMesh, combedField);

  std::cout << "Integrating" << std::endl;
  Eigen::MatrixXd cornerWholeUV;
  directional::integrate(combedField, intData, cutMesh, cutFullUV, cornerWholeUV);
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
  
  directional::readOBJ(TUTORIAL_SHARED_PATH "/bunny1k.obj", meshCoarse);
  ftbCoarse.init(meshCoarse);
 
  //Reducing curl from coarse mesh
  Eigen::VectorXi b;
  Eigen::MatrixXd bc;
  b.resize(0);
  bc.resize(0, 3);

  cout<<"Computing initial field"<<endl;
  powerField.init(ftbCoarse, directional::fieldTypeEnum::POWER_FIELD, N);
  directional::power_field(ftbCoarse, b, bc, Eigen::VectorXd::Constant(b.size(),-1),N, powerField);
  directional::power_to_raw(powerField, N, rawFieldCoarse, true);
   
  // Precompute polycurl reduction data.
  b.resize(1); b << 0;
  bc.resize(1, 6); bc << rawFieldCoarse.extField.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b << 1;
  directional::polycurl_reduction_precompute(meshCoarse, b, bc, blevel, rawFieldCoarse, pcrdata);
  
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
  
  directional::curl_matching(rawFieldCoarse,  curlCoarse);
  directional::combing(rawFieldCoarse,combedFieldCoarse);
  //directional::curl_matching(combedFieldCoarse, curlCoarse);
  double curlMax = curlCoarse.maxCoeff();
  std::cout << "Coarse optimized absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing coarse parameterization"<<endl;
  parameterize_mesh(meshCoarse, rawFieldCoarse,meshCutCoarse, cutFullUVCoarse);
  
  cout<<"Doing field subdivision"<<endl;
  //Subdividing field and mesh. Curl precision should result in very low curl on the fine mesh as well
  directional::subdivide_field(rawFieldCoarse, targetLevel,  meshFine, ftbFine, rawFieldFine);

  directional::curl_matching(rawFieldFine, curlFine);
  directional::combing(rawFieldFine,combedFieldFine);
  directional::curl_matching(combedFieldFine, curlFine);
  curlMax = curlFine.maxCoeff();
  std::cout << "Fine subdivided absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing fine parameterization"<<endl;
  parameterize_mesh(meshFine, rawFieldFine, meshCutFine, cutFullUVFine);
  
  cout<<"Done!"<<endl;
  
  //coarse mesh
  viewer.set_mesh(meshCoarse, 0);
  viewer.set_field(combedFieldCoarse,Eigen::MatrixXd(),0);
  viewer.set_seams(combedFieldCoarse.matching,0);
  viewer.set_edge_data(curlCoarse, curlCoarse.minCoeff(), curlCoarse.maxCoeff(),0);
  
  //fine mesh
  viewer.set_mesh(meshFine, 1);
  viewer.set_field(combedFieldFine,Eigen::MatrixXd(),1);
  viewer.set_seams(combedFieldFine.matching,1);
  viewer.set_edge_data(curlFine, curlCoarse.minCoeff(), curlCoarse.maxCoeff(), 1);
  
  //coarse texture mesh
  viewer.set_mesh(meshCutCoarse,2);
  viewer.set_uv(cutFullUVCoarse,2);
  viewer.toggle_texture(true,2);
  
  //coarse texture mesh
  viewer.set_mesh(meshCutFine, 3);
  viewer.set_uv(cutFullUVFine,3);
  viewer.toggle_texture(true,3);
  
  // Update view
  update_viewer();
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  return 0;
}
