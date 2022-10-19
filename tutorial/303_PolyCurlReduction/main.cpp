#include <iostream>
#include <fstream>
#include <directional/readOFF.h>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include "tutorial_shared_path.h"

using namespace std;

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawFieldOrig, rawFieldCF;
directional::CartesianField combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge
directional::DirectionalViewer viewer;

double curlMax, curlMaxOrig;
int N;
int iter=0;

// The set of parameters for calculating the curl-free fields
directional::polycurl_reduction_parameters params;

// Solver data (needed for precomputation)
directional::PolyCurlReductionSolverData pcrdata;


typedef enum {ORIGINAL_FIELD, ORIGINAL_CURL, OPTIMIZED_FIELD, OPTIMIZED_CURL} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;


void update_triangle_mesh()
{
  Eigen::VectorXd currCurl = (viewingMode==ORIGINAL_CURL ? curlOrig: curlCF);
    viewer.set_edge_data(currCurl, 0.0,curlMaxOrig);

}

void update_raw_field_mesh()
{
  using namespace std;
  using namespace Eigen;
  
  if ((viewingMode==ORIGINAL_CURL) || (viewingMode==OPTIMIZED_CURL)){
    viewer.toggle_seams(false);
    viewer.toggle_singularities(false);
    viewer.toggle_field(false);
    viewer.toggle_mesh(false);
    viewer.toggle_edge_data(true);
  } else {
    viewer.toggle_seams(true);
    viewer.toggle_singularities(true);
    viewer.toggle_field(true);
    viewer.toggle_mesh(true);
    viewer.toggle_edge_data(false);
    viewer.set_field(viewingMode==ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF,directional::DirectionalViewer::indexed_glyph_colors(viewingMode==ORIGINAL_FIELD ? combedFieldOrig.extField : combedFieldCF.extField));
    viewer.set_seams((viewingMode==ORIGINAL_FIELD ? combedFieldOrig.matching : combedFieldCF.matching));
  }
  
}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  
  if ((key >= '1') && (key <='4'))
    viewingMode = (ViewingModes)(key - '1');
  
  if (key == 'A')
  {
    //do a batch of iterations
    printf("--Improving Curl--\n");
    for (int bi = 0; bi<5; ++bi)
    {
      
      printf("\n\n **** Batch %d ****\n", iter);
      directional::polycurl_reduction_solve(pcrdata, params, rawFieldCF, iter ==0);
      iter++;
      params.wSmooth *= params.redFactor_wsmooth;
    }
    
    Eigen::VectorXi prinIndices;
    directional::curl_matching(rawFieldCF, curlCF);
    directional::combing(rawFieldCF, combedFieldCF);
    directional::curl_matching(combedFieldCF,curlCF);
    curlMax= curlCF.maxCoeff();
    std:: cout<<"curlMax optimized: "<<curlMax<<std::endl;
  }
  
  if (key == 'W'){
    if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-cf.rawfield", rawFieldCF))
      std::cout << "Saved raw field" << std::endl;
    else
      std::cout << "Unable to save raw field. " << std::endl;
  }
  
  update_triangle_mesh();
  update_raw_field_mesh();
  return false;
}

int main(int argc, char *argv[])
{
  
  std::cout <<
  "  A      Optimize 5 batches for curl reduction." << std::endl <<
  "  1      Original field" << std::endl <<
  "  2      L2 norm of original-field curl" << std::endl <<
  "  3      Curl-reduced field" << std::endl <<
  "  4      Curl of curl-reduced field." << std::endl;
  
  // Load a mesh
  directional::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", mesh);
  ftb.init(mesh);
  rawFieldOrig.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
  rawFieldCF.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka.rawfield", ftb, N, rawFieldOrig);
  
  //combing the field in a way that minimizes curl
  directional::curl_matching(rawFieldOrig,curlOrig);
  curlMaxOrig= curlOrig.maxCoeff();
  curlMax = curlMaxOrig;
  std:: cout<<"curlMax original: "<<curlMax<<std::endl;
  
  directional::combing(rawFieldOrig, combedFieldOrig);
  //directional::curl_matching(combedFieldOrig, curlOrig);
  
  //trivial constraints
  Eigen::VectorXi b; b.resize(1); b<<0;
  Eigen::MatrixXd bc; bc.resize(1,6); bc<<rawFieldOrig.extField.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b<<1;
  directional::polycurl_reduction_precompute(mesh, b, bc, blevel, rawFieldOrig , pcrdata);
  
  rawFieldCF = rawFieldOrig;
  combedFieldCF = combedFieldOrig;
  curlCF = curlOrig;
  
  //triangle mesh setup
  viewer.set_mesh(mesh);
  update_triangle_mesh();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  return 0;
}
