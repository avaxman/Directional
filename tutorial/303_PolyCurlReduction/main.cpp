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
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include "tutorial_shared_path.h"

using namespace std;

Eigen::VectorXi matchingOrig, matchingCF, combedMatchingOrig, combedMatchingCF;
Eigen::VectorXd effortOrig, effortCF, combedEffortOrig, combedEffortCF;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V, barycenters;
Eigen::MatrixXd CMesh;
Eigen::MatrixXd rawFieldOrig, rawFieldCF;
Eigen::MatrixXd combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge
directional::DirectionalViewer viewer;

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;


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
    viewer.set_edge_data(currCurl, 0.0,curlMaxOrig, EV, FE, EF);

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
    viewer.set_field(viewingMode==ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF,directional::DirectionalViewer::indexed_glyph_colors(viewingMode==ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF));
    viewer.set_singularities((viewingMode==ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode==ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF));
    viewer.set_seams(EV, FE, EF, (viewingMode==ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF));
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
    directional::curl_matching(V, F,EV, EF, FE, rawFieldCF, matchingCF, effortCF, curlCF,singVerticesCF, singIndicesCF);
    directional::combing(V,F, EV, EF, FE, rawFieldCF, matchingCF, combedFieldCF);
    directional::curl_matching(V, F,EV, EF, FE, combedFieldCF, combedMatchingCF, combedEffortCF, curlCF,singVerticesCF, singIndicesCF);
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
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka.rawfield", N, rawFieldOrig);
  igl::edge_topology(V, F, EV, FE, EF);
  
  igl::barycenter(V, F, barycenters);
  
  Eigen::VectorXi prinIndices;
  directional::curl_matching(V, F,EV, EF, FE, rawFieldOrig, matchingOrig, effortOrig, curlOrig, singVerticesOrig, singIndicesOrig);
  curlMaxOrig= curlOrig.maxCoeff();
  curlMax = curlMaxOrig;
  std:: cout<<"curlMax original: "<<curlMax<<std::endl;
  
  directional::combing(V,F, EV, EF, FE,rawFieldOrig, matchingOrig,combedFieldOrig);
  directional::curl_matching(V, F,EV, EF, FE, combedFieldOrig, combedMatchingOrig, combedEffortOrig, curlOrig, singVerticesOrig, singIndicesOrig);
  
  //trivial constraints
  Eigen::VectorXi b; b.resize(1); b<<0;
  Eigen::MatrixXd bc; bc.resize(1,6); bc<<rawFieldOrig.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b<<1;
  directional::polycurl_reduction_precompute(V, F, b, bc, blevel, rawFieldOrig , pcrdata);
  
  rawFieldCF = rawFieldOrig;
  matchingCF = matchingOrig;
  effortCF = effortOrig;
  combedFieldCF =combedFieldOrig;
  combedMatchingCF =combedMatchingOrig;
  combedEffortCF =combedEffortOrig;
  curlCF = curlOrig;
  singVerticesCF = singVerticesOrig;
  singIndicesCF = singIndicesOrig;
  
 
  //triangle mesh setup
  viewer.set_mesh(V, F);
  update_triangle_mesh();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  return 0;
}
