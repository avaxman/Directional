#include <iostream>
#include <fstream>
#include <unordered_set>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/edge_topology.h>
#include <igl/jet.h>
#include <igl/false_barycentric_subdivision.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <igl/edge_flaps.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/subdivide_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include "tutorial_shared_path.h"
#include <directional/cut_mesh_with_singularities.h>
#include <directional/directional_viewer.h>


using namespace std;

Eigen::VectorXi matchingCoarse, matchingFine, combedMatchingCoarse, combedMatchingFine, singVerticesCoarse, singVerticesFine, singIndicesCoarse, singIndicesFine;
Eigen::VectorXd effortCoarse, effortFine, combedEffortCoarse, combedEffortFine;
Eigen::MatrixXi FCoarse, FFine, FCutCoarse, FCutFine;
Eigen::MatrixXi EVCoarse, EFCoarse, FECoarse,EVFine, EFFine, FEFine;

// Separate matrices for different visualizations
Eigen::MatrixXd VCoarse, VCutCoarse, VFine, VCutFine;
Eigen::MatrixXd CMeshCoarse, CMeshFine;
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


void update_view()
{
  viewer.set_active(viewingMode==COARSE_FIELD | viewingMode==COARSE_CURL,0);
  viewer.set_active(viewingMode==COARSE_PARAMETERIZATION,1);
  viewer.set_active(viewingMode==FINE_FIELD | viewingMode==FINE_CURL,2);
  viewer.set_active(viewingMode==FINE_PARAMETERIZATION,3);
  
  viewer.toggle_field(viewingMode==COARSE_FIELD,0);
  viewer.toggle_field(viewingMode==FINE_FIELD,2);
  viewer.toggle_singularities(viewingMode==COARSE_FIELD,0);
  viewer.toggle_singularities(viewingMode==FINE_FIELD,2);
  viewer.toggle_seams(viewingMode==COARSE_FIELD,0);
  viewer.toggle_seams(viewingMode==FINE_FIELD,2);
  
  if (viewingMode==COARSE_CURL){
    Eigen::VectorXd faceCurl = AE2FCoarse*curlCoarse;
    igl::jet(faceCurl, 0.0,0.01, CMeshCoarse);
    viewer.set_mesh_colors(CMeshCoarse,0);
  } else {
    viewer.set_mesh_colors(Eigen::MatrixXd(),0);
  }
  
  if (viewingMode==FINE_CURL){
    Eigen::VectorXd faceCurl = AE2FFine*curlFine;
    igl::jet(faceCurl, 0.0,0.01, CMeshFine);
    viewer.set_mesh_colors(CMeshFine,2);
  } else {
      viewer.set_mesh_colors(Eigen::MatrixXd(),2);
  }
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
  update_view();
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
  directional::ParameterizationData pd;
  directional::cut_mesh_with_singularities(VWhole, FWhole, singVertices, pd.face2cut);
  directional::combing(VWhole, FWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField, combedMatching);
  
  std::cout << "Setting up parameterization" << std::endl;
  
  directional::setup_parameterization(directional::sign_symmetry(N),VWhole, FWhole, EV, EF, FE, combedMatching, singVertices, pd, VCut, FCut);
  
  double lengthRatio = 0.03;
  bool isInteger = false;  //do not do translational seamless.
  std::cout << "Solving parameterization" << std::endl;
  Eigen::MatrixXd cutReducedUV,  cornerWholeUV;
  directional::parameterize(VWhole, FWhole, FE, combedField, lengthRatio, pd, VCut, FCut, isInteger, cutReducedUV, cutFullUV, cornerWholeUV);
  cutFullUV = cutFullUV.block(0, 0, cutFullUV.rows(), 2).eval();
  std::cout << "Done!" << std::endl;
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
  
  //directional::subdivide_field(VCoarse, FCoarse, EVCoarse, EFCoarse, combedFieldCoarse, combedMatchingCoarse, targetLevel, VFine, FFine, EVFine, EFFine, rawFieldFine, matchingFine);

  igl::edge_topology(VFine, FFine, EVFine, FEFine, EFFine);
    
  directional::curl_matching(VFine, FFine, EVFine, EFFine, FEFine, rawFieldFine, matchingFine, effortFine, curlFine, singVerticesFine, singIndicesFine);
  directional::combing(VFine, FFine, EVFine, EFFine, FEFine, rawFieldFine, matchingFine, combedFieldFine);
  directional::curl_matching(VFine, FFine, EVFine, EFFine, FEFine, combedFieldFine, combedMatchingFine, combedEffortFine, curlFine,singVerticesFine, singIndicesFine);
  curlMax = curlFine.maxCoeff();
  std::cout << "Fine subdivided absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing fine parameterization"<<endl;
  parameterize_mesh(VFine, FFine, EVFine, EFFine, FEFine, singVerticesFine, rawFieldFine, matchingFine, VCutFine, FCutFine, cutFullUVFine);
  
  cout<<"Done!"<<endl;
  
  //coarse Whole mesh
  viewer.set_mesh(VCoarse, FCoarse,Eigen::MatrixXd(),0);
  viewer.set_field(combedFieldCoarse,directional::indexed_glyph_colors(combedFieldCoarse),0);
  viewer.set_singularities(N, singVerticesCoarse, singIndicesCoarse,0);
  viewer.set_seams(EVCoarse,combedMatchingCoarse,0);
  viewer.toggle_mesh_edges(false,0);
  
  //Coarse cut mesh
  viewer.set_mesh(VCutCoarse, FCutCoarse,Eigen::MatrixXd(),1);
  viewer.set_uv(cutFullUVCoarse,1);
  viewer.set_texture(directional::DirectionalViewer::default_texture(),1);
  viewer.toggle_mesh_edges(false,1);
  
  //Fine whole mesh
  viewer.set_mesh(VFine, FFine,Eigen::MatrixXd(),2);
  viewer.set_field(combedFieldFine,directional::indexed_glyph_colors(combedFieldFine),2);
  viewer.set_singularities(N, singVerticesFine, singIndicesFine,2);
  viewer.set_seams(EVFine,combedMatchingFine,2);
  viewer.toggle_mesh_edges(false,2);
  
  //Fine cut mesh
  viewer.set_mesh(VCutFine, FCutFine,Eigen::MatrixXd(),3);
  viewer.set_uv(cutFullUVFine,3);
  viewer.set_texture(directional::DirectionalViewer::default_texture(),3);
  viewer.toggle_mesh_edges(false,3);
    
  update_view();
  
  //creating the AE2F operator for curl viewing on faces
  std::vector<Eigen::Triplet<double> > AE2FTriplets;
  for (int i = 0; i < EFCoarse.rows(); i++) {
    AE2FTriplets.push_back(Eigen::Triplet<double>(EFCoarse(i, 0), i, 1.0));
    AE2FTriplets.push_back(Eigen::Triplet<double>(EFCoarse(i, 1), i, 1.0));
  }
  AE2FCoarse.resize(FCoarse.rows(), EFCoarse.rows());
  AE2FCoarse.setFromTriplets(AE2FTriplets.begin(), AE2FTriplets.end());
  
  AE2FTriplets.clear();
  for (int i = 0; i < EFFine.rows(); i++) {
    AE2FTriplets.push_back(Eigen::Triplet<double>(EFFine(i, 0), i, 1.0));
    AE2FTriplets.push_back(Eigen::Triplet<double>(EFFine(i, 1), i, 1.0));
  }
  AE2FFine.resize(FFine.rows(), EFFine.rows());
  AE2FFine.setFromTriplets(AE2FTriplets.begin(), AE2FTriplets.end());
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  return 0;
}
