#include <iostream>
#include <fstream>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
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
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <igl/edge_flaps.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/subdivide_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include "tutorial_shared_path.h"
#include <directional/cut_mesh_with_singularities.h>
#include <unordered_set>

using namespace std;

Eigen::VectorXi matchingCoarse, matchingFine, combedMatchingCoarse, combedMatchingFine, singVerticesCoarse, singVerticesFine, singIndicesCoarse, singIndicesFine;
Eigen::VectorXd effortCoarse, effortFine, combedEffortCoarse, combedEffortFine;
Eigen::MatrixXi FCoarse, FFine, FCutCoarse, FCutFine, FFieldCoarse, FSingsCoarse, FSeamsCoarse, FFieldFine, FSingsFine, FSeamsFine;
Eigen::MatrixXi EVCoarse, EFCoarse, FECoarse,EVFine, EFFine, FEFine;

// Separate matrices for different visualizations
Eigen::MatrixXd VCoarse, VCutCoarse, VFieldCoarse, VSingsCoarse, VSeamsCoarse, VFine, VCutFine, VFieldFine, VSingsFine, VSeamsFine, barycenters;
Eigen::MatrixXd CMeshCoarse, CFieldCoarse, CSingsCoarse, CSeamsCoarse, CMeshFine, CFieldFine, CSingsFine, CSeamsFine;
Eigen::MatrixXd rawFieldCoarse, rawFieldFine;
Eigen::MatrixXd combedFieldCoarse, combedFieldFine;
Eigen::VectorXd curlCoarse, curlFine; // norm of curl per edge

Eigen::MatrixXd cutReducedUVCoarse, cutFullUVCoarse, cornerWholeUVCoarse;
Eigen::MatrixXd cutReducedUVFine, cutFullUVFine, cornerWholeUVFine;

Eigen::MatrixXcd powerField;


// The igl viewer
igl::opengl::glfw::Viewer viewer;

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
  viewer.data_list[0].clear();
  if (viewingMode==COARSE_FIELD){
    viewer.data_list[0].set_mesh(VCoarse, FCoarse);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].show_texture=false;
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VFieldCoarse, FFieldCoarse);
    viewer.data_list[1].set_colors(CFieldCoarse);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSingsCoarse, FSingsCoarse);
    viewer.data_list[2].set_colors(CSingsCoarse);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    viewer.data_list[3].clear();
    viewer.data_list[3].set_mesh(VSeamsCoarse, FSeamsCoarse);
    viewer.data_list[3].set_colors(CSeamsCoarse);
    viewer.data_list[3].show_faces = true;
    viewer.data_list[3].show_lines = false;
    
  }
  
  if (viewingMode==COARSE_CURL){
    viewer.data_list[0].set_mesh(VCoarse, FCoarse);
    Eigen::VectorXd faceCurl = AE2FCoarse*curlCoarse;
    igl::jet(faceCurl, 0.0,0.01, CMeshCoarse);
    viewer.data_list[0].set_colors(CMeshCoarse);
    viewer.data_list[0].show_texture=false;
  }
  
  if (viewingMode==COARSE_PARAMETERIZATION){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VCutCoarse, FCutCoarse);
    viewer.data_list[0].set_uv(cutFullUVCoarse);
    viewer.data_list[0].show_texture=true;
  }
  
  if (viewingMode==FINE_FIELD){
    viewer.data_list[0].set_mesh(VFine, FFine);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].show_texture=false;
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VFieldFine, FFieldFine);
    viewer.data_list[1].set_colors(CFieldFine);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSingsFine, FSingsFine);
    viewer.data_list[2].set_colors(CSingsFine);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    viewer.data_list[3].clear();
    viewer.data_list[3].set_mesh(VSeamsFine, FSeamsFine);
    viewer.data_list[3].set_colors(CSeamsFine);
    viewer.data_list[3].show_faces = true;
    viewer.data_list[3].show_lines = false;
  }
  
  if (viewingMode==FINE_CURL){
    viewer.data_list[0].set_mesh(VFine, FFine);
    Eigen::VectorXd faceCurl = AE2FFine*curlFine;
    igl::jet(faceCurl, 0.0,0.01, CMeshFine);
    viewer.data_list[0].set_colors(CMeshFine);
    viewer.data_list[0].show_texture=false;
  }
  
  if (viewingMode==FINE_PARAMETERIZATION){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VCutFine, FCutFine);
    viewer.data_list[0].set_uv(cutFullUVFine);
    viewer.data_list[0].show_texture=true;
  }
  
  for (int i=1;i<=3;i++)  //hide all other meshes
    viewer.data_list[i].show_faces= (viewingMode==COARSE_FIELD)||(viewingMode==FINE_FIELD);
  
  
  viewer.data_list[0].show_lines=false;
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
  //directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.rawfield", N, rawFieldCoarse);
 
  igl::edge_topology(VCoarse, FCoarse, EVCoarse, FECoarse, EFCoarse);
  igl::barycenter(VCoarse, FCoarse, barycenters);
  
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
  
  directional::curl_matching(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, rawFieldCoarse, matchingCoarse, effortCoarse, curlCoarse);
  directional::effort_to_indices(VCoarse, FCoarse, EVCoarse, EFCoarse, effortCoarse, matchingCoarse, N, singVerticesCoarse, singIndicesCoarse);
  directional::combing(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, rawFieldCoarse, matchingCoarse, combedFieldCoarse);
  directional::curl_matching(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, combedFieldCoarse, combedMatchingCoarse, combedEffortCoarse, curlCoarse);
  double curlMax = curlCoarse.maxCoeff();
  std::cout << "Coarse optimized absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing coarse parameterization"<<endl;
  parameterize_mesh(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, singVerticesCoarse, rawFieldCoarse, matchingCoarse, VCutCoarse, FCutCoarse, cutFullUVCoarse);
  
  cout<<"Doing field subdivision"<<endl;
  //Subdividing field and mesh. Curl precision should result in very low curl on the fine mesh as well
  directional::subdivide_field(VCoarse, FCoarse, rawFieldCoarse, targetLevel,  VFine, FFine, rawFieldFine);
  
  //directional::subdivide_field(VCoarse, FCoarse, EVCoarse, EFCoarse, combedFieldCoarse, combedMatchingCoarse, targetLevel, VFine, FFine, EVFine, EFFine, rawFieldFine, matchingFine);

  igl::edge_topology(VFine, FFine, EVFine, FEFine, EFFine);
    
  directional::curl_matching(VFine, FFine, EVFine, EFFine, FEFine, rawFieldFine, matchingFine, effortFine, curlFine);
  directional::effort_to_indices(VFine, FFine, EVFine, EFFine, effortFine, matchingFine, N, singVerticesFine, singIndicesFine);
  directional::combing(VFine, FFine, EVFine, EFFine, FEFine, rawFieldFine, matchingFine, combedFieldFine);
  directional::curl_matching(VFine, FFine, EVFine, EFFine, FEFine, combedFieldFine, combedMatchingFine, combedEffortFine, curlFine);
  curlMax = curlFine.maxCoeff();
  std::cout << "Fine subdivided absolute curl: " << curlMax << std::endl;
  
  cout<<"Doing fine parameterization"<<endl;
  parameterize_mesh(VFine, FFine, EVFine, EFFine, FEFine, singVerticesFine, rawFieldFine, matchingFine, VCutFine, FCutFine, cutFullUVFine);
  
  cout<<"Done!"<<endl;
  
  
  //field mesh
  viewer.append_mesh();
  directional::glyph_lines_raw(VCoarse, FCoarse, combedFieldCoarse, directional::indexed_glyph_colors(combedFieldCoarse), VFieldCoarse, FFieldCoarse, CFieldCoarse,2.0);
  directional::glyph_lines_raw(VFine, FFine, combedFieldFine, directional::indexed_glyph_colors(combedFieldFine), VFieldFine, FFieldFine, CFieldFine,1.5);
  
  //singularity mesh
  viewer.append_mesh();
  directional::singularity_spheres(VCoarse, FCoarse, N, singVerticesCoarse, singIndicesCoarse, VSingsCoarse, FSingsCoarse, CSingsCoarse,1.0);
  directional::singularity_spheres(VFine, FFine, N, singVerticesFine, singIndicesFine, VSingsFine, FSingsFine, CSingsFine,2.0);
  
  //seams mesh
  viewer.append_mesh();
  directional::seam_lines(VCoarse,FCoarse,EVCoarse,combedMatchingCoarse, VSeamsCoarse,FSeamsCoarse,CSeamsCoarse,1.0);
  directional::seam_lines(VFine,FFine,EVFine,combedMatchingFine, VSeamsFine,FSeamsFine,CSeamsFine,2.0);
  
  // Update view
  update_view();
  
  viewer.selected_data_index = 0;
  
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
