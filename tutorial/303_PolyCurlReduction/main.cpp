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
#include <directional/glyph_lines_raw.h>
#include <directional/line_cylinders.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/curl_combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>

#include "tutorial_shared_path.h"

using namespace std;

Eigen::VectorXi matchingOrig, matchingCF, combedMatchingOrig, combedMatchingCF;
Eigen::VectorXd effortOrig, effortCF, combedEffortOrig, combedEffortCF;
Eigen::MatrixXi FMesh, FField, FSings, FSeams;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings, VSeams, barycenters;
Eigen::MatrixXd CMesh, CField, CSings, CSeams;
Eigen::MatrixXd rawFieldOrig, rawFieldCF;
Eigen::MatrixXd combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd glyphPrincipalColors(5,3);

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;

// "Subdivided" mesh obtained by splitting each triangle into 3 (only needed for display)
//Eigen::MatrixXd Vbs;
//Eigen::MatrixXi Fbs;

//for averagin curl to faces, for visualization
Eigen::SparseMatrix<double> AE2F;

// Scale for visualizing textures
//double uv_scale;
double curlMax;
int N;


// The set of parameters for calculating the curl-free fields
directional::polycurl_reduction_parameters params;

// Solver data (needed for precomputation)
directional::PolyCurlReductionSolverData<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd> pcrdata;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;


typedef enum {ORIGINAL_FIELD, ORIGINAL_CURL, OPTIMIZED_FIELD, OPTIMIZED_CURL} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;

int iter = 0;

// Create a texture that hides the integer translation in the parametrization
void line_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
{
  unsigned size = 128;
  unsigned size2 = size/2;
  unsigned lineWidth = 3;
  texture_R.setConstant(size, size, 255);
  for (unsigned i=0; i<size; ++i)
    for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
      texture_R(i,j) = 0;
  for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
    for (unsigned j=0; j<size; ++j)
      texture_R(i,j) = 0;
  
  texture_G = texture_R;
  texture_B = texture_R;
}


void update_triangle_mesh()
{
  if ((viewingMode ==ORIGINAL_FIELD)||(viewingMode ==OPTIMIZED_FIELD))
    viewer.data_list[0].set_colors(Eigen::RowVector3d::Constant(3,1.0));
  else{  //curl viewing - currently averaged to the face
    Eigen::VectorXd currCurl = AE2F*(viewingMode==ORIGINAL_CURL ? curlOrig: curlCF);
    igl::jet(currCurl, 0.0,curlMax, CMesh);
    viewer.data_list[0].set_colors(CMesh);
  }
}


void update_raw_field_mesh()
{
  using namespace std;
  using namespace Eigen;
  
  if ((viewingMode==ORIGINAL_CURL) || (viewingMode==OPTIMIZED_CURL)){
    for (int i=1;i<=3;i++){  //hide all other meshes
      viewer.data_list[i].show_faces=false;
      viewer.data_list[i].show_lines = false;
    }
  } else {
    Eigen::MatrixXd fullGlyphColor(FMesh.rows(),3*N);
    for (int i=0;i<FMesh.rows();i++)
      for (int j=0;j<N;j++)
        fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
    
    directional::glyph_lines_raw(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF), fullGlyphColor,VField, FField, CField);
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    //singularity mesh
    directional::singularity_spheres(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode==ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), directional::defaultSingularityColors(N), VSings, FSings, CSings);
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings, FSings);
    viewer.data_list[2].set_colors(CSings);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    //seam mesh
    double avgEdgeLength = igl::avg_edge_length(VMesh, FMesh);
    std::vector<int> seamEdges;
    for (int i=0;i<EV.rows();i++)
      if ((viewingMode==ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF)(i)!=0)
        seamEdges.push_back(i);
    
    Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
    for (int i=0;i<seamEdges.size();i++){
      P1.row(i)=VMesh.row(EV(seamEdges[i],0));
      P2.row(i)=VMesh.row(EV(seamEdges[i],1));
    }
    
    directional::line_cylinders(P1, P2, avgEdgeLength/25.0, Eigen::MatrixXd::Constant(FMesh.rows(), 3, 0.0), 6, VSeams, FSeams, CSeams);
    
    viewer.data_list[3].clear();
    viewer.data_list[3].set_mesh(VSeams, FSeams);
    viewer.data_list[3].set_colors(CSeams);
    viewer.data_list[3].show_faces = true;
    viewer.data_list[3].show_lines = false;
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
    directional::curl_matching(VMesh, FMesh,EV, EF, FE, rawFieldCF, matchingCF, effortCF, curlCF);
    directional::effort_to_indices(VMesh,FMesh,EV, EF, effortCF,matchingCF, N,prinIndices);
    
    std::vector<int> singVerticesList;
    std::vector<int> singIndicesList;
    for (int i=0;i<VMesh.rows();i++)
      if (prinIndices(i)!=0){
        singVerticesList.push_back(i);
        singIndicesList.push_back(prinIndices(i));
      }
    
    singVerticesCF.resize(singVerticesList.size());
    singIndicesCF.resize(singIndicesList.size());
    for (int i=0;i<singVerticesList.size();i++){
      singVerticesCF(i)=singVerticesList[i];
      singIndicesCF(i)=singIndicesList[i];
    }
    
    directional::curl_combing(VMesh,FMesh, EV, EF, FE,rawFieldCF, combedFieldCF, combedMatchingCF, combedEffortCF);
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
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka.rawfield", N, rawFieldOrig);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  //global_scale =  .2*igl::avg_edge_length(VMesh, FMesh);
  //uv_scale = 0.6/igl::avg_edge_length(V, F);
  igl::barycenter(VMesh, FMesh, barycenters);
  
  Eigen::VectorXi prinIndices;
  directional::curl_matching(VMesh, FMesh,EV, EF, FE, rawFieldOrig, matchingOrig, effortOrig, curlOrig);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effortOrig,matchingOrig, N,prinIndices);
  
  curlMax= curlOrig.maxCoeff();
  std:: cout<<"curlMax original: "<<curlMax<<std::endl;
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<VMesh.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVerticesOrig.resize(singVerticesList.size());
  singIndicesOrig.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVerticesOrig(i)=singVerticesList[i];
    singIndicesOrig(i)=singIndicesList[i];
  }
  
  directional::curl_combing(VMesh,FMesh, EV, EF, FE,rawFieldOrig, combedFieldOrig, combedMatchingOrig, combedEffortOrig);
  
  //trivial constraints
  Eigen::VectorXi b; b.resize(1); b<<0;
  Eigen::MatrixXd bc; bc.resize(1,6); bc<<rawFieldOrig.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b<<1;
  directional::polycurl_reduction_precompute(VMesh, FMesh, b, bc, blevel, rawFieldOrig , pcrdata);
  
  rawFieldCF = rawFieldOrig;
  matchingCF = matchingOrig;
  effortCF = effortOrig;
  combedFieldCF =combedFieldOrig;
  combedMatchingCF =combedMatchingOrig;
  combedEffortCF =combedEffortOrig;
  curlCF = curlOrig;
  singVerticesCF = singVerticesOrig;
  singIndicesCF = singIndicesOrig;
  
  // Replace the standard texture with an integer shift invariant texture
  line_texture(texture_R, texture_G, texture_B);
  
  //triangle mesh setup
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(Eigen::RowVector3d::Constant(3,1.0));
  
  //apending and updating raw field mesh
  viewer.append_mesh();
  
  //singularity mesh
  viewer.append_mesh();
  
  //seam mesh
  viewer.append_mesh();
  
  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.selected_data_index=0;
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0;
  
  //creating the AE2F operator
  std::vector<Eigen::Triplet<double> > AE2FTriplets;
  for (int i=0;i<EF.rows();i++){
    AE2FTriplets.push_back(Eigen::Triplet<double>(EF(i,0), i,1.0));
    AE2FTriplets.push_back(Eigen::Triplet<double>(EF(i,1), i,1.0));
  }
  AE2F.resize(FMesh.rows(), EF.rows());
  AE2F.setFromTriplets(AE2FTriplets.begin(), AE2FTriplets.end());
  
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.launch();
  
  return 0;
}
