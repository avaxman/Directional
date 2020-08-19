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
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/subdivide_field.h>

#include "tutorial_shared_path.h"
#include "directional/cut_mesh_with_singularities.h"
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


// The igl viewer
igl::opengl::glfw::Viewer viewer;

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;



//for averaging curl to faces, for visualization
Eigen::SparseMatrix<double> AE2FCoarse, AE2FFine;

int N;
int targetLevel = 1;
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
    igl::jet(faceCurl, 0.0,faceCurl.maxCoeff(), CMeshCoarse);
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
    igl::jet(faceCurl, 0.0,faceCurl.maxCoeff(), CMeshFine);
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

void update_raw_field_mesh()
{
  
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
  
  Eigen::MatrixXi symmFunc(4, 2);
  symmFunc << 1, 0,
  0, 1,
  -1, 0,
  0, -1;
  
  directional::setup_parameterization(symmFunc, VWhole, FWhole, EV, EF, FE, combedMatching, singVertices, pd, VCut, FCut);
  
  double lengthRatio = 0.01;
  bool isInteger = false;  //do not do translational seamless.
  std::cout << "Solving parameterization" << std::endl;
  Eigen::MatrixXd cutReducedUV,  cornerWholeUV;
  directional::parameterize(VWhole, FWhole, FE, combedField, lengthRatio, pd, VCut, FCut, isInteger, cutReducedUV, cutFullUV, cornerWholeUV);
  cutFullUV = cutFullUV.block(0, 0, cutFullUV.rows(), 2);
  std::cout << "Done!" << std::endl;
  
  /*output.clear();
   output.set_mesh(VMeshCut, FMeshCut);
   output.set_uv(cutFullUV);
   output.show_texture = true;
   output.show_lines = false;*/
}

// Compute curl statistics on directional
/*void computeDirectionalCurl(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& EV, const Eigen::MatrixXi& EF, const Eigen::MatrixXd& rawField, const Eigen::VectorXi& matching,
 Eigen::VectorXd& curlColumn, double& maxSqNorm, Eigen::VectorXd& l2CurlOnFaces)
 {
 Eigen::MatrixXi EI_local, SFE_local;
 directional::shm_edge_topology(F, EV, EF, EI_local, SFE_local);
 Eigen::SparseMatrix<double> C, columndirectional_to_g2, avgToFaces;
 directional::Matched_Curl(EF, SFE_local, EI_local, matching, F.rows(), N, C);
 directional::columndirectional_to_gamma2_matrix(V, F, EV, SFE_local, EF, N, columndirectional_to_g2);
 Eigen::VectorXd columnDirectionalFine;
 directional::rawfield_to_columndirectional(rawField, N, columnDirectionalFine);
 // Output
 curlColumn = C * columndirectional_to_g2 * columnDirectionalFine;
 
 // Compute l2 norm
 Eigen::MatrixXd curlPerEdge(curlColumn.rows() / N, N);
 for(int n = 0; n < N; ++n)
 {
 curlPerEdge.col(n) = curlColumn.segment(n*curlPerEdge.rows(), curlPerEdge.rows());
 }
 auto sqNorm = curlPerEdge.rowwise().norm();
 maxSqNorm = sqNorm.maxCoeff();
 // Copy to faces
 directional::Edge_To_Face_Average(SFE_local, EF.rows(), avgToFaces);
 l2CurlOnFaces = avgToFaces * sqNorm;
 }*/

// Compute subdivision of coarse optimized field and apply parameterization
/*void construct_fine_level_parameterization()
 {
 {
 Eigen::VectorXd curlColumn, l2CurlOnFaces;
 double maxSqNorm;
 computeDirectionalCurl(VCoarse, FCoarse, EVCoarse, EFCoarse, combedFieldCF, combedMatchingCF, curlColumn, maxSqNorm, l2CurlOnFaces);
 
 std::cout << "Coarse max l2 curl before: " << maxSqNorm << std::endl;
 }
 directional::subdivide_field(VCoarse, FCoarse, EVCoarse, EFCoarse, combedFieldCF, combedMatchingCF, targetLevel, VMesh_fine, FMesh_fine,
 EV_fine, EF_fine, rawfield_fine, matchingCF_fine);
 
 {
 Eigen::VectorXd curlColumn;
 double maxSqNorm;
 computeDirectionalCurl(VMesh_fine, FMesh_fine, EV_fine, EF_fine, rawfield_fine, matchingCF_fine, curlColumn, maxSqNorm, l2Curl_fine);
 
 std::cout << "Fine max l2 curl: " << maxSqNorm << std::endl;
 }
 
 // Set the fine mesh
 resetMeshData(viewer.data_list[DisplayMeshes::FineMesh], VMesh_fine, FMesh_fine, directional::default_mesh_color(), false, false);
 
 // Reconstruct FE where FE(f,0) is now the edge that has vertices F(f,0) and F(f,1) (igl::edge_topology style...)
 Eigen::MatrixXi FE_fine(FMesh_fine.rows(), 3);
 for (int e = 0; e < EV_fine.rows(); ++e)
 {
 for (int s = 0; s < 2; ++s)
 {
 int targetF = EF_fine(e, s);
 for (int j = 0; j < 3; ++j)
 {
 const int otherV = FMesh_fine(targetF, (j + 2) % 3);
 if ( EV_fine(e, 0) != otherV && EV_fine(e, 1) != otherV)
 {
 FE_fine(targetF, j) = e;
 break;
 }
 }
 }
 }
 // Parameterize in the fine level and output to libigl viewer mesh
 // We can reuse singVerticesCF since the coarse level vertex IDs remained unchanged and we assume the singularities to be stationary under the subdivision
 parameterize_mesh(VMesh_fine, FMesh_fine, EV_fine, EF_fine, FE_fine, singVerticesCF, rawfield_fine, matchingCF_fine, viewer.data_list[DisplayMeshes::FineParam]);
 }*/

/*
 * Update meshes for view and make the appropriate meshes visible
 */
/*void update_view()
 {
 switch(viewingMode)
 {
 case ORIGINAL_FIELD:
 {
 directional::glyph_lines_raw(VCoarse, FCoarse, (viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF),
 directional::indexed_glyph_colors((viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF)), VField, FField, CField, 2.0);
 // Reset the mesh data for the field
 resetMeshData(viewer.data_list[1], VField, FField, CField, true, false);
 
 //singularity mesh
 directional::singularity_spheres(VCoarse, FCoarse, N, (viewingMode == ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode == ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), VSings, FSings, CSings, 1.5);
 
 resetMeshData(viewer.data_list[DisplayMeshes::Singularities], VSings, FSings, CSings, true, false);
 
 //seam mesh
 directional::seam_lines(VCoarse, FCoarse, EVCoarse, (viewingMode == ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF), VSeams, FSeams, CSeams, 2.0);
 
 resetMeshData(viewer.data_list[DisplayMeshes::Seams], VSeams, FSeams, CSeams, true, false);
 viewer.data_list[Mesh].set_colors(directional::default_mesh_color());
 break;
 }
 case ORIGINAL_CURL:
 {
 Eigen::VectorXd currCurl = AE2F * curlOrig;
 igl::jet(currCurl, 0.0, curlMaxOrig, CMesh);
 viewer.data_list[DisplayMeshes::Mesh].set_colors(CMesh);
 }
 break;
 case OPTIMIZED_FIELD:
 {
 directional::glyph_lines_raw(VCoarse, FCoarse, (viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF),
 directional::indexed_glyph_colors((viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF)), VField, FField, CField, 2.0);
 // Reset the mesh data for the field
 resetMeshData(viewer.data_list[1], VField, FField, CField, true, false);
 
 //singularity mesh
 directional::singularity_spheres(VCoarse, FCoarse, N, (viewingMode == ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode == ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), VSings, FSings, CSings, 1.5);
 
 resetMeshData(viewer.data_list[DisplayMeshes::Singularities], VSings, FSings, CSings, true, false);
 
 //seam mesh
 directional::seam_lines(VCoarse, FCoarse, EVCoarse, (viewingMode == ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF), VSeams, FSeams, CSeams, 2.0);
 
 resetMeshData(viewer.data_list[DisplayMeshes::Seams], VSeams, FSeams, CSeams, true, false);
 viewer.data_list[Mesh].set_colors(directional::default_mesh_color());
 break;
 }
 case OPTIMIZED_CURL:
 {
 Eigen::VectorXd currCurl = AE2F * curlCF;
 igl::jet(currCurl, 0.0, curlMaxOrig, CMesh);
 viewer.data_list[DisplayMeshes::Mesh].set_colors(CMesh);
 }
 break;
 case FINE_FIELD:
 {
 directional::glyph_lines_raw(VMesh_fine, FMesh_fine, rawfield_fine,
 directional::indexed_glyph_colors(rawfield_fine), VField, FField, CField, 2.0);
 // Reset the mesh data for the field
 resetMeshData(viewer.data_list[FineField], VField, FField, CField, true, false);
 viewer.data_list[FineMesh].set_colors(directional::default_mesh_color());
 }
 break;
 case FINE_CURL:
 {
 Eigen::MatrixXd Cs;
 igl::jet(l2Curl_fine, 0.0, curlMaxOrig * 0.25, Cs);
 viewer.data_list[DisplayMeshes::FineMesh].set_colors(Cs);
 break;
 }
 /*case COARSE_PARAM:
 {
 parameterize_cf(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, singVerticesCF, combedFieldCF, combedMatchingCF, viewer.data_list[CoarseParam]);
 break;
 }
 case FINE_PARAM:
 viewer.data_list[FineParam].show_faces = true;
 break;
 /*
 COARSE_PARAM,
 FINE_PARAM
 default:
 break;
 }
 for(DisplayMeshes m = Mesh; m < DisplayMeshCount; m = static_cast<DisplayMeshes>(m+1))
 {
 bool show = meshesPerView[viewingMode].find(m) != meshesPerView[viewingMode].end();
 viewer.data_list[m].show_faces = show;
 }
 
 }*/

/*bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
 {
 // Set the view mode
 if ((key >= '1') && (key <= '6'))
 {
 viewingMode = static_cast<ViewingModes>(key - '1');
 }
 
 
 
 //do a batch of iterations
 /*if (key == 'A')
 {
 printf("--Improving Curl--\n");
 // Batch of 5 iterations
 for (int bi = 0; bi < 60; ++bi)
 {
 
 printf("\n\n **** Batch %d ****\n", iter);
 // Solve
 directional::polycurl_reduction_solve(pcrdata, params, rawFieldCF, iter == 0);
 ++iter;
 // Adjust the smoothness weight
 params.wSmooth *= params.redFactor_wsmooth;
 }
 
 // Reconstruct matching that minimizes the curl
 directional::curl_matching(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, rawFieldCF, matchingCF, effortCF, curlCF);
 // Get indices from the effort for the matching
 directional::effort_to_indices(VCoarse, FCoarse, EVCoarse, EFCoarse, effortCF, matchingCF, N, singVerticesCF, singIndicesCF);
 // Comb the field
 directional::combing(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, rawFieldCF, matchingCF, combedFieldCF);
 // Apply curl matching on the combed fields
 directional::curl_matching(VCoarse, FCoarse, EVCoarse, EFCoarse, FECoarse, combedFieldCF, combedMatchingCF, combedEffortCF, curlCF);
 curlMax = curlCF.maxCoeff();
 std::cout << "curlMax optimized: " << curlMax << std::endl;
 }
 if(key == 'S')
 {
 std::cout << "Subdividing and parameterizing " << std::endl;
 construct_fine_level_parameterization();
 }
 
 // Write the current coarse field to file.
 if (key == 'W') {
 if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision-cf.rawfield", rawFieldCF))
 std::cout << "Saved raw field" << std::endl;
 else
 std::cout << "Unable to save raw field. " << std::endl;
 }
 
 //update_triangle_mesh();
 update_view();
 return false;
 }*/

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
  
  // Load a mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.off", VCoarse, FCoarse);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.rawfield", N, rawFieldCoarse);
  
  igl::edge_topology(VCoarse, FCoarse, EVCoarse, FECoarse, EFCoarse);
  
  // compute barycenters
  igl::barycenter(VCoarse, FCoarse, barycenters);
  
  //Reducing curl from coarse mesh
  Eigen::VectorXi b; b.resize(1); b << 0;
  Eigen::MatrixXd bc; bc.resize(1, 6); bc << rawFieldCoarse.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b << 1;
  // Precompute polycurl reduction data.
  directional::polycurl_reduction_precompute(VCoarse, FCoarse, b, bc, blevel, rawFieldCoarse, pcrdata);
  
  printf("--Improving Curl--\n");
  //Performing 60 Batches of 5 iterations of curl reduction
  for (int bi = 0; bi < 5; ++bi)
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
  directional::subdivide_field(VCoarse, FCoarse, rawFieldCoarse, 1,  VFine, FFine, rawFieldFine);
  
  //directional::subdivide_field(VCoarse, FCoarse, EVCoarse, EFCoarse, combedFieldCoarse, combedMatchingCoarse, targetLevel, VFine, FFine, EVFine, EFFine, rawFieldFine, matchingFine);

  
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
  directional::glyph_lines_raw(VCoarse, FCoarse, combedFieldCoarse, directional::indexed_glyph_colors(combedFieldCoarse), VFieldCoarse, FFieldCoarse, CFieldCoarse,2.5);
  directional::glyph_lines_raw(VFine, FFine, combedFieldFine, directional::indexed_glyph_colors(combedFieldFine), VFieldFine, FFieldFine, CFieldFine,2.5);
  
  //singularity mesh
  viewer.append_mesh();
  directional::singularity_spheres(VCoarse, FCoarse, N, singVerticesCoarse, singIndicesCoarse, VSingsCoarse, FSingsCoarse, CSingsCoarse,2.5);
  directional::singularity_spheres(VFine, FFine, N, singVerticesFine, singIndicesFine, VSingsFine, FSingsFine, CSingsFine,2.5);
  
  //seams mesh
  viewer.append_mesh();
  directional::seam_lines(VCoarse,FCoarse,EVCoarse,combedMatchingCoarse, VSeamsCoarse,FSeamsCoarse,CSeamsCoarse,2.5);
  directional::seam_lines(VFine,FFine,EVFine,combedMatchingFine, VSeamsFine,FSeamsFine,CSeamsFine,2.5);
  
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
