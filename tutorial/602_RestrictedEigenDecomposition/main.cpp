#include <iostream>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/edge_topology.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/singularity_spheres.h>
#include <directional/combing.h>
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/power_field.h>
#include <directional/cut_mesh_with_singularities.h>
#include <igl/grad.h>
#include "directional/get_directional_subdivision_suite.h"
// Import matlab interface for eigen computation
#ifdef LIBIGL_WITH_MATLAB
#include <igl/matlab/matlabinterface.h>
#endif

using namespace std;

int subdivisionLevel = 3;

int numberOfEigenValues = 15;

Eigen::VectorXi matchingCoarse, matchingFine, combedMatchingCoarse, combedMatchingFine, singVerticesCoarse, singVerticesFine, singIndicesCoarse, singIndicesFine;
Eigen::MatrixXi FCoarse, FFine, FCutCoarse, FCutFine, FFieldCoarse, FSingsCoarse, FSeamsCoarse, FFieldFine, FSingsFine, FSeamsFine;
Eigen::MatrixXi EVCoarse, EFCoarse, FECoarse,EVFine, EFFine, FEFine;

// Separate matrices for different visualizations
Eigen::MatrixXd VCoarse, VCutCoarse, VFieldCoarse, VSingsCoarse, VSeamsCoarse, VFine, VCutFine, VFieldFine, VSingsFine, VSeamsFine, barycenters;
Eigen::MatrixXd CMeshCoarse, CFieldCoarse, CSingsCoarse, CSeamsCoarse, CMeshFine, CFieldFine, CSingsFine, CSeamsFine;
Eigen::MatrixXd rawFieldCoarse, rawFieldFine;
Eigen::MatrixXd combedFieldCoarse, combedFieldFine;
Eigen::VectorXd curlCoarse, curlFine; // norm of curl per edge

// The igl viewer
igl::opengl::glfw::Viewer viewer;

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
  //directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.rawfield", N, rawFieldCoarse);
 
  igl::edge_topology(VCoarse, FCoarse, EVCoarse, FECoarse, EFCoarse);
  igl::barycenter(VCoarse, FCoarse, barycenters);
  
  Eigen::SparseMatrix<double> PCoarse, S_Gamma, PInvFine, S_epsstar, S_0, S_1, S_2,WCoarse, WInvFine;
  // Get subdivision matrix
  directional::get_pcvf_subdivision_suite(VCoarse, FCoarse, EVCoarse, subdivisionLevel, S_epsstar, S_0, S_1, S_2, S_Gamma, WCoarse, PCoarse, EVFine, FFine, WInvFine, PInvFine);

  // Compute the SHM Laplacian and Co-Laplacian
  {
      // Alternative construction of same Laplacian
      Eigen::SparseMatrix<double> G;
      // Gradient/Divergence
      igl::grad(VCoarse, FCoarse, G);
      // Diagonal per-triangle "mass matrix"
      Eigen::VectorXd dblAFine;
      igl::doublearea(VFine, FFine, dblAFine);
      // Place areas along diagonal #dim times. (const auto& ?, does that work?)
      const auto & T = 1.*(dblAFine.replicate(3, 1)*0.5).asDiagonal();

      // SHM the mass matrix
      Eigen::SparseMatrix<double> mass = (PInvFine * S_Gamma).transpose() * T * (PInvFine * S_Gamma);

      // Laplacian K built as discrete divergence of gradient or equivalently
      // discrete Dirichelet energy Hessian
      Eigen::SparseMatrix<double> L = -G.transpose() * mass * G;
      // Compute mass matrix for generalized solutions
      Eigen::SparseMatrix<double> M_v;

      igl::massmatrix(VFine, FFine, igl::MASSMATRIX_TYPE_BARYCENTRIC, M_v);


      // We solve the generalized eigen value problem  \lambda M \phi = L \phi for eigenvector \phi and eigenvalue \lambda, where M is a mass matrix.
#ifdef LIBIGL_WITH_MATLAB
      // Matlab instance
      Engine* engine;
      // Launch MATLAB
      igl::matlab::mlinit(&engine);
      // Send Laplacian matrix to matlab
      igl::matlab::mlsetmatrix(&engine, "L", L);
      // Send mass matrix
      igl::matlab::mlsetmatrix(&engine, "L", L);

      std::stringstream code;
      code << "[eigenValues, eigenVectors] = eigs(-L,M," << numberOfEigenValues << ", 1e-5)";

      // Extract the first 10 eigenvectors
      igl::matlab::mleval(&engine, code.str().c_str());

      // Retrieve the result
      igl::matlab::mlgetmatrix(&engine, "EV", EV);
#else 
      // We actually need the largest eigenvalues, so IGL support is somewhat limited...
#endif
  }
  
  
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
