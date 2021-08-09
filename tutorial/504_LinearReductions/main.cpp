#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/line_cylinders.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/branched_isolines.h>


int N;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams;
Eigen::MatrixXd CField, CSeams, CSings;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd NFunctionSign, NFunctionTri, NCornerFunc;
igl::opengl::glfw::Viewer viewer;

typedef enum {FIELD, SIGN_SYMMETRY, TRI_SYMMETRY} ViewingModes;
ViewingModes viewingMode=FIELD;


void update_triangle_mesh()
{
  if (viewingMode==FIELD){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_face_based(false);
    viewer.data_list[0].show_lines=false;
  } else if ((viewingMode==SIGN_SYMMETRY) || (viewingMode==TRI_SYMMETRY)){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshCut, FMeshCut);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_face_based(false);
    viewer.data_list[0].show_lines=false;
  }
  
  viewer.data_list[4].show_faces=(viewingMode==SIGN_SYMMETRY);
  viewer.data_list[5].show_faces=(viewingMode==TRI_SYMMETRY);
}

void update_raw_field_mesh()
{
  for (int i=1;i<=3;i++)  //hide all other meshes
    viewer.data_list[i].show_faces=(viewingMode==FIELD);
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = SIGN_SYMMETRY; break;
    case '3': viewingMode = TRI_SYMMETRY; break;
  }
  update_triangle_mesh();
  update_raw_field_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show only sign-symmetric integrated functions" << std::endl <<
  "  3  Show triangular-symmetric integrated functions" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/dome.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/dome-6.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  //combing and cutting
  directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,singVertices, singIndices);
  
  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  
  intData.verbose=false;
  intData.integralSeamless=true;

  std::cout<<"Free (sign-symmetric) Integrating..."<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut, NFunctionSign, NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  
  std::cout<<"Solving triangular-constrained integration..."<<std::endl;
  intData.set_triangular_symmetry(N);
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut,  NFunctionTri, NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  //raw field mesh
  directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, directional::indexed_glyph_colors(combedField), VField, FField, CField,1.0);
  viewer.append_mesh();
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
  
  //singularity mesh
  directional::singularity_spheres(VMeshWhole, FMeshWhole, N, singVertices, singIndices, VSings, FSings, CSings,2.5);
  viewer.append_mesh();
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].set_colors(CSings);
  viewer.data_list[2].show_faces = true;
  viewer.data_list[2].show_lines = false;
  
  //seams mesh
  Eigen::VectorXi isSeam=Eigen::VectorXi::Zero(EV.rows());
  for (int i=0;i<FE.rows();i++)
    for (int j=0;j<3;j++)
      if (intData.face2cut(i,j))
        isSeam(FE(i,j))=1;
  directional::seam_lines(VMeshWhole, FMeshWhole, EV, combedMatching, VSeams, FSeams, CSeams,2.5);
  
  viewer.append_mesh();
  viewer.data_list[3].clear();
  viewer.data_list[3].set_mesh(VSeams, FSeams);
  viewer.data_list[3].set_colors(CSeams);
  viewer.data_list[3].show_faces = true;
  viewer.data_list[3].show_lines = false;
  
  //sign-symmetric isolines mesh

  viewer.append_mesh();
  Eigen::MatrixXd VIsoLines, CIsoLines;
  Eigen::MatrixXi FIsoLines;
  directional::branched_isolines(VMeshCut, FMeshCut, NFunctionSign, VIsoLines, FIsoLines, CIsoLines);
  viewer.data_list[4].clear();
  viewer.data_list[4].set_mesh(VIsoLines, FIsoLines);
  viewer.data_list[4].set_colors(CIsoLines);
  viewer.data_list[4].show_faces = false;
  viewer.data_list[4].show_lines = false;
  
  //tri-symmetric isolines mesh
  viewer.append_mesh();
  directional::branched_isolines(VMeshCut, FMeshCut, NFunctionTri, VIsoLines, FIsoLines, CIsoLines);
  viewer.data_list[5].clear();
  viewer.data_list[5].set_mesh(VIsoLines, FIsoLines);
  viewer.data_list[5].set_colors(CIsoLines);
  viewer.data_list[5].show_faces = false;
  viewer.data_list[5].show_lines = false;
  
  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.data_list[0].show_lines=false;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


