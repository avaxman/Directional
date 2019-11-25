#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <igl/png/readPNG.h>
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
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/cut_mesh_with_singularities.h>


int N;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams;
Eigen::MatrixXd CField, CSeams, CSings;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd cutUVFull, cutUVRot;
igl::opengl::glfw::Viewer viewer;

typedef enum {FIELD, ROT_PARAMETERIZATION, FULL_PARAMETERIZATION, UV_COORDS} ViewingModes;
ViewingModes viewingMode=FIELD;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;


// Create a texture that hides the integer translation in the parametrization
void hexline_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                     Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                     Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
{
   Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A;
   igl::png::readPNG(TUTORIAL_SHARED_PATH "/hexpattern2.png", texture_R, texture_G, texture_B, A);
}


void update_triangle_mesh()
{
  if (viewingMode==FIELD){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_face_based(false);
    viewer.data_list[0].show_texture=false;
    viewer.data_list[0].show_lines=false;
  } else if ((viewingMode==ROT_PARAMETERIZATION) || (viewingMode==FULL_PARAMETERIZATION)){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshCut, FMeshCut);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_uv(viewingMode==ROT_PARAMETERIZATION ? cutUVRot : cutUVFull);
    viewer.data_list[0].set_texture(texture_R, texture_G, texture_B);
    viewer.data_list[0].set_face_based(true);
    viewer.data_list[0].show_texture=true;
    viewer.data_list[0].show_lines=false;
  } else {
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(cutUVFull, FMeshCut);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_face_based(false);
    viewer.data_list[0].show_texture=false;
    viewer.data_list[0].show_lines=true;
  }
}

void update_raw_field_mesh()
{
  for (int i=1;i<=4;i++)  //hide all other meshes
    viewer.data_list[i].show_faces=(viewingMode==FIELD);
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = ROT_PARAMETERIZATION; break;
    case '3': viewingMode = FULL_PARAMETERIZATION; break;
    case '4': viewingMode = UV_COORDS; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-rot-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVRot, FMeshCut);
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-full-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVFull, FMeshCut);
      break;
  }
  update_triangle_mesh();
  update_raw_field_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show textured rotationally-seamless parameterization mesh" << std::endl <<
  "  3  Show textured fully-seamless parameterization mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  
  hexline_texture(texture_R, texture_G, texture_B);
  igl::readOFF(TUTORIAL_SHARED_PATH "/aqua-center.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/aqua-center.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  //combing and cutting
  Eigen::VectorXd curlNorm;
  directional::curl_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort, curlNorm);
  
  directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,singVertices, singIndices);
  
  directional::ParameterizationData pd;
  directional::cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
  directional::combing(VMeshWhole,FMeshWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField, combedMatching);
  //directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, combedField, combedMatching, combedEffort);
  std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;
  
  std::cout<<"Setting up parameterization"<<std::endl;
  
  directional::setup_parameterization(N, VMeshWhole, FMeshWhole, EV, EF, FE, combedMatching, singVertices, pd, VMeshCut, FMeshCut);
  
  double lengthRatio=0.007;
  bool isInteger = false;  //do not do translational seamless.
  std::cout<<"Solving rotationally-seamless parameterization"<<std::endl;
  directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, cutUVRot);
  std::cout<<"Done!"<<std::endl;
  
  isInteger = true;  //do not do translational seamless.
  std::cout<<"Solving fully-seamless parameterization"<<std::endl;
  directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger,  cutUVFull);
  
  // convert UVs for opengl
  Eigen::Matrix2d c;
  c << 1., -1. / 2., 0. , -3. / 2.;

  for(int i = 0; i < cutUVFull.rows(); i++)
  {
    Eigen::Vector2d t = cutUVFull.row(i);
    cutUVFull.row(i) = c * t;
  }

  for(int i = 0; i < cutUVRot.rows(); i++)
  {
    Eigen::Vector2d t = cutUVRot.row(i);
    cutUVRot.row(i) = c * t;
  }

  std::cout<<"Done!"<<std::endl;
  
  //raw field mesh
  viewer.append_mesh();
  directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, directional::indexed_glyph_colors(combedField), VField, FField, CField,1.0);
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
  
  //singularity mesh
  viewer.append_mesh();
  directional::singularity_spheres(VMeshWhole, FMeshWhole, N, singVertices, singIndices, VSings, FSings, CSings,2.5);
  
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].set_colors(CSings);
  viewer.data_list[2].show_faces = true;
  viewer.data_list[2].show_lines = false;
  
  //seams mesh
  viewer.append_mesh();
  
  Eigen::VectorXi isSeam=Eigen::VectorXi::Zero(EV.rows());
  for (int i=0;i<FE.rows();i++)
    for (int j=0;j<3;j++)
      if (pd.face2cut(i,j))
        isSeam(FE(i,j))=1;
  directional::seam_lines(VMeshWhole, FMeshWhole, EV, combedMatching, VSeams, FSeams, CSeams,2.5);
  
  viewer.data_list[3].clear();
  viewer.data_list[3].set_mesh(VSeams, FSeams);
  viewer.data_list[3].set_colors(CSeams);
  viewer.data_list[3].show_faces = true;
  viewer.data_list[3].show_lines = false;
  
  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.data_list[0].show_lines=false;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


