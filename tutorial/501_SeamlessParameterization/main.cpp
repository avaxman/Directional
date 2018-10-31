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
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/principal_matching.h>
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
Eigen::RowVector3d rawGlyphColor;
igl::opengl::glfw::Viewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd cutUV;


typedef enum {FIELD, PARAMETERIZATION} ViewingModes;
ViewingModes viewingMode=FIELD;


void update_triangle_mesh()
{
  if (viewingMode==FIELD){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].show_texture=false;
    viewer.data_list[0].show_lines=false;
  } else {
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshCut, FMeshCut);
    viewer.data_list[0].set_uv(cutUV);
    viewer.data_list[0].show_texture=true;
    viewer.data_list[0].show_lines=false;
  }
}

void update_raw_field_mesh()
{
  for (int i=1;i<=3;i++)  //hide all other meshes
    viewer.data_list[i].show_faces= (viewingMode==FIELD);
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = PARAMETERIZATION; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/botijo-param.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUV, FMeshCut);
      break;
    //case '2': viewer.data().show_texture=!viewer.data().show_texture; update_mesh(); break;
  }
  update_triangle_mesh();
  update_raw_field_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show textured mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  
  igl::readOBJ(TUTORIAL_SHARED_PATH "/botijo.obj", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/botijo.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  //combing and cutting
  directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,singVertices, singIndices);
  
  directional::ParameterizationData pd;
  directional::cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
  directional::combing(VMeshWhole,FMeshWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField);
   directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, combedField,  combedMatching, combedEffort);
 
  std::cout<<"Setting up parameterization"<<std::endl;
  
  directional::setup_parameterization(N, VMeshWhole, FMeshWhole,  EV, EF, FE, combedMatching, singVertices, pd, VMeshCut, FMeshCut);
  
  double lengthRatio=0.01;
  bool isInteger = false;  //do not do translational seamless.
  std::cout<<"Solving parameterization"<<std::endl;
  directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, cutUV);
  std::cout<<"Done!"<<std::endl;
  
  //apending and updating raw field mesh
  viewer.append_mesh();
  directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, directional::indexed_glyph_colors(combedField), VField, FField, CField,2.5);
  
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
  directional::seam_lines(VMeshWhole,FMeshWhole,EV,combedMatching, VSeams,FSeams,CSeams,2.5);
  
  viewer.data_list[3].clear();
  viewer.data_list[3].set_mesh(VSeams, FSeams);
  viewer.data_list[3].set_colors(CSeams);
  viewer.data_list[3].show_faces = true;
  viewer.data_list[3].show_lines = false;

  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


