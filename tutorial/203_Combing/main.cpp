#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/principal_combing.h>
#include <directional/line_cylinders.h>


int currF=0, N;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXd VMesh, VField, VSings, CSings, CField;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::RowVector3d rawGlyphColor;
igl::opengl::glfw::Viewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd glyphPrincipalColors(5,3);


void update_raw_field_mesh()
{
  
  Eigen::MatrixXd fullGlyphColor(FMesh.rows(),3*N);
  for (int i=0;i<FMesh.rows();i++)
    for (int j=0;j<N;j++)
      fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);

  directional::glyph_lines_raw(VMesh, FMesh, (viewer.data_list[3].show_faces ? combedField : rawField), fullGlyphColor,VField, FField, CField);
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_lines = false;
  
  
}
  


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewer.data_list[3].show_faces = !viewer.data_list[3].show_faces; update_raw_field_mesh(); break;
    case '2': viewer.data_list[2].show_faces = !viewer.data_list[2].show_faces; break;
  }
  return true;
}


int main()
{
  std::cout <<
  "  1        Toggle raw field/Combed field" << std::endl <<
  "  2        Show/hide singularities" << std::endl;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", N, rawField);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  //computing
  directional::principal_combing(VMesh,FMesh, EV, EF, FE, rawField, combedField, combedMatching, combedEffort);
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N,prinIndices);
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<VMesh.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVertices.resize(singVerticesList.size());
  singIndices.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVertices(i)=singVerticesList[i];
    singIndices(i)=singIndicesList[i];
  }
  
  igl::barycenter(VMesh, FMesh, barycenters);
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0;
  
  //triangle mesh setup
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(Eigen::RowVector3d::Constant(3,1.0));

  //apending and updating raw field mesh
  viewer.append_mesh();
  update_raw_field_mesh();
  
  //singularity mesh
  viewer.append_mesh();
  directional::singularity_spheres(VMesh, FMesh, singVertices, singIndices, directional::defaultSingularityColors(N), VSings, FSings, CSings);
  viewer.data().set_mesh(VSings, FSings);
  viewer.data().set_colors(CSings);
  viewer.data_list[2].show_faces = true;
  viewer.data_list[2].show_lines = false;
  
  
  //seam mesh
  double l = igl::avg_edge_length(VMesh, FMesh);
  std::vector<int> seamEdges;
  for (int i=0;i<EV.rows();i++)
    if (combedMatching(i)!=0)
      seamEdges.push_back(i);
    
  Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
  for (int i=0;i<seamEdges.size();i++){
    P1.row(i)=VMesh.row(EV(seamEdges[i],0));
    P2.row(i)=VMesh.row(EV(seamEdges[i],1));
  }
  
  Eigen::MatrixXd VSeam, CSeam;
  Eigen::MatrixXi FSeam;
  directional::line_cylinders(P1, P2, l/25.0, Eigen::MatrixXd::Constant(FMesh.rows(), 3, 0.0), 6, VSeam, FSeam, CSeam);
  
  viewer.append_mesh();
  viewer.data().set_mesh(VSeam, FSeam);
  viewer.data().set_colors(CSeam);
  viewer.data_list[3].show_faces = false;
  viewer.data_list[3].show_lines = false;
  
  viewer.selected_data_index=0;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


