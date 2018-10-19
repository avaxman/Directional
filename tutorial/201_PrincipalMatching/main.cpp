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


int currF=0, N;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField, barycenters;
Eigen::VectorXd effort;
Eigen::RowVector3d defaultGlyphColor;
igl::opengl::glfw::Viewer viewer;
Eigen::VectorXi matching, indices;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singVertices;

Eigen::MatrixXd glyphPrincipalColors(5,3);

bool zeroPressed=false;



void update_triangle_mesh()
{
  
  CMesh=Eigen::MatrixXd::Constant(FMesh.rows(), 3, 1.0);
  CMesh.row(currF)<<0.5,0.1,0.1;
  
  viewer.data_list[0].set_colors(CMesh);
}

void update_raw_field_mesh()
{
  Eigen::Vector3i otherFaces;
  Eigen::Vector3i zeroInFace;
  for (int i=0;i<3;i++){
    otherFaces(i)=(EF(FE(currF,i),0)==currF ? EF(FE(currF,i),1) : EF(FE(currF,i),0));
    zeroInFace(i)=(EF(FE(currF,i),0)==currF ? matching(FE(currF,i)) : -matching(FE(currF,i)));
  }
  
  Eigen::MatrixXd rawGlyphColor(FMesh.rows(),3*N);
  for (int i=0;i<FMesh.rows();i++)
    for (int j=0;j<N;j++)
      rawGlyphColor.block(i,3*j,1,3)<<defaultGlyphColor;
  
  //for marked face and its neighbors, coloring adjacent matchings
  for (int i=0;i<N;i++){
    rawGlyphColor.block(currF,3*i,1,3)<<glyphPrincipalColors.row(i);
    for (int j=0;j<3;j++)
      if (otherFaces(j)!=-1)  //boundary edge
        rawGlyphColor.block(otherFaces(j),3*((i+zeroInFace(j)+N)%N),1,3)<<glyphPrincipalColors.row(i);
  }
  
  directional::glyph_lines_raw(VMesh, FMesh, rawField, rawGlyphColor, VField, FField, CField);
  
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_lines=false;
}



bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=false; break;
  }
  return true;
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
    case GLFW_KEY_SPACE:
      viewer.data_list[2].show_faces=!viewer.data_list[2].show_faces;
      break;
  }
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifiers)
{
  if (!zeroPressed)
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
                               viewer.core.proj, viewer.core.viewport, VMesh, FMesh, fid, bc))
  {
    
    //choosing face
    if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
      currF=fid;
      update_triangle_mesh();
      update_raw_field_mesh();
      return true;
    }
  }
  return false;
};

int main()
{
  std::cout <<
  "  0+Left button    Choose face" << std::endl <<
  "  <space>  Show/hide singularities" << std::endl;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", N, rawField);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  //computing
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching,N,prinIndices);
  
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
  
  defaultGlyphColor <<0.0, 0.2, 1.0;
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0;
  
  //triangle mesh setup
  viewer.data().set_mesh(VMesh, FMesh);
  update_triangle_mesh();
  
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
  
  viewer.selected_data_index=0;
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}

