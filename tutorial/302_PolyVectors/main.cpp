#include <iostream>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/glyph_lines_raw.h>
#include <directional/singularity_spheres.h>
#include <directional/write_raw_field.h>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>

int currF, currVec;
Eigen::VectorXi b, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings, barycenters;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField,representative, bc;
Eigen::MatrixXcd pvField;
igl::opengl::glfw::Viewer viewer;

int N = 3;

//User input variables
int cur = 0;
bool normalized = false;
bool zeroPressed = false;

void update_triangle_mesh()
{
  Eigen::MatrixXd CMesh=directional::default_mesh_color().replicate(FMesh.rows(),1);
  for (int i = 0; i < b.rows(); i++)
    CMesh.row(b(i)) = directional::selected_face_color();
  
  viewer.data_list[0].set_colors(CMesh);
}

void recompute_field()
{
  Eigen::VectorXi bcSoft;
  Eigen::MatrixXd wSoft;
  Eigen::MatrixXd bSoft;
  directional::polyvector_field(VMesh, FMesh, b, bc, bcSoft, wSoft, bSoft, N, pvField);
}

void update_raw_field_mesh()
{
  directional::polyvector_to_raw(VMesh, FMesh, pvField, N, rawField);
  if (normalized)
    for(int n = 0; n < N; n++)
      rawField.middleCols(n*3, 3).rowwise().normalize();
  
  directional::principal_matching(VMesh, FMesh, EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N, singVertices, singIndices);
  
  directional::singularity_spheres(VMesh, FMesh, N, singVertices, singIndices, VSings, FSings, CSings);
  Eigen::MatrixXd glyphColors=directional::default_glyph_color().replicate(FMesh.rows(),N);
  if (b.rows()!=0){
    glyphColors.row(b(b.rows()-1))=directional::selected_face_glyph_color().replicate(1,N);
    glyphColors.block(b(b.rows()-1),3*currVec,1,3)=directional::selected_vector_glyph_color();
  }
  
  directional::glyph_lines_raw(VMesh, FMesh, rawField, glyphColors,  VField, FField, CField);
  
  if (viewer.data_list.size()<2){
    
    //apending and updating raw field mesh
    viewer.append_mesh();
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    viewer.append_mesh();
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    viewer.selected_data_index = 0;
  }
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].set_colors(CSings);
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
    case '1':
      currVec = (currVec+1)%N;
      update_raw_field_mesh();
      break;
      // Reset the constraints
    case 'R':
      b.resize(0);
      bc.resize(0, 3*N);
      recompute_field();
      update_raw_field_mesh();
      update_triangle_mesh();
      break;
      
      // Toggle normalization
    case 'N':
      normalized = !normalized;
      update_raw_field_mesh();
      break;
      
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/fandisk.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. " << std::endl;
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
  Eigen::Vector3d baryInFace;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
                               viewer.core.proj, viewer.core.viewport, VMesh, FMesh, fid, baryInFace))
  {
    
    //checking if face already exists
    int currConst;
    for (currConst=0; currConst<b.rows(); currConst++)
      if (b(currConst) == fid)
        break;
    
    //choosing face
    if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
      
     
      if (currConst == b.rows())  //new face
      {
        b.conservativeResize(b.rows() + 1);
        bc.conservativeResize(bc.rows() + 1, 3*N);
        b(currConst) = fid;
        bc.row(currConst)=rawField.row(fid);   //copying existing information
      }
      
      currF=fid;
      update_triangle_mesh();
      update_raw_field_mesh();
      return true;
    }
    
    //moving vector within face
    if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Right)&&(currConst!=b.rows())){
      // Calculate direction from the center of the face to the mouse
      Eigen::RowVector3d newVec =(VMesh.row(FMesh(fid, 0)) * baryInFace(0) +
                                  VMesh.row(FMesh(fid, 1)) * baryInFace(1) +
                                  VMesh.row(FMesh(fid, 2)) * baryInFace(2) - barycenters.row(fid)).normalized();
      
      bc.block(currConst, currVec*3, 1,3)=newVec;
      recompute_field();
      update_raw_field_mesh();
      return true;
      
    }
  }
  return false;
};

int main()
{
  
  std::cout <<
  "  0+L-bttn   Choose face" << std::endl <<
  "  0+R-bttn   Edit vector in current face" << std::endl<<
  "  1          Choose vector in current face." << std::endl <<
  "  R          Reset the constraints" << std::endl <<
  "  N          Toggle field normalization" << std::endl;
  
  // Load mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/fandisk.off", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  igl::barycenter(VMesh, FMesh, barycenters);
  
  b.resize(0);
  bc.resize(0, 3*N);
  
  //triangle mesh setup
  viewer.data_list[0].set_mesh(VMesh, FMesh);
  viewer.data_list[0].set_colors(directional::default_mesh_color());
  
  viewer.selected_data_index = 0;
  recompute_field();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
