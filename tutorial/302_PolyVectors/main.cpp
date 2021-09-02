#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

int currF, currVec;
Eigen::VectorXi b, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V, barycenters;
Eigen::MatrixXd CMesh;
Eigen::MatrixXd rawField,representative, bc;
Eigen::MatrixXcd pvField;
directional::DirectionalViewer viewer;

int N = 3;

//User input variables
int cur = 0;
bool normalized = false;
bool zeroPressed = false;

void update_triangle_mesh()
{
  Eigen::MatrixXd CMesh=directional::DirectionalViewer::default_mesh_color().replicate(F.rows(),1);
  for (int i = 0; i < b.rows(); i++)
    CMesh.row(b(i)) = directional::DirectionalViewer::selected_face_color();
  
  viewer.set_mesh_colors(CMesh);
}

void recompute_field()
{
  directional::polyvector_field(V, F, b, bc, N, pvField);
}

void update_raw_field_mesh()
{
  directional::polyvector_to_raw(V, F, pvField, N, rawField);
  if (normalized)
    for(int n = 0; n < N; n++)
      rawField.middleCols(n*3, 3).rowwise().normalize();
  
  directional::principal_matching(V, F, EV, EF, FE, rawField, matching, effort, singVertices, singIndices);

  
  /*Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(F.rows(),N);
  if (b.rows()!=0){
    glyphColors.row(b(b.rows()-1))=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,N);
    glyphColors.block(b(b.rows()-1),3*currVec,1,3)=directional::DirectionalViewer::selected_vector_glyph_color();
  }*/
  
  viewer.set_field(rawField);
  viewer.set_singularities(singVertices, singIndices);
  if (b.size()!=0){
    viewer.set_selected_faces(b);
    viewer.set_selected_vector(b(b.rows()-1), currVec);
  }
  
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
  double y = viewer.core().viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                               viewer.core().proj, viewer.core().viewport, V, F, fid, baryInFace))
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
      Eigen::RowVector3d newVec =(V.row(F(fid, 0)) * baryInFace(0) +
                                  V.row(F(fid, 1)) * baryInFace(1) +
                                  V.row(F(fid, 2)) * baryInFace(2) - barycenters.row(fid)).normalized();
      
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
  igl::readOFF(TUTORIAL_SHARED_PATH "/fandisk.off", V, F);
  igl::edge_topology(V, F, EV, FE, EF);
  igl::barycenter(V, F, barycenters);
  
  b.resize(0);
  bc.resize(0, 3*N);
  
  //triangle mesh setup
  viewer.set_mesh(V, F);
  recompute_field();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
