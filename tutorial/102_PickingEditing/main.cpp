#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>

int currF, currVec, N;
Eigen::MatrixXi FMesh, FField;
Eigen::MatrixXd VMesh, VField,  barycenters;
Eigen::MatrixXd CMesh, CField;
Eigen::MatrixXd rawField, rawGlyphColor;
Eigen::RowVector3d defaultGlyphColor;
Eigen::RowVector3d selectedFaceGlyphColor;
Eigen::RowVector3d selectedVectorGlyphColor;
igl::opengl::glfw::Viewer viewer;

//User input variables
bool drag = false;
bool zeroPressed = false;


void update_triangle_mesh()
{
  
  CMesh=Eigen::MatrixXd::Constant(FMesh.rows(), 3, 1.0);
  CMesh.row(currF)<<0.5,0.1,0.1;
  
  viewer.data_list[0].set_colors(CMesh);
}

void update_raw_field_mesh()
{
  Eigen::MatrixXd rawGlyphColor(FMesh.rows(),3*N);
  for (int i=0;i<FMesh.rows();i++){
    for (int j=0;j<N;j++){
      if (i==currF)
        rawGlyphColor.block(i,3*j,1,3)<<(j==currVec ? selectedVectorGlyphColor : selectedFaceGlyphColor);
      else
        rawGlyphColor.block(i,3*j,1,3)<<defaultGlyphColor;
    }
  }
  
  directional::glyph_lines_raw(VMesh, FMesh, rawField, rawGlyphColor, VField, FField, CField);
  
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_lines=false;
}

/*void update_mesh()
{
  
  Eigen::MatrixXd fullC(F.rows(),3);
  fullC.col(0)=Eigen::VectorXd::Constant(F.rows(),1.0);
  fullC.col(1)=Eigen::VectorXd::Constant(F.rows(),1.0);
  fullC.col(2)=Eigen::VectorXd::Constant(F.rows(),1.0);
  
  fullC.row(currF)<<0.5,0.1,0.1;
  
  Eigen::MatrixXd fullV=V;
  Eigen::MatrixXi fullF=F;
  
  Eigen::MatrixXd fullGlyphColor(F.rows(),3*N);
  for (int i=0;i<F.rows();i++){
    for (int j=0;j<N;j++){
      if (i==currF)
        fullGlyphColor.block(i,3*j,1,3)<<(j==currVec ? selectedVectorGlyphColor : selectedFaceGlyphColor);
      else
        fullGlyphColor.block(i,3*j,1,3)<<defaultGlyphColor;
    }
  }
  
  directional::glyph_lines_raw(V, F, rawField, fullGlyphColor, false, true, fullV, fullF, fullC);
  
  viewer.data().clear();
  viewer.data().set_face_based(true);
  viewer.data().set_mesh(fullV, fullF);
  viewer.data().set_colors(fullC);
}*/

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
  int borders;
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
      currVec = key - '1';
      update_raw_field_mesh();
      break;
   
      // Toggle field drawing for easier rotation
    case 'D':
      drag = !drag;
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
    //choosing face
    if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Right)&&(fid==currF)){
      // Calculate direction from the center of the face to the mouse
      Eigen::RowVector3d newVec =(VMesh.row(FMesh(fid, 0)) * bc(0) +
                                  VMesh.row(FMesh(fid, 1)) * bc(1) +
                                  VMesh.row(FMesh(fid, 2)) * bc(2) - barycenters.row(fid)).normalized();
      
      rawField.block(currF, currVec*3, 1,3)=newVec;
      update_raw_field_mesh();
      return true;
      
    }
  }
  return false;
};

int main()
{
  igl::readOBJ(TUTORIAL_SHARED_PATH "/torus.obj", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/torus.rawfield", N, rawField);
  
  std::cout <<
  "  1-"<< N <<"  Choose vector in current face." << std::endl <<
  "  0+Left button    Choose face" << std::endl <<
  "  0+Right button   Edit vector in current face" << std::endl;
  
  
  igl::barycenter(VMesh, FMesh, barycenters);
  
  defaultGlyphColor<<0.0, 0.2, 1.0;
  selectedFaceGlyphColor<<223.0/255.0, 210.0/255.0, 16.0/255.0;
  selectedVectorGlyphColor<<0.0,1.0,0.5;
  
  //triangle mesh setup
  viewer.data().set_face_based(true);
  viewer.data().set_mesh(VMesh, FMesh);
  update_triangle_mesh();
  //apending and updating raw field mesh
  viewer.append_mesh();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
