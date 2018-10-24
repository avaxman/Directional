#include <iostream>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/glyph_lines_raw.h>
#include <directional/singularity_spheres.h>
#include <directional/write_raw_field.h>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>


Eigen::VectorXi cIDs, matching, indices;
Eigen::VectorXd effort;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField,representative, cValues;
Eigen::MatrixXcd powerField;
igl::opengl::glfw::Viewer viewer;

int N = 4;
bool normalized = false;
bool zeroPressed = false;

void update_triangle_mesh()
{
  
  CMesh=Eigen::MatrixXd::Constant(FMesh.rows(), 3, 1.0);
  for (int i = 0; i < cIDs.rows(); i++)
    CMesh.row(cIDs(i)) = Eigen::RowVector3d(0.5,0.1,0.1);
  
  viewer.data_list[0].set_colors(CMesh);
}

void update_raw_field_mesh()
{
  directional::power_field(VMesh, FMesh, cIDs, cValues, N, powerField);
  directional::power_to_representative(VMesh, FMesh, powerField, N, representative);
  
  if (normalized)
    representative.rowwise().normalize();
  
  directional::representative_to_raw(VMesh,FMesh,representative, N, rawField);
  
  //if (cIDs.rows()!=0){
    directional::principal_matching(VMesh, FMesh, EV, EF, FE, rawField, matching, effort);
    
    directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N, indices);
    std::vector<int> singIndicesList,singVerticesList;
    for (int i=0;i<VMesh.rows();i++)
      if (indices(i)!=0){
        singIndicesList.push_back(indices(i));
        singVerticesList.push_back(i);
      }
    
    Eigen::VectorXi singIndices(singIndicesList.size());
    Eigen::VectorXi singVertices(singVerticesList.size());
    for (int i=0;i<singIndicesList.size();i++){
      singIndices(i)=singIndicesList[i];
      singVertices(i)=singVerticesList[i];
    }
    
    directional::singularity_spheres(VMesh, FMesh, singVertices, singIndices, directional::defaultSingularityColors(N), VSings, FSings, CSings);
    
      directional::glyph_lines_raw(VMesh, FMesh, rawField, Eigen::RowVector3d(0, 0, 1),  VField, FField, CField);
    
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
  //}
  
}

bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': zeroPressed=false; break;
  }
  return true;
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Toggle field drawing for easier rotation
      
    case '0': zeroPressed=true; break;
      
      // Reset the constraints
    case 'R':
      cIDs.resize(0);
      cValues.resize(0, 6);
      update_triangle_mesh();
      update_raw_field_mesh();
      break;
      
      // Toggle normalization
    case 'N':
      normalized = !normalized;
      update_triangle_mesh();
      update_raw_field_mesh();
      break;
      
      
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/horsers.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. Error: " << errno << std::endl;
      
  }
  
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
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
    //Remove constraint
    if (key == 2)
    {
      if (cIDs.size()==0)  //nothing to remove
        return false;
      int i;
      for (i = 0; i < cIDs.rows(); i++)
        if (cIDs(i) == fid)
          break;
      if (i == cIDs.rows())
        return false;
      cIDs(i) = cIDs(cIDs.size() - 1);
      cIDs.conservativeResize(cIDs.rows() - 1);
      cValues.row(i) = cValues.row(cValues.rows() - 1);
      cValues.conservativeResize(cValues.rows() - 1, 3);
      update_triangle_mesh();
      update_raw_field_mesh();
      return true;
    }
    
    int i;
    for (i = 0; i < cIDs.rows(); i++)
      if (cIDs(i) == fid)
        break;
    if (i == cIDs.rows())
    {
      cIDs.conservativeResize(cIDs.rows() + 1);
      cIDs(i) = fid;
      cValues.conservativeResize(cValues.rows() + 1, 3);
    }
    
    // Compute direction from the center of the face to the mouse
    cValues.row(i) =
    (VMesh.row(FMesh(fid, 0)) * bc(0) +
     VMesh.row(FMesh(fid, 1)) * bc(1) +
     VMesh.row(FMesh(fid, 2)) * bc(2) -
     (VMesh.row(FMesh(fid, 0)) +
      VMesh.row(FMesh(fid, 1)) +
      VMesh.row(FMesh(fid, 2))) / 3).normalized();
    update_triangle_mesh();
    update_raw_field_mesh();
    return true;
  }
  return false;
};

int main()
{
  
  std::cout <<
  "  R        Reset the constraints" << std::endl <<
  "  N        Toggle field normalization" << std::endl <<
  "  W        Save raw field" << std::endl <<
  "  0+L-bttn Place constraint pointing from the center of face to the cursor" << std::endl <<
  "  0+R-bttn Remove constraint" << std::endl;
  
  // Load mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/horsers.off", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV,FE,EF);
  
  cIDs.resize(0);
  cValues.resize(0, 3);
  
  //triangle mesh setup
  viewer.data_list[0].set_mesh(VMesh, FMesh);
  viewer.data_list[0].set_colors(Eigen::RowVector3d::Constant(3,1.0));
  
  viewer.selected_data_index = 0;
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
