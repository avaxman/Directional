#include <iostream>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
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

Eigen::VectorXi cIDs, matching, indices;
Eigen::VectorXd effort;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField,representative, cValues;
Eigen::MatrixXcd pvField;
igl::opengl::glfw::Viewer viewer;

int N = 3;

//User input variables
int cur = 0;
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
  
  // Compute the field
  directional::polyvector_field(VMesh, FMesh, cIDs, cValues, N, pvField);
  
  // Convert it so it can be drawn
  directional::polyvector_to_raw(VMesh, FMesh, pvField, N, rawField);
  
  if (normalized)
    for(int n = 0; n < N; n++)
      rawField.middleCols(n*3, 3).rowwise().normalize();
  
  if (cIDs.rows()!=0){
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
    }
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings, FSings);
    viewer.data_list[2].set_colors(CSings);
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
  int borders;
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
    case '1':
      cur = (cur+1)%N;
      break;
      // Reset the constraints
    case 'R':
      cIDs.resize(0);
      cValues.resize(0, 6);
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
  
  ;
  
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
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view ,
                               viewer.core.proj, viewer.core.viewport, VMesh, FMesh, fid, bc))
  {
    //Remove constraint
    if (key == 2)
    {
      int i;
      for (i = 0; i < cIDs.rows(); i++)
        if (cIDs(i) == fid)
          break;
      if (i == cIDs.rows())
        return false;
      cIDs(i) = cIDs(cIDs.size()-1);
      cIDs.conservativeResize(cIDs.rows() - 1);
      cValues.row(i) = cValues.row(cValues.rows() - 1);
      cValues.conservativeResize(cValues.rows() - 1, 3 * N);
      update_triangle_mesh();
      update_raw_field_mesh();
      return true;
    }
    
    if (key == 0)
    {
      int i;
      for (i = 0; i < cIDs.rows(); i++)
        if (cIDs(i) == fid)
          break;
      
      // Calculate direction from the center of the face to the mouse
      Eigen::RowVector3d rep =
      (VMesh.row(FMesh(fid, 0)) * bc(0) +
       VMesh.row(FMesh(fid, 1)) * bc(1) +
       VMesh.row(FMesh(fid, 2)) * bc(2) -
       (VMesh.row(FMesh(fid, 0)) +
        VMesh.row(FMesh(fid, 1)) +
        VMesh.row(FMesh(fid, 2))) / 3).normalized();
      
      // Add new entry
      if (i == cIDs.rows())
      {
        cIDs.conservativeResize(cIDs.rows() + 1);
        cIDs(i) = fid;
        cValues.conservativeResize(cValues.rows() + 1, 3 * N);
        
        //Create n-rosy for initial constraint
        Eigen::MatrixXd raw;
        Eigen::MatrixXd norm = Eigen::RowVector3d(VMesh.row(FMesh(fid, 1)) - VMesh.row(FMesh(fid, 0))).cross(Eigen::RowVector3d(VMesh.row(FMesh(fid, 2)) - VMesh.row(FMesh(fid, 0)))).normalized();
        directional::representative_to_raw(norm, rep, N, raw);
        
        // Rotate columns so first row is at current position and add them to the matrix
        cValues.block(i,0,1,N * 3 - cur * 3)=raw.block(0,cur*3,1,N * 3 - cur * 3);
        cValues.block(i,N * 3 - cur * 3, 1, cur*3)=raw.block(0,0,1,cur * 3);
        update_triangle_mesh();
        update_raw_field_mesh();
        return true;
      }
      
      // Calculate direction from the center of the face to the mouse
      cValues.block<1, 3>(i, cur * 3) = rep;
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
  "  R       Reset the constraints" << std::endl <<
  "  N       Toggle field normalization" << std::endl <<
  "  0+L-bttn  Place constraint pointing from the center of face to the cursor" << std::endl <<
  "  0+R-bttn  Remove constraint" << std::endl <<
  "  1      Toggle specific vector in face." << std::endl <<
  
  // Load mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/fandisk.off", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  cIDs.resize(0);
  cValues.resize(0, 3*N);
  
  //triangle mesh setup
  viewer.data_list[0].set_mesh(VMesh, FMesh);
  viewer.data_list[0].set_colors(Eigen::RowVector3d::Constant(3,1.0));
  
  viewer.selected_data_index = 0;
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
