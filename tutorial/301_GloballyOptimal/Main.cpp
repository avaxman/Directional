#include <iostream>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/dual_cycles.h>
#include <directional/get_indices.h>
#include <directional/glyph_lines_raw.h>
#include <directional/singularity_spheres.h>
//#include <directional/write_raw_field.h>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>


Eigen::VectorXi cIDs, matching;
Eigen::MatrixXi F;
Eigen::MatrixXd V, C, rawField,representative, cValues;
Eigen::MatrixXcd powerField;
igl::viewer::Viewer viewer;

Eigen::MatrixXd positiveIndexColors(4, 3), negativeIndexColors(4, 3);

int N = 4;
bool drag = false;
bool normalized = false;
bool showSing = false;
bool onePressed = false;

void draw_field()
{
  
  directional::power_field(V, F, cIDs, cValues, N, powerField);
  directional::power_to_representative(V, F, powerField, N, representative);
  
  
  if (normalized)
    representative.rowwise().normalize();
  
  directional::representative_to_raw(V,F,representative, N, rawField);
  
  C = Eigen::RowVector3d(1, 1, 1).replicate(F.rows(), 1);
  
  for (int i = 0; i < cIDs.rows(); i++)
    C.row(cIDs(i)) = Eigen::RowVector3d(1, 0, 0);
  
  Eigen::VectorXi indices;
  Eigen::VectorXd effort;
  
  Eigen::MatrixXd fullV=V;
  Eigen::MatrixXi fullF=F;
  Eigen::MatrixXd fullC=C;
  
  if (cIDs.rows()!=0){
    directional::principal_matching(V, F, representative, N, matching, effort);
    
    directional::get_indices(V,F,effort,N, indices);
    std::vector<int> singIndicesList,singPositionsList;
    for (int i=0;i<V.rows();i++)
      if (indices(i)!=0){
        singIndicesList.push_back(indices(i));
        singPositionsList.push_back(i);
      }
    
    Eigen::VectorXi singIndices(singIndicesList.size());
    Eigen::VectorXi singPositions(singPositionsList.size());
    for (int i=0;i<singIndicesList.size();i++){
      singIndices(i)=singIndicesList[i];
      singPositions(i)=singPositionsList[i];
    }
    
    directional::singularity_spheres(V, F, singPositions, singIndices, positiveIndexColors, negativeIndexColors, false, true, fullV, fullF, fullC);
  }
  directional::glyph_lines_raw(V, F, rawField, Eigen::RowVector3d(0, 0, 1), false, true, fullV, fullF, fullC);
  
  // Update the viewer
  viewer.data.clear();
  viewer.data.set_face_based(true);
  viewer.data.set_mesh(fullV, fullF);
  viewer.data.set_colors(fullC);
}

bool key_up(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '1': onePressed=false; break;
  }
  return true;
}


bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Toggle field drawing for easier rotation
      
    case '1': onePressed=true; break;
    case 'D':
      drag = !drag;
      break;
      
      // Toggle singularities
    case 'S':
      showSing = !showSing;
      draw_field();
      break;
      
      // Reset the constraints
    case 'R':
      cIDs.resize(0);
      cValues.resize(0, 6);
      draw_field();
      break;
      
      // Toggle normalization
    case 'N':
      normalized = !normalized;
      draw_field();
      break;
      
      /*case 'W':
       if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/torus.rawfield", rawField))
       std::cout << "Saved raw field" << std::endl;
       else
       std::cout << "Unable to save raw field. Error: " << errno << std::endl;
       break;*/
      
  }
  
  
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  if (drag || (key != 0 && key != 2) || !onePressed)
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
                               viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
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
      cIDs(i) = cIDs(cIDs.size() - 1);
      cIDs.conservativeResize(cIDs.rows() - 1);
      cValues.row(i) = cValues.row(cValues.rows() - 1);
      cValues.conservativeResize(cValues.rows() - 1, 3);
      draw_field();
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
    (V.row(F(fid, 0)) * bc(0) +
     V.row(F(fid, 1)) * bc(1) +
     V.row(F(fid, 2)) * bc(2) -
     (V.row(F(fid, 0)) +
      V.row(F(fid, 1)) +
      V.row(F(fid, 2))) / 3).normalized();
    draw_field();
    return true;
  }
  return false;
};

int main()
{
  
  std::cout <<
  "  R       Reset the constraints" << std::endl <<
  "  N       Toggle field normalization" << std::endl <<
  "  1+L-bttn  Place constraint pointing from the center of face to the cursor" << std::endl <<
  "  1+R-bttn  Remove constraint" << std::endl <<
  
  // Set colors for Singularities
  positiveIndexColors << .25, 0, 0,
  .5, 0, 0,
  .75, 0, 0,
  1, 0, 0;
  
  negativeIndexColors << 0, .25, 0,
  0, .5, 0,
  0, .75, 0,
  0, 1, 0;
  
  // Load mesh
  igl::readOBJ(TUTORIAL_SHARED_PATH "/torus.obj", V, F);
  
  cIDs.resize(0);
  cValues.resize(0, 3);
  
  draw_field();
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
