#include <iostream>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
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
#include <igl/edge_topology.h>


Eigen::VectorXi cIDs, matching, indices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V, rawField,representative, cValues;
Eigen::MatrixXcd powerField;
igl::viewer::Viewer viewer;

Eigen::MatrixXd positiveIndexColors(4, 3), negativeIndexColors(4, 3);

int N = 5;
bool normalized = false;
bool onePressed = false;

void update_mesh()
{
  
  Eigen::MatrixXd fullC(F.rows(),3);
  fullC = Eigen::RowVector3d(1, 1, 1).replicate(F.rows(), 1);
  
  for (int i = 0; i < cIDs.rows(); i++)
    fullC.row(cIDs(i)) = Eigen::RowVector3d(1, 0, 0);
  
  Eigen::MatrixXd fullV=V;
  Eigen::MatrixXi fullF=F;
  
  directional::power_field(V, F, cIDs, cValues, N, powerField);
  directional::power_to_representative(V, F, powerField, N, representative);
  
  if (normalized)
    representative.rowwise().normalize();
  
  directional::representative_to_raw(V,F,representative, N, rawField);
  
  if (cIDs.rows()!=0){
    directional::principal_matching(V, F, EV, EF, FE, rawField, matching, effort);
    
    directional::effort_to_indices(V,F,EV, EF, effort,N, indices);
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
    
      // Reset the constraints
    case 'R':
      cIDs.resize(0);
      cValues.resize(0, 6);
      update_mesh();
      break;
      
      // Toggle normalization
    case 'N':
      normalized = !normalized;
      update_mesh();
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
  if ((key != 0 && key != 2) || !onePressed)
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
      update_mesh();
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
    update_mesh();
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
  "  1+R-bttn  Remove constraint" << std::endl;
  
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
  igl::edge_topology(V, F, EV,FE,EF);
  
  cIDs.resize(0);
  cValues.resize(0, 3);
  
  update_mesh();
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
