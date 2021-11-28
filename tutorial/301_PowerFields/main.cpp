#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>


Eigen::VectorXi constFaces, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField,representative, constVectors, barycenters;
Eigen::MatrixXcd powerFieldHard, powerFieldSoft;
directional::DirectionalViewer viewer;
Eigen::VectorXd alignWeights;

int N = 4;
bool normalized = false;
bool zeroPressed = false;
bool viewFieldHard = true;

void recompute_field(){
  directional::power_field(V, F, constFaces, constVectors, Eigen::VectorXd::Constant(constFaces.size(),-1.0), N, powerFieldHard);
  directional::power_field(V, F, constFaces, constVectors, alignWeights, N, powerFieldSoft);
}

void update_visualization()
{
  directional::power_to_representative(V, F, (viewFieldHard ? powerFieldHard : powerFieldSoft), N, representative);
  if (normalized)
    representative.rowwise().normalize();
  
  directional::representative_to_raw(V,F,representative, N, rawField);
  directional::principal_matching(V, F, EV, EF, FE, rawField, matching, effort, singVertices, singIndices);
  viewer.set_field(rawField,Eigen::MatrixXd(),0);
  viewer.set_singularities(singVertices, singIndices,0);
  
  //Ghost mesh just showing field, to compare against constraints
  Eigen::VectorXcd constraintField = Eigen::VectorXcd::Zero(powerFieldHard.rows());
  for (int i=0;i<constFaces.size();i++)
    constraintField(constFaces(i))=powerFieldHard(constFaces(i));
  directional::power_to_raw(V, F,constraintField, N, rawField);
  viewer.set_field(rawField,Eigen::MatrixXd(), 1);
  viewer.toggle_mesh(false,0);
  viewer.toggle_field(true,0);
  viewer.toggle_field(true,1);
  viewer.set_selected_faces(constFaces,1);
  
}

bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': zeroPressed=false; break;
  }
  return true;
}


bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directionalViewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
      // Toggle field drawing for easier rotation
      
    case '0': zeroPressed=true; break;
      
      // Reset the constraints
    case 'R':
      constFaces.resize(0);
      constVectors.resize(0, 3);
      recompute_field();
      break;
      
      // Toggle normalization
    case 'N':
      normalized = !normalized;
      break;
      
    case 'P':
      alignWeights.array()*=2.0;
      if (alignWeights.size()!=0)
        std::cout<<"alignWeights: "<<alignWeights(0)<<std::endl;
      recompute_field();
      break;
      
    case 'M':
      alignWeights.array()/=2.0;
      if (alignWeights.size()!=0)
        std::cout<<"alignWeights: "<<alignWeights(0)<<std::endl;
      recompute_field();
      break;

    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/rocker-arm2500.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. Error: " << errno << std::endl;
      
    case 'H':
      viewFieldHard = !viewFieldHard;
      if (viewFieldHard)
        std::cout<<"Viewing hard-constrained field"<<std::endl;
      else
        std::cout<<"Viewing Soft-constrained field"<<std::endl;
      break;
  }
  
  update_visualization();
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directionalViewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
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
   
    int i;
    for (i = 0; i < constFaces.rows(); i++)
      if (constFaces(i) == fid)
        break;
    if (i == constFaces.rows())
    {
      constFaces.conservativeResize(constFaces.rows() + 1);
      constFaces(i) = fid;
      constVectors.conservativeResize(constVectors.rows() + 1, 3);
      alignWeights.conservativeResize(alignWeights.size()+1);
      if (alignWeights.size()==1)
        alignWeights(0)=1.0;
      else
        alignWeights(alignWeights.size()-1)=alignWeights(alignWeights.size()-2);
    }
    
    // Compute direction from the center of the face to the mouse
    constVectors.row(i) =(V.row(F(fid, 0)) * baryInFace(0) +
                          V.row(F(fid, 1)) * baryInFace(1) +
                          V.row(F(fid, 2)) * baryInFace(2) - barycenters.row(fid)).normalized();
    recompute_field();
    update_visualization();
    return true;
  }
  return false;
};

int main()
{
  
  std::cout <<
  "  R        Reset the constraints" << std::endl <<
  "  N        Toggle field normalization" << std::endl <<
  "  P/M      Increase/Decrease soft-alignment weights"<< std::endl <<
  "  H        Toggle hard/soft alignment"<< std::endl <<
  "  W        Save raw field" << std::endl <<
  "  0+L-bttn Place constraint pointing from the center of face to the cursor" << std::endl;
  
  // Load mesh
  igl::readOBJ(TUTORIAL_SHARED_PATH "/rocker-arm2500.obj", V, F);
  igl::edge_topology(V, F, EV,FE,EF);
  igl::barycenter(V, F, barycenters);
  
  constFaces.resize(0);
  constVectors.resize(0, 3);
  
  viewer.set_mesh(V,F,0);
  
  //ghost mesh only for constraints
  viewer.set_mesh(V,F, 1);
  recompute_field();
  update_visualization();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
