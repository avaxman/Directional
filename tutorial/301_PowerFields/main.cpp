#include <iostream>
#include <Eigen/Core>
#include <igl/unproject_onto_mesh.h>
#include <directional/TriMesh.h>
#include <directional/readOBJ.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

directional::TriMesh mesh;
directional::IntrinsicVertexTangentBundle vtb;
directional::CartesianField rawField,powerFieldHard, powerFieldSoft;
Eigen::VectorXi constVertices;
Eigen::MatrixXd constVectors;
Eigen::VectorXd alignWeights;

directional::DirectionalViewer viewer;

int N = 4;
bool normalized = false;
bool zeroPressed = false;
bool viewFieldHard = true;

void recompute_field(){
  directional::power_field(vtb, constVertices, constVectors, Eigen::VectorXd::Constant(constVertices.size(),-1.0), N,powerFieldHard);
  directional::power_field(vtb, constVertices, constVectors, alignWeights, N,powerFieldSoft);
}

void update_visualization()
{
  directional::power_to_raw((viewFieldHard ? powerFieldHard : powerFieldSoft), N, rawField,normalized);
  
  directional::principal_matching(rawField);
  viewer.set_field(rawField,Eigen::MatrixXd(),0, 0.9, 0, 0.5);
  
  //Ghost mesh just showing field, to compare against constraints
  Eigen::MatrixXd constraintIntField = Eigen::MatrixXd::Zero(powerFieldHard.intField.rows(),2);
  for (int i=0;i<constVertices.size();i++)
    constraintIntField.row(constVertices(i))=powerFieldHard.intField.row(constVertices(i));
  
  directional::CartesianField constraintRawField, constraintPowerField;
  constraintPowerField.init(*(rawField.tb),POWER_FIELD,N);
  constraintPowerField.set_intrinsic_field(constraintIntField);
  
  directional::power_to_raw(constraintPowerField, N, constraintRawField);
  viewer.set_field(constraintRawField,Eigen::MatrixXd(), 1,0.9, 0, 0.3);
  viewer.toggle_mesh(false,0);
  viewer.toggle_field(true,0);
  viewer.toggle_field(true,1);
  viewer.set_selected_vertices(constVertices,1);

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
      constVertices.resize(0);
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
                                 viewer.core().proj, viewer.core().viewport, mesh.V, mesh.F, fid, baryInFace))
    {
        Eigen::Vector3d::Index maxCol;
        baryInFace.maxCoeff(&maxCol);
        int currVertex=mesh.F(fid, maxCol);
        int i;
        for (i = 0; i < constVertices.rows(); i++)
            if (constVertices(i) == currVertex)
                break;
        if (i == constVertices.rows())
        {
            constVertices.conservativeResize(constVertices.rows() + 1);
            constVertices(i) = currVertex;
            constVectors.conservativeResize(constVectors.rows() + 1, 3);
            alignWeights.conservativeResize(alignWeights.size()+1);
            if (alignWeights.size()==1)
                alignWeights(0)=1.0;
            else
                alignWeights(alignWeights.size()-1)=alignWeights(alignWeights.size()-2);
        }

        // Compute direction from the center of the face to the mouse
        //TODO: change way direction is chosen
        constVectors.row(i) =(mesh.V.row(mesh.F(fid, 0)) * baryInFace(0) +
                              mesh.V.row(mesh.F(fid, 1)) * baryInFace(1) +
                              mesh.V.row(mesh.F(fid, 2)) * baryInFace(2) -
                              mesh.V.row(currVertex)).normalized();
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
  directional::readOBJ(TUTORIAL_SHARED_PATH "/rocker-arm2500.obj", mesh);
  vtb.init(mesh);
  powerFieldHard.init(vtb, POWER_FIELD, N);
  powerFieldSoft.init(vtb, POWER_FIELD, N);
  constVertices.resize(0);
  constVectors.resize(0, 3);

  viewer.set_mesh(mesh,0);

  //ghost mesh only for constraints
  viewer.set_mesh(mesh, 1);
  recompute_field();
  update_visualization();

  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
