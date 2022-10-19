#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/principal_matching.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>

Eigen::VectorXi constFaces;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField pvFieldHard, pvFieldSoft, rawFieldHard, rawFieldSoft,constraintsField;
Eigen::MatrixXd constVectors;
Eigen::VectorXd alignWeights;

double smoothWeight, roSyWeight;

directional::DirectionalViewer viewer;

int N = 4;

typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;

bool alterRoSyWeight=false;

void recompute_field()
{
  directional::polyvector_field(ftb, constFaces, constVectors, smoothWeight, roSyWeight, Eigen::VectorXd::Constant(constFaces.size(),-1), N, pvFieldHard);
  directional::polyvector_field(ftb, constFaces, constVectors, smoothWeight, roSyWeight, alignWeights, N, pvFieldSoft);
}

void update_visualization()
{
  viewer.set_field(constraintsField,Eigen::MatrixXd(), 1,0.9,0,0.2);
  viewer.set_selected_faces(constFaces,1);
  viewer.toggle_field(true,1);
  viewer.toggle_mesh(false,0);
  if (viewingMode==CONSTRAINTS)
    viewer.toggle_field(false,0);
    
  if (viewingMode==HARD_PRESCRIPTION){
    directional::polyvector_to_raw(pvFieldHard, rawFieldHard, N%2==0);
    directional::principal_matching(rawFieldHard);
    viewer.set_field(rawFieldHard,Eigen::MatrixXd(), 0,0.9,0,2.0);
    viewer.toggle_field(true,0);
  }
  
  if (viewingMode==SOFT_PRESCRIPTION){
    directional::polyvector_to_raw(pvFieldSoft, rawFieldSoft, N%2==0);
    directional::principal_matching(rawFieldSoft);
    viewer.set_field(rawFieldSoft,Eigen::MatrixXd(), 0,0.9,0,2.0);
    viewer.toggle_field(true,0);
  }
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1':
      alterRoSyWeight = !alterRoSyWeight;
      if (alterRoSyWeight)
        std::cout<<"Altering RoSy weight"<<std::endl;
      else
        std::cout<<"Altering alignment weight"<<std::endl;
      break;
      
    case 'P':
      if (alterRoSyWeight){
        roSyWeight*=2.0;
        std::cout<<"roSyWeight: "<<roSyWeight<<std::endl;
      }else{
        alignWeights.array()*=2.0;
        std::cout<<"alignWeights: "<<alignWeights(1)<<std::endl;
      }
      recompute_field();
      break;
      
    case 'M':
      if (alterRoSyWeight){
        roSyWeight/=2.0;
        std::cout<<"roSyWeight: "<<roSyWeight<<std::endl;
      }else{
        alignWeights.array()/=2.0;
        std::cout<<"alignWeights: "<<alignWeights(1)<<std::endl;
      }
      recompute_field();
      break;
      
    case '2':
      viewingMode=CONSTRAINTS;
      break;
      
    case '3':
      viewingMode=HARD_PRESCRIPTION;
      break;
      
    case '4':
      viewingMode=SOFT_PRESCRIPTION;
      break;
    
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/fandisk.rawfield", (viewingMode==HARD_PRESCRIPTION ? rawFieldHard : rawFieldSoft)))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. " << std::endl;
      break;
      
      
  }
  
  update_visualization();
  return true;
}


int main()
{
  
  std::cout <<
  "  1    Choose either RoSy or Alignment weight altering (default: Alignment)" << std::endl <<
  "  P/M  Increase/Decrease (RoSy or Alignment) weight" << std::endl<<
  "  2    Show only constraints" << std::endl <<
  "  3    Show fixed-alignment result" << std::endl <<
  "  4    Shoe soft-alignment result" << std::endl <<
  "  W    Save field"<< std::endl;
  
  // Load mesh
  directional::readOFF(TUTORIAL_SHARED_PATH "/fandisk.off", mesh);
  ftb.init(mesh);
  pvFieldHard.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
  pvFieldSoft.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
  
  //discovering and constraining sharp edges
  std::vector<int> constFaceslist;
  std::vector<Eigen::Vector3d> constVectorslist;
  for (int i=0;i<mesh.EF.rows();i++){
    if (mesh.faceNormals.row(mesh.EF(i,0)).dot(mesh.faceNormals.row(mesh.EF(i,1)))<0.5){
      constFaceslist.push_back(mesh.EF(i,0));
      constFaceslist.push_back(mesh.EF(i,1));
      constVectorslist.push_back((mesh.V.row(mesh.EV(i,0))-mesh.V.row(mesh.EV(i,1))).normalized());
      constVectorslist.push_back((mesh.V.row(mesh.EV(i,0))-mesh.V.row(mesh.EV(i,1))).normalized());
    }
  }
  
  constFaces.resize(constFaceslist.size());
  constVectors.resize(constVectorslist.size(),3);
  for (int i=0;i<constFaces.size();i++){
    constFaces(i)=constFaceslist[i];
    constVectors.row(i)=constVectorslist[i];
  }
                                
  //generating the viewing fields
  Eigen::MatrixXd rawFieldConstraints=Eigen::MatrixXd::Zero(mesh.F.rows(),N*3);
  Eigen::VectorXi posInFace=Eigen::VectorXi::Zero(mesh.F.rows());
  for (int i=0;i<constFaces.size();i++){
    rawFieldConstraints.block(constFaces(i),3*posInFace(constFaces(i)), 1,3)=constVectors.row(i);
    posInFace(constFaces(i))++;
  }
  
  //Just to show the other direction if N is even, since we are by default constraining it
  if (N%2==0)
    rawFieldConstraints.middleCols(rawFieldConstraints.cols()/2, rawFieldConstraints.cols()/2)=-rawFieldConstraints.middleCols(0, rawFieldConstraints.cols()/2);
  
  
  constraintsField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
  constraintsField.set_extrinsic_field(rawFieldConstraints);
  
  smoothWeight = 1.0;
  roSyWeight = 1.0;
  alignWeights = Eigen::VectorXd::Constant(constFaces.size(),1.0);
  
  //triangle mesh setup
  viewer.set_mesh(mesh,0);
  //ghost mesh only for constraints
  viewer.set_mesh(mesh, 1);
  recompute_field();
  update_visualization();
  
  viewer.callback_key_down = &key_down;

  viewer.launch();
}
