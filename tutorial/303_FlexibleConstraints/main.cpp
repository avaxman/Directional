#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/local_basis.h>
#include <igl/sharp_edges.h>
#include <igl/edge_topology.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_to_raw_companion.h>
#include <directional/polyvector_field.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

Eigen::VectorXi constFaces, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V, B1, B2;
Eigen::MatrixXd normals,constVectors;
Eigen::MatrixXd rawFieldConstraints, rawFieldHard, rawFieldSoft;
Eigen::VectorXd alignWeights;
Eigen::MatrixXcd pvFieldHard, pvFieldSoft;
double smoothWeight, roSyWeight;

directional::DirectionalViewer viewer;

int N = 4;

typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;

bool alterRoSyWeight=false;

void recompute_field()
{
  directional::polyvector_field(V, F, constFaces, constVectors, smoothWeight, roSyWeight, Eigen::VectorXd::Constant(constFaces.size(),-1), N, pvFieldHard);

  directional::polyvector_field(V, F, constFaces, constVectors, smoothWeight, roSyWeight, alignWeights, N, pvFieldSoft);
}

void update_visualization()
{
  
  if (viewingMode==CONSTRAINTS){
    viewer.set_selected_faces(constFaces);
    viewer.set_field(rawFieldConstraints);
  }
  if (viewingMode==HARD_PRESCRIPTION){
    viewer.set_selected_faces(Eigen::VectorXi());
    directional::polyvector_to_raw(B1,B2, pvFieldHard, N, rawFieldHard, N%2==0);
   
    directional::principal_matching(V, F, EV, EF, FE, rawFieldHard, matching, effort, singVertices, singIndices);
    viewer.set_field(rawFieldHard);
    viewer.set_singularities(singVertices, singIndices);
  }
  
  if (viewingMode==SOFT_PRESCRIPTION){
    viewer.set_selected_faces(Eigen::VectorXi());
    directional::polyvector_to_raw(B1, B2, pvFieldSoft, N, rawFieldSoft);
    
    directional::principal_matching(V, F, EV, EF, FE, rawFieldSoft, matching, effort, singVertices, singIndices);
    viewer.set_field(rawFieldSoft);
    viewer.set_singularities(singVertices, singIndices);

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
  "  1    Choose either RoSy or Alignment weight altering (default: RoSy)" << std::endl <<
  "  +-   Alter (RoSy or Alignment) weight" << std::endl<<
  "  2    Show constraints" << std::endl <<
  "  3    Show fixed-alignment result" << std::endl <<
  "  4    Shoe soft-alignment result" << std::endl <<
  "  W    Save field"<< std::endl;
  
  // Load mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/fandisk.off", V, F);
  igl::edge_topology(V, F, EV, FE, EF);
  igl::local_basis(V, F, B1, B2, normals);
  
  //discovering and constraining sharp edges
  std::vector<int> constFaceslist;
  std::vector<Eigen::Vector3d> constVectorslist;
  for (int i=0;i<EF.rows();i++){
    if (normals.row(EF(i,0)).dot(normals.row(EF(i,1)))<0.5){
      constFaceslist.push_back(EF(i,0));
      constFaceslist.push_back(EF(i,1));
      constVectorslist.push_back((V.row(EV(i,0))-V.row(EV(i,1))).normalized());
      constVectorslist.push_back((V.row(EV(i,0))-V.row(EV(i,1))).normalized());
    }
  }
  
  constFaces.resize(constFaceslist.size());
  constVectors.resize(constVectorslist.size(),3);
  for (int i=0;i<constFaces.size();i++){
    constFaces(i)=constFaceslist[i];
    constVectors.row(i)=constVectorslist[i];
  }
                                
  //generating the viewing fields
  rawFieldConstraints=Eigen::MatrixXd::Zero(F.rows(),N*3);
  Eigen::VectorXi posInFace=Eigen::VectorXi::Zero(F.rows());
  for (int i=0;i<constFaces.size();i++){
    rawFieldConstraints.block(constFaces(i),3*posInFace(constFaces(i)), 1,3)=constVectors.row(i);
    posInFace(constFaces(i))++;
  }
  
  smoothWeight = 1.0;
  roSyWeight = 1.0;
  alignWeights = Eigen::VectorXd::Constant(constFaces.size(),1.0);
  
  //triangle mesh setup
  viewer.set_mesh(V, F);
  recompute_field();
  update_visualization();
  
  viewer.callback_key_down = &key_down;

  viewer.launch();
}
