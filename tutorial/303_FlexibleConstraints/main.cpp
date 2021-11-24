#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/sharp_edges.h>
#include <igl/edge_topology.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

Eigen::VectorXi constFaces, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V;
Eigen::MatrixXd normals,constVectors;
Eigen::MatrixXd rawFieldConstraints, rawFieldHard, rawFieldSoft;
Eigen::VectorXd alignWeights;
Eigen::MatrixXcd pvFieldHard, pvFieldSoft;
double smoothWeight, roSyWeight;

directional::DirectionalViewer viewer;

int N = 4;

typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;

bool normalized=false;

void recompute_field()
{
  directional::polyvector_field(V, F, constFaces, constVectors, smoothWeight, roSyWeight, Eigen::VectorXd::Constant(constFaces.size(),-1), N, pvFieldHard);
  //std::cout<<"Done computing field!"<<std::endl;
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
    //std::cout<<"before polyvector_to_raw()"<<std::endl;
    //std::cout<<"pvFieldHard.rows(): "<<pvFieldHard.rows()<<std::endl;
    //std::cout<<"pvFieldHard.cols(): "<<pvFieldHard.cols()<<std::endl;
    //std::cout<<"F.rows(): "<<F.rows()<<std::endl;
    directional::polyvector_to_raw(V, F, pvFieldHard, N, rawFieldHard, N%2==0);
    //std::cout<<"after polyvector_to_raw()"<<std::endl;
    if (normalized)
      for(int n = 0; n < N; n++)
        rawFieldHard.middleCols(n*3, 3).rowwise().normalize();
    
    //std::cout<<"after normalize()"<<std::endl;
    
    directional::principal_matching(V, F, EV, EF, FE, rawFieldHard, matching, effort, singVertices, singIndices);
    //std::cout<<"after principal_matching()"<<std::endl;
    viewer.set_field(rawFieldHard);
    viewer.set_singularities(singVertices, singIndices);
  }
  
  if (viewingMode==SOFT_PRESCRIPTION){
    viewer.set_selected_faces(Eigen::VectorXi());
    directional::polyvector_to_raw(V, F, pvFieldSoft, N, rawFieldSoft, N%2==0);
    if (normalized)
      for(int n = 0; n < N; n++)
        rawFieldSoft.middleCols(n*3, 3).rowwise().normalize();
    
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
      roSyWeight*=2.0;
      std::cout<<"roSyWeight: "<<roSyWeight<<std::endl;
      recompute_field();
      break;
    case '2':
      roSyWeight/=2.0;
      roSyWeight = (roSyWeight < 0.0 ? 0.0 : roSyWeight);
      std::cout<<"roSyWeight: "<<roSyWeight<<std::endl;
      recompute_field();
      break;
      
    case '3':
      alignWeights.array()*=2.0;
      std::cout<<"alignWeights: "<<alignWeights(1)<<std::endl;
      recompute_field();
      break;
      
    case '4':
      alignWeights.array()/=2.0;
      alignWeights = (alignWeights(1) < 0.0 ? Eigen::VectorXd::Zero(alignWeights.size()) : alignWeights);
      std::cout<<"alignWeights: "<<alignWeights(1)<<std::endl;
      recompute_field();
      break;
      
    case '5':
      viewingMode=CONSTRAINTS;
      break;
      
    case '6':
      viewingMode=SOFT_PRESCRIPTION;
      break;
      
    case '7':
      viewingMode=HARD_PRESCRIPTION;
      break;
    
    case 'N':
      normalized = !normalized;
      break;
      
    /*case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/fandisk.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. " << std::endl;
      break;*/
      
      
  }
  
  update_visualization();
  return true;
}


int main()
{
  
  std::cout <<
  "  0+L-bttn   Choose face" << std::endl <<
  "  0+R-bttn   Edit vector in current face" << std::endl<<
  "  1          Choose vector in current face." << std::endl <<
  "  R          Reset the constraints" << std::endl <<
  "  N          Toggle field normalization" << std::endl;
  
  // Load mesh
  igl::readOBJ(TUTORIAL_SHARED_PATH "/spherers.obj", V, F);
  igl::edge_topology(V, F, EV, FE, EF);
  
  igl::per_face_normals(V, F, normals);
  
  //discovering sharp edges
  /*std::vector<int> blist;
  std::vector<Eigen::Vector3d> bclist;
  for (int i=0;i<EF.rows();i++){
    if (normals.row(EF(i,0)).dot(normals.row(EF(i,1)))<0.5){
      blist.push_back(EF(i,0));
      blist.push_back(EF(i,1));
      bclist.push_back((V.row(EV(i,0))-V.row(EV(i,1))).normalized());
      bclist.push_back((V.row(EV(i,0))-V.row(EV(i,1))).normalized());
    }
  }
  
  b.resize(blist.size());
  bc.resize(bclist.size(),3);
  for (int i=0;i<blist.size();i++){
    b(i)=blist[i];
    bc.row(i)=bclist[i];
  }
                                
  w=Eigen::VectorXd::Constant(F.rows(),1.0);  //Equal weight with the smoothness*/
  
  //Inserting a few hard constraints
  constFaces.resize(5);
  constFaces<<0,250,500, 1000,0;
  constVectors.resize(constFaces.rows(),3);
  for (int i=0;i<constFaces.size()-1;i++)
    constVectors.row(i)=(V.row(F(constFaces(i),1))-V.row(F(constFaces(i),0))).normalized();
  
  constVectors.row(constFaces.size()-1)=(V.row(F(constFaces(constFaces.size()-1),2))-V.row(F(constFaces(constFaces.size()-1),1))).normalized();
  
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
