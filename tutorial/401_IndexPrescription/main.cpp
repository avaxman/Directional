#include <iostream>
#include <Eigen/Core>
#include <igl/unproject_onto_mesh.h>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/index_prescription.h>
#include <directional/rotation_to_raw.h>
#include <directional/write_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/directional_viewer.h>


directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field;
Eigen::VectorXi cycleIndices, presSingVertices, presSingIndices;
std::vector<Eigen::VectorXi> cycleFaces;
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltSolver;
int currCycle=0;

directional::DirectionalViewer viewer;

int N = 2;

bool drag = false;
bool _select = false;

double globalRotation=0;

void update_raw_field()
{
  using namespace Eigen;
  VectorXd rotationAngles;
  double linfError;
  
  int sum = round(cycleIndices.head(cycleIndices.size() - mesh.numGenerators).sum());
  if (mesh.eulerChar*N != sum)
  {
    std::cout << "Warning: All non-generator singularities should add up to N * the Euler characteristic."<<std::endl;
    std::cout << "Total indices: " << sum <<"/"<<N<< std::endl;
    std::cout << "Expected: " << mesh.eulerChar*N<<"/"<<N<< std::endl;
  }
  
  directional::index_prescription(cycleIndices, N,globalRotation, ldltSolver, field, rotationAngles, linfError);
  std::cout<<"Index prescription linfError: "<<linfError<<std::endl;
  
  viewer.set_field(field);
  viewer.set_selected_faces(cycleFaces[currCycle]);

}

void update_singularities()
{
  Eigen::VectorXi singVertices, singIndices;
  std::vector<int> singVerticesList, singIndicesList;
  for (int i=0;i<field.tb->local2Cycle.rows();i++)
    if (cycleIndices(field.tb->local2Cycle(i))){
      singVerticesList.push_back(i);
      singIndicesList.push_back(cycleIndices(field.tb->local2Cycle(i)));
    }
  
  singVertices.resize(singVerticesList.size());
  singIndices.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVertices(i)=singVerticesList[i];
    singIndices(i)=singIndicesList[i];
  }
  
  viewer.set_singularities(singVertices, singIndices);
}


bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': _select=false; break;
  }
  return true;
}
bool key_down(igl::opengl::glfw::Viewer& _viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': _select=true; break;
    case '1':
      globalRotation+=0.314;
      update_raw_field();
      break;
    case '-':
    case '_':
      cycleIndices(currCycle)--;
      update_raw_field();
      update_singularities();
      break;
    case '+':
    case '=':
      cycleIndices(currCycle)++;
      update_raw_field();
      update_singularities();
      break;

    case 'B':
      if (mesh.boundaryLoops.size())
      {
        //Loop through the boundary cycles.
        if (currCycle >= field.tb->cycles.rows()-mesh.boundaryLoops.size()-mesh.numGenerators && currCycle < field.tb->cycles.rows()-mesh.numGenerators-1)
          currCycle++;
        else
          currCycle = field.tb->cycles.rows()-mesh.boundaryLoops.size()-mesh.numGenerators;
          viewer.set_selected_faces(cycleFaces[currCycle]);
      }
      break;
    case 'G':
      if (mesh.numGenerators)
      {
        //Loop through the generators cycles.
        if (currCycle >= field.tb->cycles.rows() - mesh.numGenerators && currCycle < field.tb->cycles.rows() - 1)
          currCycle++;
        else
          currCycle = field.tb->cycles.rows() - mesh.numGenerators;
        viewer.set_selected_faces(cycleFaces[currCycle]);
      }
      break;
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/fertility.rawfield", field))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. Error: " << errno << std::endl;
  }
  return true;
}


bool mouse_down(igl::opengl::glfw::Viewer& _viewer, int key, int modifiers)
{
  if ((key != 0)||(!_select))
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core().viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                               viewer.core().proj, viewer.core().viewport, mesh.V, mesh.F, fid, bc))
  {
    Eigen::Vector3d::Index maxCol;
    bc.maxCoeff(&maxCol);
    int currVertex=mesh.F(fid, maxCol);
    currCycle=field.tb->local2Cycle(currVertex);
    viewer.set_selected_faces(cycleFaces[currCycle]);
    return true;
  }
  return false;
};


int main()
{
  
  
  std::cout <<
  "  The field will appear if indices are correct " << std::endl<<
  "  0+ L-bttn  Select vertex cycle" << std::endl <<
  "  B          Loop through boundary cycles" << std::endl <<
  "  G          Loop through generator cycles" << std::endl <<
  "  +          Increase index of current cycle" << std::endl <<
  "  -          Decrease index  of current cycle" << std::endl <<
  "  1          rotate field globally" << std::endl;
  
  directional::readOBJ(TUTORIAL_SHARED_PATH "/fertility.obj",mesh);
  ftb.init(mesh);
  field.init(ftb, RAW_FIELD, N);

  cycleIndices=Eigen::VectorXi::Constant(field.tb->cycles.rows(),0);
  
  //loading singularities
  directional::read_singularities(TUTORIAL_SHARED_PATH "/fertility.sings", N,presSingVertices,presSingIndices);
  field.set_singularities(presSingVertices,presSingIndices);
  
  std::cout<<"Euler characteristic: "<<mesh.eulerChar<<std::endl;
  std::cout<<"#generators: "<<mesh.numGenerators<<std::endl;
  std::cout<<"#boundaries: "<<mesh.boundaryLoops.size()<<std::endl;
  
  //collecting cycle faces for visualization
  std::vector<std::vector<int> > cycleFaceList(field.tb->cycles.rows());
  for (int k=0; k<field.tb->cycles.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(field.tb->cycles,k); it; ++it){
      int f1=mesh.EF(mesh.innerEdges(it.col()),0);
      int f2=mesh.EF(mesh.innerEdges(it.col()),1);
      if (f1!=-1)
        cycleFaceList[it.row()].push_back(f1);
      if (f2!=-1)
        cycleFaceList[it.row()].push_back(f2);
    }
  }
  
  for (int i=0;i<presSingVertices.size();i++)
    cycleIndices(field.tb->local2Cycle(presSingVertices(i)))=presSingIndices(i);

  cycleFaces.resize(field.tb->cycles.rows());
  for (int i=0;i<cycleFaceList.size();i++)
    cycleFaces[i] = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(cycleFaceList[i].data(), cycleFaceList[i].size());
  
  //triangle mesh setup
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  viewer.set_selected_faces(cycleFaces[currCycle]);
  update_raw_field();
  update_singularities();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
