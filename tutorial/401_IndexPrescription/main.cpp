#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <directional/dual_cycles.h>
#include <directional/index_prescription.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/power_to_representative.h>
#include <directional/power_field.h>
#include <directional/write_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/directional_viewer.h>


Eigen::VectorXi cycleIndices;
Eigen::VectorXd cycleCurvature;
Eigen::SparseMatrix<double> basisCycles;
Eigen::VectorXi vertex2cycle, innerEdges;
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltSolver;

Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, BC, FN;
Eigen::MatrixXd CMesh;
Eigen::MatrixXd rawField;
Eigen::VectorXd rotationField;
std::vector<std::vector<int>> cycleFaces;
int currCycle=0;

directional::DirectionalViewer viewer;

int eulerChar, numGenerators, numBoundaries;

int N = 2;

bool drag = false;
bool _select = false;

double globalRotation=0;

void update_triangle_mesh()
{
  CMesh=directional::DirectionalViewer::default_mesh_color().replicate(F.rows(),1);
  
  for (int i=0;i<cycleFaces[currCycle].size();i++)
    CMesh.row(cycleFaces[currCycle][i])<<directional::DirectionalViewer::selected_face_color();
  
  viewer.set_mesh_colors(CMesh);
}

void update_raw_field()
{
  using namespace Eigen;
  VectorXd rotationAngles;
  double linfError;
  
  int sum = round(cycleIndices.head(cycleIndices.size() - numGenerators).sum());
  if (eulerChar*N != sum)
  {
    std::cout << "Warning: All non-generator singularities should add up to N * the Euler characteristic."<<std::endl;
    std::cout << "Total indices: " << sum <<"/"<<N<< std::endl;
    std::cout << "Expected: " << eulerChar*N<<"/"<<N<< std::endl;
  }
  
  directional::index_prescription(V,F,EV, innerEdges, basisCycles,cycleCurvature, cycleIndices,ldltSolver, N,rotationAngles, linfError);
  std::cout<<"Index prescription linfError: "<<linfError<<std::endl;
  
  Eigen::MatrixXd representative;
  directional::rotation_to_representative(V, F,EV,EF,rotationAngles,N,globalRotation, representative);
  directional::representative_to_raw(V,F,representative,N, rawField);
  
  viewer.set_field(rawField);
}

void update_singularities()
{
  Eigen::VectorXi singVertices, singIndices;
  std::vector<int> singVerticesList, singIndicesList;
  for (int i=0;i<V.rows();i++)
    if (cycleIndices(vertex2cycle(i))){
      singVerticesList.push_back(i);
      singIndicesList.push_back(cycleIndices(vertex2cycle(i)));
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
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
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
      if (numBoundaries)
      {
        //Loop through the boundary cycles.
        if (currCycle >= basisCycles.rows()-numBoundaries-numGenerators && currCycle < basisCycles.rows()-numGenerators-1)
          currCycle++;
        else
          currCycle = basisCycles.rows()-numBoundaries-numGenerators;
        update_triangle_mesh();
      }
      break;
    case 'G':
      if (numGenerators)
      {
        //Loop through the generators cycles.
        if (currCycle >= basisCycles.rows() - numGenerators && currCycle < basisCycles.rows() - 1)
          currCycle++;
        else
          currCycle = basisCycles.rows() - numGenerators;
        update_triangle_mesh();
      }
      break;
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/fertility.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. Error: " << errno << std::endl;
  }
  return true;
}


bool mouse_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  if ((key != 0)||(!_select))
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core().viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                               viewer.core().proj, viewer.core().viewport, V, F, fid, bc))
  {
    Eigen::Vector3d::Index maxCol;
    bc.maxCoeff(&maxCol);
    int currVertex=F(fid, maxCol);
    currCycle=vertex2cycle(currVertex);
    update_triangle_mesh();
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
  
  igl::readOBJ(TUTORIAL_SHARED_PATH "/fertility.obj", V, F);
  igl::edge_topology(V, F, EV,FE,EF);
  
  directional::dual_cycles(V, F,EV, EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
  cycleIndices=Eigen::VectorXi::Constant(basisCycles.rows(),0);
  
  //loading singularities
  Eigen::VectorXi singVertices, singIndices;
  directional::read_singularities(TUTORIAL_SHARED_PATH "/fertility.sings", N, singVertices, singIndices);
  
  for (int i=0;i<singVertices.size();i++)
    cycleIndices(vertex2cycle(singVertices(i)))=singIndices(i);
  
  std::vector<std::vector<int>> boundaryLoops;
  igl::boundary_loop(F, boundaryLoops);
  numBoundaries=boundaryLoops.size();
  eulerChar = V.rows() - EV.rows() + F.rows();
  numGenerators = 2 - eulerChar - boundaryLoops.size();
  
  std::cout<<"Euler characteristic: "<<eulerChar<<std::endl;
  std::cout<<"#generators: "<<numGenerators<<std::endl;
  std::cout<<"#boundaries: "<<numBoundaries<<std::endl;
  
  //collecting cycle faces for visualization
  cycleFaces.resize(basisCycles.rows());
  for (int k=0; k<basisCycles.outerSize(); ++k){
    for (Eigen::SparseMatrix<double>::InnerIterator it(basisCycles,k); it; ++it){
      int f1=EF(innerEdges(it.col()),0);
      int f2=EF(innerEdges(it.col()),1);
      if (f1!=-1)
        cycleFaces[it.row()].push_back(f1);
      if (f2!=-1)
        cycleFaces[it.row()].push_back(f2);
    }
  }
  
  //triangle mesh setup
  viewer.set_mesh(V, F);  
  update_triangle_mesh();
  update_raw_field();
  update_singularities();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
