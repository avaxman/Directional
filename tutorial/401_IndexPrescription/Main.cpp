#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
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
#include <directional/singularity_spheres.h>
#include <directional/glyph_lines_raw.h>
#include <directional/write_raw_field.h>


Eigen::VectorXi cycleIndices;
Eigen::VectorXd cycleCurvature;
Eigen::SparseMatrix<double> basisCycles;
Eigen::VectorXi vertex2cycle, innerEdges;
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltSolver;

Eigen::RowVector3d rawGlyphColor;
Eigen::MatrixXi FMesh, FSings, FField, EV, FE, EF;
Eigen::MatrixXd VMesh, VSings, VField, BC, FN;
Eigen::MatrixXd CMesh, CSings, CField;
Eigen::MatrixXd rawField;
Eigen::VectorXd rotationField;
std::vector<std::vector<int>> cycleFaces;
int currCycle=0;

igl::opengl::glfw::Viewer viewer;

int eulerChar, numGenerators, numBoundaries;

int N = 2;

bool drag = false;
bool select = false;

double globalRotation=0;


void update_triangle_mesh()
{
  
  CMesh=Eigen::MatrixXd::Constant(FMesh.rows(), 3, 1.0);
  
  for (int i=0;i<cycleFaces[currCycle].size();i++)
    CMesh.row(cycleFaces[currCycle][i])<<1.0,0.0,0.0;
  
  viewer.data_list[0].set_colors(CMesh);
}

void update_raw_field_mesh()
{
  using namespace Eigen;
  VectorXd rotationAngles;
  double linfError;
  
  int sum = round(cycleIndices.head(cycleIndices.size() - numGenerators).sum());
  if (eulerChar*N != sum)
  {
    std::cout << "Warning: All non-generator singularities should add up to N * the Euler characteristic."<<std::endl;
    std::cout << "Total indices: " << sum << std::endl;
    std::cout << "Expected: " << eulerChar*N << std::endl;
  }
  
  directional::index_prescription(VMesh,FMesh,EV, innerEdges, basisCycles,cycleCurvature, cycleIndices,ldltSolver, N,rotationAngles, linfError);
  std::cout<<"field linfError: "<<linfError<<std::endl;
  
  Eigen::MatrixXd representative;
  directional::rotation_to_representative(VMesh, FMesh,EV,EF,rotationAngles,N,globalRotation, representative);
  directional::representative_to_raw(VMesh,FMesh,representative,N, rawField);
  
  directional::glyph_lines_raw(VMesh, FMesh, rawField, rawGlyphColor,  VField, FField, CField);
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
  viewer.data_list[1].set_colors(CField);
  
}

void update_singularities_mesh()
{
  Eigen::VectorXi singVertices, singIndices;
  std::vector<int> singVerticesList, singIndicesList;
  for (int i=0;i<VMesh.rows();i++)
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
  
  directional::singularity_spheres(VMesh, FMesh, singVertices, singIndices, directional::defaultSingularityColors(N), VSings, FSings, CSings);
  
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].show_faces = true;
  viewer.data_list[2].show_lines = false;
  viewer.data_list[2].set_colors(CSings);
}


bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': select=false; break;
  }
  return true;
}
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': select=true; break;
    case '1':
      globalRotation+=0.314;
      update_raw_field_mesh();
      break;
    case '-':
    case '_':
      cycleIndices(currCycle)--;
      update_raw_field_mesh();
      update_singularities_mesh();
      break;
    case '+':
    case '=':
      cycleIndices(currCycle)++;
      update_raw_field_mesh();
      update_singularities_mesh();
      break;

    case 'B':
      if (numBoundaries)
      {
        //Loop through the boundary cycles.
        if (currCycle >= basisCycles.rows()-numBoundaries-numGenerators && currCycle < basisCycles.rows()-numBoundaries-numGenerators)
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
  if ((key != 0)||(!select))
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
                               viewer.core.proj, viewer.core.viewport, VMesh, FMesh, fid, bc))
  {
    Eigen::Vector3d::Index maxCol;
    bc.maxCoeff(&maxCol);
    int currVertex=FMesh(fid, maxCol);
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
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV,FE,EF);
  
  directional::dual_cycles(VMesh, FMesh,EV, EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
  cycleIndices=Eigen::VectorXi::Constant(basisCycles.rows(),0);
  
  std::vector<std::vector<int>> boundaryLoops;
  igl::boundary_loop(FMesh, boundaryLoops);
  numBoundaries=boundaryLoops.size();
  eulerChar = VMesh.rows() - EV.rows() + FMesh.rows();
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
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(Eigen::RowVector3d::Constant(3,1.0));
  viewer.data().show_lines = false;
  
  //apending and updating raw field mesh
  viewer.append_mesh();
  
  //singularity mesh
  viewer.append_mesh();
  
  update_triangle_mesh();
  update_raw_field_mesh();
  update_singularities_mesh();
  viewer.selected_data_index=0;
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
