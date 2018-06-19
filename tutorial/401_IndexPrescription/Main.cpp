#include <iostream>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
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
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, BC, FN;
Eigen::MatrixXd rawField;
Eigen::MatrixXd positiveIndexColors(4,3), negativeIndexColors(4,3);
Eigen::VectorXd rotationField;
std::vector<std::vector<int>> cycleFaces;
int currCycle=0;

igl::viewer::Viewer viewer;

int eulerChar, numGenerators, numBoundaries;

int N = 4;

bool drag = false;
bool select = false;

double globalRotation=0;

void update_mesh()
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
  
  directional::index_prescription(V,F,EV, innerEdges, basisCycles,cycleCurvature, cycleIndices,ldltSolver, N,rotationAngles, linfError);
  std::cout<<"field linfError: "<<linfError<<std::endl;
  
  Eigen::MatrixXd representative;
  directional::rotation_to_representative(V, F,EV,EF,rotationAngles,N,globalRotation, representative);
  directional::representative_to_raw(V,F,representative,N, rawField);
  
  Eigen::MatrixXd C(F.rows(),3);
  C.col(0)=Eigen::VectorXd::Constant(F.rows(),1.0);
  C.col(1)=Eigen::VectorXd::Constant(F.rows(),1.0);
  C.col(2)=Eigen::VectorXd::Constant(F.rows(),1.0);
  

  //cycle colors
  for (int i=0;i<cycleFaces[currCycle].size();i++)
    C.row(cycleFaces[currCycle][i])<<1.0,0.0,0.0;
  
  
  Eigen::MatrixXd fullV=V;
  Eigen::MatrixXi fullF=F;
  Eigen::MatrixXd fullC=C;
  
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
  
  directional::singularity_spheres(V, F, singVertices, singIndices, directional::defaultSingularityColors(N), false, true, fullV, fullF, fullC);
  
  directional::glyph_lines_raw(V, F, rawField, rawGlyphColor, false, true, fullV, fullF, fullC);
  
  viewer.data.clear();
  viewer.data.set_face_based(true);
  viewer.data.set_mesh(fullV, fullF);
  viewer.data.set_colors(fullC);
  
}


bool key_up(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '1': select=false; break;
  }
  return true;
}
bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '1': select=true; break;
    case '2':
      globalRotation+=0.314;
      update_mesh();
      break;
    case '-':
    case '_':
      cycleIndices(currCycle)--;
      update_mesh();
      break;
    case '+':
    case '=':
      cycleIndices(currCycle)++;
      update_mesh();
      break;

    case 'B':
      if (numBoundaries)
      {
        //Loop through the boundary cycles.
        if (currCycle >= basisCycles.rows()-numBoundaries-numGenerators && currCycle < basisCycles.rows()-numBoundaries-numGenerators)
          currCycle++;
        else
          currCycle = basisCycles.rows()-numBoundaries-numGenerators;
        update_mesh();
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
        update_mesh();
      }
      break;
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/torus.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. Error: " << errno << std::endl;
      
      /*if (directional::write_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singIndices, singPositions))
        std::cout << "Saved singularities" << std::endl;
      else
        std::cout << "Unable to save singularities. Error: " << errno << std::endl;
      break;
    case 'R':
      double x;
      //directional::read_trivial_field("../../data/field/trivial", meshV, meshF, indices, N, x);
      update_mesh();
      calculate_field();
      draw_field();
      break;*/
  }
  return true;
}


bool mouse_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  if ((key != 0)||(!select))
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
                               viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
  {
    Eigen::Vector3d::Index maxCol;
    bc.maxCoeff(&maxCol);
    int currVertex=F(fid, maxCol);
    currCycle=vertex2cycle(currVertex);
    update_mesh();
    return true;
  }
  return false;
};


int main()
{
  
  
  std::cout <<
  "  The field will appear if indices are correct " << std::endl<<
  "  1+ L-bttn  Select vertex cycle" << std::endl <<
  "  B          Loop through boundary cycles" << std::endl <<
  "  G          Loop through generator cycles" << std::endl <<
  "  +          Increase index of current cycle" << std::endl <<
  "  -          Decrease index  of current cycle" << std::endl <<
  "  2          rotate field globally" << std::endl;
  
  igl::readOBJ(TUTORIAL_SHARED_PATH "/torus.obj", V, F);
  igl::edge_topology(V, F, EV,FE,EF);
  
  directional::dual_cycles(V, F,EV, EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
  cycleIndices=Eigen::VectorXi::Constant(basisCycles.rows(),0);
  
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
     //std::cout<<"k: "<<k<< std::endl;
    for (Eigen::SparseMatrix<double>::InnerIterator it(basisCycles,k); it; ++it){
      int f1=EF(innerEdges(it.col()),0);
      int f2=EF(innerEdges(it.col()),1);
      //std::cout<<"it.col():"<<it.col()<<std::endl;
      //std::cout<<"it.row():"<<it.row()<<std::endl;
      //std::cout<<"f1, f2: "<<f1<<","<<f2<<std::endl;
      if (f1!=-1)
        cycleFaces[it.row()].push_back(f1);
      if (f2!=-1)
        cycleFaces[it.row()].push_back(f2);
    }
  }
  
  update_mesh();
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
