#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/branched_isolines.h>
#include <directional/directional_viewer.h>

#define NUM_N 4

int N[NUM_N];
int currN = 0;
Eigen::MatrixXi FMeshWhole, FMeshCut[NUM_N];
Eigen::MatrixXd VMeshWhole, VMeshCut[NUM_N];
Eigen::MatrixXd rawField[NUM_N], combedField[NUM_N];
Eigen::VectorXd effort[NUM_N], combedEffort[NUM_N];
Eigen::VectorXi matching[NUM_N], combedMatching[NUM_N];
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices[NUM_N], singVertices[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
directional::DirectionalViewer viewer;



typedef enum {FIELD, INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_viewer()
{
  for (int i=0;i<NUM_N;i++){
    viewer.toggle_field(false,i);
    viewer.toggle_singularities(false,i);
    viewer.toggle_seams(false,i);
    viewer.toggle_isolines(false,i);
  }
  if (viewingMode==FIELD){
    viewer.toggle_field(true,currN);
    viewer.toggle_singularities(true,currN);
    viewer.toggle_seams(true,currN);
    viewer.toggle_isolines(false,currN);
  } else {
    viewer.toggle_field(false,currN);
    viewer.toggle_singularities(false,currN);
    viewer.toggle_seams(false,currN);
    viewer.toggle_isolines(true,currN);
  }
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = INTEGRATION; break;
    case '3': currN=(currN+1)%NUM_N; break;
  }
  update_viewer();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show isoline mesh" << std::endl <<
  "  3  change between different N" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/vase.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-2.rawfield", N[0], rawField[0]);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-4.rawfield", N[1], rawField[1]);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-7.rawfield", N[2], rawField[2]);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/vase-11.rawfield", N[3], rawField[3]);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
 
  //combing and cutting
  for (int i=0;i<NUM_N;i++){
    directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField[i], matching[i], effort[i],singVertices[i], singIndices[i]);
    
    
    directional::IntegrationData intData(N[i]);
    std::cout<<"Setting up Integration N="<<N[i]<<std::endl;
    directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField[i], matching[i], singVertices[i], intData, VMeshCut[i], FMeshCut[i], combedField[i], combedMatching[i]);
    
    intData.verbose=false;
    intData.integralSeamless=true;
    intData.roundSeams=false;
  
    std::cout<<"Solving integration N=" << N[i]<<std::endl;
    directional::integrate(VMeshWhole, FMeshWhole, FE, combedField[i],  intData, VMeshCut[i], FMeshCut[i], NFunction[i],NCornerFunction[i]);
    
    std::cout<<"Done!"<<std::endl;
    
    viewer.set_mesh(VMeshWhole, FMeshWhole,directional::DirectionalViewer::default_mesh_color(),i);
    viewer.set_field(combedField[i], directional::DirectionalViewer::indexed_glyph_colors(combedField[i]), i);
    viewer.set_singularities(singVertices[i], singIndices[i],i);
    viewer.set_seams(EV, FE, EF, combedMatching[i], i);
    viewer.set_isolines(VMeshCut[i], FMeshCut[i],NFunction[i],i);

    
  }
  
  update_viewer();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


