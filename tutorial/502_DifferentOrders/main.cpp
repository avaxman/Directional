#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOFF.h>
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
directional::TriMesh meshWhole, meshCut[NUM_N];
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawField[NUM_N], combedField[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
directional::DirectionalViewer viewer;

//typedef enum {FIELD, INTEGRATION} ViewingModes;
//ViewingModes viewingMode=FIELD;

/*void update_viewer()
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
}*/


// Handle keyboard input
/*bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
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
}*/


int main()
{
  directional::readOFF(TUTORIAL_DATA_PATH "/vase.off", meshWhole);
  ftb.init(meshWhole);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-2.rawfield", ftb, N[0], rawField[0]);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-4.rawfield", ftb, N[1], rawField[1]);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-7.rawfield", ftb, N[2], rawField[2]);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-11.rawfield", ftb, N[3], rawField[3]);
  
  //combing and cutting
  viewer.init();
  for (int i=0;i<NUM_N;i++){
    directional::principal_matching(rawField[i]);

    directional::IntegrationData intData(N[i]);
    std::cout<<"Setting up Integration N="<<N[i]<<std::endl;
    directional::setup_integration(rawField[i], intData, meshCut[i], combedField[i]);
      Eigen::VectorXi seams = Eigen::VectorXi::Zero(meshWhole.EV.rows());
      for (int i=0;i<meshWhole.F.rows();i++)
          for (int j=0;j<3;j++)
              seams(meshWhole.HE(meshWhole.FH(i,j)))=intData.face2cut(i,j);
    
    intData.verbose=false;
    intData.integralSeamless=true;
    intData.roundSeams=false;
  
    std::cout<<"Solving integration N=" << N[i]<<std::endl;
    directional::integrate(combedField[i],  intData, meshCut[i], NFunction[i],NCornerFunction[i]);

    //Eigen::MatrixXd cutUVFull=NFunction[i].block(0,0,NFunction[i].rows(),2);
    
    std::cout<<"Done!"<<std::endl;
    
    viewer.set_mesh(meshWhole);
    viewer.set_field(combedField[i], "", 0, i);
    viewer.set_seams(seams, 0, i);
    viewer.set_isolines(meshCut[i],NFunction[i],0, i);

    viewer.set_mesh(meshCut[i], 1);
    //viewer.set_uv(cutUVFull, 1);
  }

  viewer.launch();
}


