#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/mesh_function_isolines.h>
#include <directional/directional_viewer.h>
#include "polygonal_write_OFF.h"

#define NUM_N 1

int N[NUM_N];
int currN = 0;
directional::TriMesh meshWhole, meshCut[NUM_N];
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawField[NUM_N], combedField[NUM_N];
Eigen::VectorXi DPolyMesh[NUM_N];
Eigen::MatrixXi FPolyMesh[NUM_N];
Eigen::MatrixXd VPolyMesh[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
directional::DirectionalViewer viewer;

typedef enum {FIELD, INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;



int main()
{

  directional::readOFF(TUTORIAL_DATA_PATH "/vase.off",meshWhole);
  ftb.init(meshWhole);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-4.rawfield", ftb, N[0], rawField[0]);
  //directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-7.rawfield", ftb, N[1], rawField[1]);
  //directional::read_raw_field(TUTORIAL_DATA_PATH "/vase-11.rawfield", ftb, N[2], rawField[2]);

  bool verbose=false;
  
  //combing and cutting
  viewer.init();
  for (int i=0;i<NUM_N;i++){
    directional::principal_matching(rawField[i]);

    directional::IntegrationData intData(N[i]);
    std::cout<<"Setting up Integration #"<<i<<std::endl;
    directional::setup_integration(rawField[i], intData, meshCut[i], combedField[i]);
    
    intData.verbose=false;
    intData.integralSeamless=true;
    intData.roundSeams=false;
  
    std::cout<<"Solving integration for N="<<N[i]<<std::endl;
    directional::integrate(combedField[i],  intData, meshCut[i], NFunction[i],NCornerFunction[i]);
    
    std::cout<<"Done!"<<std::endl;
    
    //setting up mesh data from integration data
    directional::MeshFunctionIsolinesData mfiData;
    directional::setup_mesh_function_isolines(meshCut[i], intData, mfiData);
    
    //meshing and saving
    directional::mesh_function_isolines(meshWhole, mfiData,  verbose, VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);

    viewer.set_mesh(meshWhole);
    viewer.set_field(combedField[i], "", 0, i);
    viewer.set_seams(combedField[i].matching, 0, i);
    viewer.set_isolines(meshCut[i],NFunction[i],0, i);

    //Viewing polygonal mesh by direct call to Polyscope (Directional Viewer doesn't natively support polygonal meshes)
    std::vector<std::vector<int>> psF; psF.resize(DPolyMesh[i].size());
    for (int j=0;j<DPolyMesh[i].size();j++){
        psF[j].resize(DPolyMesh[i](j));
        for (int k=0;k<DPolyMesh[i](j);k++)
            psF[j][k]=FPolyMesh[i](j,k);
    }
    polyscope::registerSurfaceMesh("/vase-"+std::to_string(N[i])+"-generated.off", VPolyMesh[i], psF);

    //Also saving file
    hedra::polygonal_write_OFF(TUTORIAL_DATA_PATH "/vase-"+std::to_string(N[i])+"-generated.off", VPolyMesh[i], DPolyMesh[i], FPolyMesh[i]);

  }
  viewer.launch();
}


