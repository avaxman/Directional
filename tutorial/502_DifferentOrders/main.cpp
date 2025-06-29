#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
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
directional::PCFaceTangentBundle ftb;
directional::CartesianField rawField[NUM_N], combedField[NUM_N];
Eigen::MatrixXd NFunction[NUM_N], NCornerFunction[NUM_N];
directional::DirectionalViewer viewer;

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
    viewer.set_surface_mesh(meshWhole);
    for (int i=0;i<NUM_N;i++){
        directional::principal_matching(rawField[i]);
        
        directional::IntegrationData intData(N[i]);
        std::cout<<"Setting up Integration N="<<N[i]<<std::endl;
        directional::setup_integration(rawField[i], intData, meshCut[i], combedField[i]);
        std::vector<int> seamsList;
        for (int i=0;i<meshWhole.F.rows();i++){
            //int hebegin = meshWhole.FH(i);
            //int heiterate = hebegin;
            for (int j=0;j<3;j++)
                if (intData.face2cut(i, j))
                    seamsList.push_back(meshWhole.FE(i,j));
            //heiterate = meshWhole.nextH(heiterate);
        }
        
        Eigen::VectorXi seams = Eigen::VectorXi::Map(seamsList.data(), seamsList.size());
        
        
        intData.verbose=false;
        intData.integralSeamless=true;
        intData.roundSeams=false;
        
        std::cout<<"Solving integration N=" << N[i]<<std::endl;
        directional::integrate(combedField[i],  intData, meshCut[i], NFunction[i],NCornerFunction[i]);
        
        std::cout<<"Done!"<<std::endl;
        viewer.set_cartesian_field(combedField[i], std::to_string(N[i]) + "-field", i);
        //viewer.highlight_edges(seams, "Seams", i);
        viewer.set_isolines(meshCut[i],NFunction[i],std::to_string(N[i]) + "-function",  i, 0.05);
        if (i!=0){
            viewer.toggle_cartesian_field(false, i);
            viewer.toggle_singularities(false, i);
            viewer.toggle_isolines(false, i);
        }
        //viewer.set_surface_mesh(meshCut[i], 1);
    }
    
    viewer.launch();
}


