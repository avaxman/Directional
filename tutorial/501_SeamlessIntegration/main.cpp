#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOFF.h>
#include <directional/writeOBJ.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/directional_viewer.h>


int N;
directional::TriMesh meshWhole, meshCut;
directional::PCFaceTangentBundle ftb;
directional::CartesianField rawField, combedField;
Eigen::MatrixXd cutUVFull, cutUVRot, cornerWholeUV;
directional::DirectionalViewer viewer;

typedef enum {FIELD, ROT_INTEGRATION, FULL_INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;

void callbackFunc(){
    ImGui::PushItemWidth(100);

    if (ImGui::Button("Save Textured OBJ File")){
        Eigen::MatrixXd emptyMat;
        directional::writeOBJ(TUTORIAL_DATA_PATH "/horsers-param-rot-seamless.obj", meshCut, cutUVRot, meshCut.F, "squares.mtl", "material_0");
        directional::writeOBJ(TUTORIAL_DATA_PATH "/horsers-param-full-seamless.obj", meshCut, cutUVFull, meshCut.F, "squares.mtl", "material_0");
    }

    ImGui::PopItemWidth();
}


int main()
{
    //setup_line_texture();
    directional::readOFF(TUTORIAL_DATA_PATH "/horsers.off", meshWhole);
    ftb.init(meshWhole);
    directional::read_raw_field(TUTORIAL_DATA_PATH "/horsers-cf.rawfield", ftb, N, rawField);

    //combing and cutting
    Eigen::VectorXd curlNorm;
    directional::curl_matching(rawField, curlNorm);
    std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;

    directional::IntegrationData intData(N);
    //cut_mesh_with_singularities(meshWhole, rawField.singLocalCycles, intData.face2cut);

    std::cout<<"Setting up Integration"<<std::endl;
    directional::setup_integration(rawField, intData, meshCut, combedField);
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


    intData.verbose=true;
    intData.integralSeamless=false;

    std::cout<<"Solving for permutationally-seamless integration"<<std::endl;
    directional::integrate(combedField, intData, meshCut, cutUVRot ,cornerWholeUV);
    //Extracting the UV from [U,V,-U, -V];
    cutUVRot=cutUVRot.block(0,0,cutUVRot.rows(),2);
    std::cout<<"Done!"<<std::endl;

    intData.integralSeamless = true;  //do not do translational seamless.
    std::cout<<"Solving for integrally-seamless integration"<<std::endl;
    directional::integrate(combedField,  intData, meshCut, cutUVFull,cornerWholeUV);
    cutUVFull=cutUVFull.block(0,0,cutUVFull.rows(),2);
    std::cout<<"Done!"<<std::endl;

    //viewer cut (texture) and whole (field) meshes
    viewer.init();
    viewer.set_surface_mesh(meshWhole, 0, "Original Mesh");
    viewer.set_surface_mesh(meshCut, 1, "Cut Mesh");
    viewer.set_cartesian_field(combedField);
    viewer.highlight_edges(seams, "Seams");
    viewer.set_uv(cutUVRot, "Rotationally Seamless UV", 1);
    viewer.set_uv(cutUVFull, "Full Seamless UV", 1);
    viewer.set_callback(callbackFunc);
    viewer.launch();
}


