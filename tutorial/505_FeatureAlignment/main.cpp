#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/branched_isolines.h>
#include <directional/polyvector_field.h>
#include <directional/directional_viewer.h>

int N = 4;
directional::TriMesh meshWhole, meshCut;
directional::PCFaceTangentBundle ftb;
directional::CartesianField rawField, combedField, pvField;
Eigen::MatrixXd NFunctionSign, NFunctionNoAlign, NFunctionRot, NFunctionFull, NCornerFunc;
directional::DirectionalViewer viewer;

typedef enum {NO_ALIGNMENT, ALIGNMENT_ROTATIONAL, ALIGNMENT_FULL} ViewingModes;
ViewingModes viewingMode=NO_ALIGNMENT;


void callbackFunc(){
    
    const char* items[] = {"No Alignment", "Alignment+Rotational Seamless","Alignment+Full Seamless"};
    static const char* current_item = NULL;
    
    ImGui::PushItemWidth(300);
    static float combo_width = 0.0f;
    if (combo_width == 0.0f) {
        ImGuiStyle& style = ImGui::GetStyle();
        for (auto& item : items)
            combo_width = std::max(combo_width, ImGui::CalcTextSize(item).x);
        combo_width += style.FramePadding.x * 5.0 + ImGui::GetFontSize() + style.ItemInnerSpacing.x;
    }
    
    ImGui::PushItemWidth(combo_width);
    if (ImGui::BeginCombo("Viewing Mode", current_item))
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]);
            if (ImGui::Selectable(items[n], is_selected)) {
                switch (n) {
                    case 0:
                        viewingMode = NO_ALIGNMENT;
                        viewer.set_isolines(meshCut, NFunctionNoAlign);
                        break;
                    case 1:
                        viewingMode = ALIGNMENT_ROTATIONAL;
                        viewer.set_isolines(meshCut, NFunctionRot);
                        break;
                    case 2:
                        viewingMode = ALIGNMENT_FULL;
                        viewer.set_isolines(meshCut, NFunctionFull);
                        break;
                }
                
                current_item = items[n];
            }
            if (is_selected)
                ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }
    
    ImGui::PopItemWidth();
}

int main()
{
    
    directional::readOFF(TUTORIAL_DATA_PATH "/cube_oversampled.off", meshWhole);
    ftb.init(meshWhole);
    
    std::vector<int> constIndicesList;
    std::vector<Eigen::RowVector3d> constVectorsList;
    std::vector<int> featureIndicesList;
    for (int i=0;i<meshWhole.EF.rows();i++){
        Eigen::RowVector3d nLeft = meshWhole.faceNormals.row(meshWhole.EF(i,0));
        Eigen::RowVector3d nRight = meshWhole.faceNormals.row(meshWhole.EF(i,1));
        Eigen::RowVector3d normEdgeVector = (meshWhole.V.row(meshWhole.F(meshWhole.EF(i,0),(meshWhole.EFi(i,0)+1)%3)) - meshWhole.V.row(meshWhole.F(meshWhole.EF(i,0),meshWhole.EFi(i,0)))).normalized();
        if (nLeft.dot(nRight)<=0.5){  //sharp edges
            constIndicesList.push_back(meshWhole.EF(i,0));
            constVectorsList.push_back(normEdgeVector);
            constIndicesList.push_back(meshWhole.EF(i,1));
            constVectorsList.push_back(normEdgeVector);
            featureIndicesList.push_back(i);
        }
    }
    Eigen::VectorXi constIndices(constIndicesList.size());
    Eigen::VectorXi featureIndices(featureIndicesList.size());
    Eigen::MatrixXd constVectors(constVectorsList.size(),3);
    for (int i=0;i<constIndicesList.size();i++){
        constIndices(i) = constIndicesList[i];
        constVectors.row(i) = constVectorsList[i];
    }
    for (int i=0;i<featureIndicesList.size();i++)
        featureIndices(i) = featureIndicesList[i];
    
    //Computing a field with the same const indices as the feature lines
    directional::PolyVectorData pvData;
    pvData.N = N;
    pvData.tb = &ftb;
    pvData.verbose = false;
    pvData.constSpaces = constIndices;
    pvData.constVectors = constVectors;
    pvData.wSmooth = 1.0;
    pvData.wRoSy = 1.0;
    pvData.wAlignment =Eigen::VectorXd::Constant(constIndices.size(),-1);
    directional::polyvector_field(pvData, pvField);
    directional::polyvector_to_raw(pvField, rawField);
    directional::principal_matching(rawField);
    
    directional::IntegrationData intDataNoAlign(N);
    intDataNoAlign.verbose = true;
    intDataNoAlign.integralSeamless = true;
    std::cout<<"Doing no-alignment integration"<<std::endl;
    directional::setup_integration(rawField, intDataNoAlign, meshCut, combedField);
    directional::integrate(combedField, intDataNoAlign, meshCut, NFunctionNoAlign, NCornerFunc);
    std::cout<<"Done with no-alignment integration"<<std::endl;
    
    directional::IntegrationData intDataRot(N);
    intDataRot.verbose = true;
    intDataRot.integralSeamless = false;
    intDataRot.featureAlignment = true;
    intDataRot.autoFeatureFunc = true;
    intDataRot.featureIndices = featureIndices;
    std::cout<<"Doing aligned rotational-seamless integration"<<std::endl;
    directional::setup_integration(rawField, intDataRot, meshCut, combedField);
    directional::integrate(combedField, intDataRot, meshCut, NFunctionRot, NCornerFunc);
    std::cout<<"Done with aligned rotational-seamless integration"<<std::endl;
    
    directional::IntegrationData intDataFull(N);
    intDataFull.verbose = true;
    intDataFull.integralSeamless = true;
    intDataFull.featureAlignment = true;
    intDataFull.autoFeatureFunc = true;
    intDataFull.featureIndices = featureIndices;
    std::cout<<"Doing aligned full-seamless integration"<<std::endl;
    directional::setup_integration(rawField, intDataFull, meshCut, combedField);
    directional::integrate(combedField, intDataFull, meshCut, NFunctionRot, NCornerFunc);
    std::cout<<"Done with aligned full-seamless integration"<<std::endl;
    
    viewer.init();
    viewer.set_surface_mesh(meshWhole);
    viewer.set_cartesian_field(rawField);
    std::vector<int> seamsList;
    for (int i=0;i<meshWhole.F.rows();i++)
        for (int j=0;j<3;j++)
            if (intDataFull.face2cut(i, j))
                seamsList.push_back(meshWhole.FE(i,j));
       
    Eigen::VectorXi seams = Eigen::VectorXi::Map(seamsList.data(), seamsList.size());
    viewer.highlight_edges(seams, "Seams");
    viewer.set_isolines(meshCut, NFunctionNoAlign);
    viewer.set_callback(callbackFunc);
    viewer.launch();
    
}


