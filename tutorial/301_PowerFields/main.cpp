#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/readOBJ.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

directional::TriMesh mesh;
directional::IntrinsicVertexTangentBundle vtb;
directional::CartesianField rawField,powerFieldHard, powerFieldSoft;
Eigen::VectorXi constVertices;
Eigen::MatrixXd constVectors;
Eigen::VectorXd alignWeights;

directional::DirectionalViewer viewer;

int N = 4;
bool normalizeField = false;
bool viewFieldHard = true;

void recompute_field(){
    directional::power_field(vtb, constVertices, constVectors, Eigen::VectorXd::Constant(constVertices.size(),-1.0), N,powerFieldHard, normalizeField);
    directional::power_field(vtb, constVertices, constVectors, alignWeights, N,powerFieldSoft, normalizeField);
}

void update_visualization()
{
    directional::power_to_raw((viewFieldHard ? powerFieldHard : powerFieldSoft), N, rawField);
    
    directional::principal_matching(rawField);
    viewer.set_cartesian_field(rawField,"Power field");
    
    //Ghost mesh just showing field, to compare against constraints
    Eigen::MatrixXd constraintIntField = Eigen::MatrixXd::Zero(powerFieldHard.intField.rows(),2);
    for (int i=0;i<constVertices.size();i++)
        constraintIntField.row(constVertices(i))=powerFieldHard.intField.row(constVertices(i));
    
    directional::CartesianField constraintRawField, constraintPowerField;
    constraintPowerField.init(*(rawField.tb),directional::fieldTypeEnum::POWER_FIELD,N);
    constraintPowerField.set_intrinsic_field(constraintIntField);
    
    directional::power_to_raw(constraintPowerField, N, constraintRawField);
    viewer.set_cartesian_field(constraintRawField, "Constraints", 1);
    viewer.highlight_vertices(constVertices);
    
}

void callbackFunc()
{
    viewer.psSurfaceMeshList[0]->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);
    
    if (ImGui::Checkbox("Normalize Field", &normalizeField)){
        recompute_field();
        update_visualization();
    }
    
    ImGuiIO& io = ImGui::GetIO();
    if (io.MouseClicked[0]) { // if clicked
        glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
        polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);

        if(pickResult.isHit && pickResult.structure == viewer.psSurfaceMeshList[0]) {
            polyscope::SurfaceMeshPickResult meshPickResult =
            viewer.psSurfaceMeshList[0]->interpretPickResult(pickResult);
            
            //Polyscope doesn't allow different color vectors, so we are setting two fields - one combed one not
            if(meshPickResult.elementType == polyscope::MeshElement::FACE) {
                
                Eigen::Vector3d::Index maxCol;
                int currF = meshPickResult.index;
                Eigen::RowVector3d baryInFace; baryInFace<<meshPickResult.baryCoords[0], meshPickResult.baryCoords[1], meshPickResult.baryCoords[2];
                baryInFace.maxCoeff(&maxCol);
                int currVertex=mesh.F(currF, maxCol);
                int i;
                for (i = 0; i < constVertices.rows(); i++)
                    if (constVertices(i) == currVertex)
                        break;
                if (i == constVertices.rows())
                {
                    constVertices.conservativeResize(constVertices.rows() + 1);
                    constVertices(i) = currVertex;
                    constVectors.conservativeResize(constVectors.rows() + 1, 3);
                    alignWeights.conservativeResize(alignWeights.size()+1);
                    if (alignWeights.size()==1)
                        alignWeights(0)=1.0;
                    else
                        alignWeights(alignWeights.size()-1)=alignWeights(alignWeights.size()-2);
                }
                
                // Compute direction from the center of the face to the mouse
                constVectors.row(i) =(mesh.V.row(mesh.F(currF, 0)) * baryInFace(0) +
                                      mesh.V.row(mesh.F(currF, 1)) * baryInFace(1) +
                                      mesh.V.row(mesh.F(currF, 2)) * baryInFace(2) -
                                      mesh.V.row(currVertex)).normalized();
                
                recompute_field();
                update_visualization();
                
            }
        }
    }
}


int main()
{
    directional::readOBJ(TUTORIAL_DATA_PATH "/rocker-arm2500.obj", mesh);
    vtb.init(mesh);
    powerFieldHard.init(vtb, directional::fieldTypeEnum::POWER_FIELD, N);
    powerFieldSoft.init(vtb, directional::fieldTypeEnum::POWER_FIELD, N);
    constVertices.resize(0);
    constVectors.resize(0, 3);
    
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_callback(callbackFunc);
    recompute_field();
    update_visualization();
    viewer.launch();
}
