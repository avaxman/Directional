#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOBJ.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/directional_viewer.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/pick.h>

int currF=0, N;
directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField field;
directional::DirectionalViewer viewer;


/*void update_raw_field_mesh()
 {
 //configuring just part of the faces to have combed coloring
 Eigen::Vector3i otherFaces;
 Eigen::Vector3i zeroInFace;
 for (int i=0;i<3;i++){
 otherFaces(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? mesh.EF(mesh.FE(currF,i),1) : mesh.EF(mesh.FE(currF,i),0));
 zeroInFace(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? field.matching(mesh.FE(currF,i)) : -field.matching(mesh.FE(currF,i)));
 }
 
 Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(mesh.F.rows(),N);
 glyphColors.row(currF)=directional::DirectionalViewer::indexed_glyph_colors(field.extField.row(currF), false);
 for (int i=0;i<N;i++)
 for (int j=0;j<3;j++)
 glyphColors.block(otherFaces(j),3*((i+zeroInFace(j)+N)%N),1,3)<<glyphColors.row(currF).segment(3*i,3);
 
 viewer.set_field(field, glyphColors);
 
 }*/

void callbackFunc()
{
    viewer.psSurfaceMeshList[0]->setSelectionMode(polyscope::MeshSelectionMode::FacesOnly);
    
    // get the mouse location from ImGui
    ImGuiIO& io = ImGui::GetIO();
    if (io.MouseClicked[0]) { // if clicked
        glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
        polyscope::PickResult pickResult = polyscope::pickAtScreenCoords(screenCoords);
        
        // check out pickResult.isHit, pickResult.structureName, pickResult.depth, etc
        
        // get additional information if we clicked on a mesh
        if(pickResult.isHit && pickResult.structure == viewer.psSurfaceMeshList[0]) {
            polyscope::SurfaceMeshPickResult meshPickResult =
            viewer.psSurfaceMeshList[0]->interpretPickResult(pickResult);
            
            //Polyscope doesn't allow different color vectors, so we are setting two fields - one combed one not
            if(meshPickResult.elementType == polyscope::MeshElement::FACE) {
                std::cout << "clicked face " << meshPickResult.index << std::endl;
                int currF = meshPickResult.index;
                Eigen::MatrixXd rawFieldCombed(4,3*N);
                Eigen::MatrixXd sourcesCombed(4,3);
                Eigen::VectorXi combingMask = Eigen::VectorXi::Zero(field.extField.rows());
                Eigen::Vector3i otherFaces;
                Eigen::Vector3i zeroInFace;
                sourcesCombed.row(0) = mesh.barycenters.row(currF)+0.001*mesh.faceNormals.row(currF);
                rawFieldCombed.row(0) = field.extField.row(currF);
                for (int i=0;i<3;i++){
                    otherFaces(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? mesh.EF(mesh.FE(currF,i),1) : mesh.EF(mesh.FE(currF,i),0));
                    zeroInFace(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? field.matching(mesh.FE(currF,i)) : -field.matching(mesh.FE(currF,i)));
                    zeroInFace(i) = (zeroInFace(i)+field.N)%N;
                    //reordering the vectors in the rawfield
                    sourcesCombed.row(i+1)<<mesh.barycenters.row(otherFaces(i))+0.001*mesh.faceNormals.row(otherFaces(i));  //to avoid z-fighting. A bit hackish!
                    rawFieldCombed.row(i+1)<<field.extField.row(otherFaces(i)).segment(3*zeroInFace(i),3*field.N-3*zeroInFace(i)),field.extField.row(otherFaces(i)).segment(0,3*zeroInFace(i));
                    
                }
                viewer.set_raw_field(sourcesCombed, rawFieldCombed, mesh.avgEdgeLength, "Combed Field",  1);
                viewer.toggle_combed_colors(true, false, 1);
            }
        }
    }
}

int main()
{
    directional::readOBJ(TUTORIAL_DATA_PATH "/lilium.obj", mesh);
    ftb.init(mesh);
    viewer.init();
    directional::read_raw_field(TUTORIAL_DATA_PATH "/lilium.rawfield", ftb, N, field);
    directional::principal_matching(field);
    
    //triangle mesh setup
    viewer.set_surface_mesh(mesh);
    viewer.set_cartesian_field(field);
    viewer.set_callback(callbackFunc);
    viewer.launch();
}

