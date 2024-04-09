#include <vector>
#include <cstdlib>
#include <directional/readDMAT.h>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/conjugate_frame_fields.h>
#include <directional/ConjugateFFSolverData.h>
#include <directional/directional_viewer.h>

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawFieldOrig, rawFieldConjugate;
directional::DirectionalViewer viewer;
Eigen::VectorXd conjugacyOrig, conjugacyConjugate;

int N=4;
double conjMaxOrig;

// Input constraints
Eigen::VectorXi b;
Eigen::MatrixXd bc;

typedef enum {ORIGINAL_FIELD, ORIGINAL_CONJUGACY, OPTIMIZED_FIELD, OPTIMIZED_CONJUGACY} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;


void update_visualization()
{
    if ((viewingMode==ORIGINAL_FIELD)||(viewingMode==OPTIMIZED_FIELD)){
        viewer.highlight_faces(b);
    }else{
        Eigen::VectorXd currConjugacy = (viewingMode==ORIGINAL_CONJUGACY ? conjugacyOrig: conjugacyConjugate);
        viewer.set_face_data(currConjugacy,0.0,conjMaxOrig);
        viewer.highlight_faces(Eigen::VectorXi());
    }

}


void callbackFunc() {
    ImGui::PushItemWidth(100);


    const char* items[] = { "Original field", "Original conjugacy", "Optimized field", "Optimized conjugacy"};
    static const char* current_item = NULL;

    if (ImGui::BeginCombo("##combo", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)){
                switch (n){
                    case 0:
                        viewingMode = ORIGINAL_FIELD;
                        break;
                    case 1:
                        viewingMode = ORIGINAL_CONJUGACY;
                        break;
                    case 2:
                        viewingMode = OPTIMIZED_FIELD;
                        break;
                    case 3:
                        viewingMode = OPTIMIZED_CONJUGACY;
                        break;
                }
                update_visualization();
            }
            current_item = items[n];
            if (is_selected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

    if (ImGui::Button("Save optimized Field"))
        directional::write_raw_field(TUTORIAL_DATA_PATH "/conjugacy.rawfield", combedFieldCF);

    ImGui::PopItemWidth();
}



int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;


    // Load a mesh in OBJ format
    directional::readOBJ(TUTORIAL_DATA_PATH "/inspired_mesh.obj", mesh);
    ftb.init(mesh);

    // Load constraints
    Eigen::VectorXi bFull;
    directional::readDMAT(TUTORIAL_DATA_PATH "/inspired_mesh_b.dmat",bFull);
    Eigen::MatrixXd bcFull;
    directional::readDMAT(TUTORIAL_DATA_PATH "/inspired_mesh_bc.dmat",bcFull);

    bcFull.conservativeResize(bcFull.rows(), 3*N);
    bcFull.block(0,6,bcFull.rows(),6) = -bcFull.block(0,0,bcFull.rows(),6);

    //putting the constraints in single-vector matrices
    bc.resize(N*bcFull.rows(),3);
    b.resize(N*bFull.rows());
    for (int i=0;i<bFull.rows();i++){
        for (int j=0;j<N;j++){
            b(i*N+j)=bFull(i);
            bc.row(i*N+j)=bcFull.block(i,3*j,1,3);
        }
    }

    //initial solution
    directional::CartesianField pvField;
    directional::polyvector_field(ftb, b, bc, N, pvField);
    directional::polyvector_to_raw(pvField, rawFieldOrig);

    directional::principal_matching(rawFieldOrig);

    // Initialize conjugate field with smooth field
    directional::ConjugateFFSolverData csdata(mesh);
    directional::conjugate_frame_fields(csdata, b, rawFieldOrig, rawFieldConjugate);
    directional::principal_matching(rawFieldConjugate);
    csdata.evaluateConjugacy(rawFieldOrig, conjugacyOrig);
    conjMaxOrig = conjugacyOrig.lpNorm<Infinity>();
    csdata.evaluateConjugacy(rawFieldConjugate, conjugacyConjugate);

    viewer.init();
    viewer.set_mesh(mesh);
    viewer.set_field(rawFieldOrig, "", 0,0);
    viewer.set_field(rawFieldConjugate, "", 0, 1);
    update_visualization();

    viewer.set_callback(callbackFunc);
    viewer.launch();
}
