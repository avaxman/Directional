#include <directional/readOFF.h>
#include <directional/streamlines.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>

directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField field, powerField;
directional::DirectionalViewer viewer;

bool advanceTracing = false;

int N = 3;
int anim_t = 0;
int anim_t_dir = 1;

void callbackFunc()
{
    ImGui::PushItemWidth(100);
    if (advanceTracing)
        viewer.advance_streamlines(0.5);
    if (ImGui::Button("stop tracing"))
        advanceTracing=false;
    else
        if (ImGui::Button("resume tracing"))
            advanceTracing=true;
    
    ImGui::PopItemWidth();
}


int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;
    
    directional::readOFF(TUTORIAL_DATA_PATH "/lion.off", mesh);
    viewer.init();
    viewer.set_callback(&callbackFunc);
    ftb.init(mesh);
    
    // Create a power field
    Eigen::VectorXi constFaces(1); constFaces(0) = 0;
    Eigen::MatrixXd constVectors(1, 3); constVectors.row(0) <<(mesh.V.row(mesh.F(0, 1)) - mesh.V.row(mesh.F(0, 0))).normalized();
    Eigen::VectorXd alignWeights(1); alignWeights(0) = -1.0;
    directional::power_field(ftb, constFaces, constVectors, alignWeights ,N, powerField);
    directional::power_to_raw(powerField,N,field, true);
    
    viewer.set_surface_mesh(mesh);
    viewer.set_cartesian_field(field);
    viewer.init_streamlines(field);
    viewer.advance_streamlines(0.5);  //to get the initial step
    viewer.launch();
}
