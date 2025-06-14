#include <directional/CartesianField.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/TriMesh.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/principal_matching.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/readOBJ.h>

int N=2;
directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField field, powerField;
int sparsity=0;

directional::DirectionalViewer viewer;

void callbackFunc() {
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
    ImGui::Text("Sparsity: %d", sparsity);
    ImGui::SameLine();
    if (ImGui::Button("+")) {
        sparsity++;
        viewer.set_cartesian_field(field,"",0, sparsity);
        viewer.set_glyph_length(0.3*((double)sparsity+1.0));
    }
    ImGui::SameLine();
    if (ImGui::Button("-")) {
        sparsity--;
        viewer.set_cartesian_field(field,"",0,sparsity);
        viewer.set_glyph_length(0.3*((double)sparsity+1.0));
    }

    ImGui::PopItemWidth();
}


int main()
{

  directional::readOBJ(TUTORIAL_DATA_PATH "/armadillo.obj",mesh);
  ftb.init(mesh);
  viewer.init();
  viewer.set_callback(&callbackFunc);
  Eigen::VectorXi bc(1); bc(0)=0;
  Eigen::MatrixXd b(1,3); b.row(0)=mesh.V.row(mesh.F(0,1))-mesh.V.row(mesh.F(0,2));
  directional::power_field(ftb, bc,b, Eigen::VectorXd::Constant(bc.size(),-1), N, powerField);
  directional::power_to_raw(powerField,N,field, true);
  
  viewer.set_surface_mesh(mesh);
  viewer.set_cartesian_field(field);
  viewer.launch();
}

