#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/CochainComplex.h>
#include <directional/gradient_matrices.h>
#include <directional/curl_matrices.h>
#include <directional/mass_matrices.h>
#include <directional/extrinsic_intrinsic_matrices.h>

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField harmField;
directional::DirectionalViewer viewer;
Eigen::MatrixXd harmBasis;
Eigen::SparseMatrix<double> IE;
int currHarmBasis = 0;

void callbackFunc() {
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
    ImGui::Text("Harmonic basis function: %d", currHarmBasis);
    ImGui::SameLine();
    if (ImGui::Button("+")) {
        currHarmBasis = (currHarmBasis + 1)%harmBasis.cols();
        //std::cout<<"currHarmBasis: "<<currHarmBasis<<std::endl;
        harmField.set_extrinsic_field(IE*harmBasis.col(currHarmBasis));
        viewer.set_field(harmField);
        //std::cout<<"after set field"<<std::endl;
    }
    ImGui::SameLine();
    if (ImGui::Button("-")) {
        currHarmBasis = (currHarmBasis + harmBasis.cols() - 1)%harmBasis.cols();
        //std::cout<<"currHarmBasis: "<<currHarmBasis<<std::endl;
        harmField.set_extrinsic_field(IE*harmBasis.col(currHarmBasis));
        viewer.set_field(harmField);
    }

    ImGui::PopItemWidth();
}


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/1146164.off",mesh);
    ftb.init(mesh);

    //TODO: create the actual field

    //Must use intrinsic since otherwise the harmonic field will have spurious normal components
    Eigen::SparseMatrix<double> G = directional::conf_gradient_matrix_2D<double>(&mesh, true);
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(&mesh, Eigen::VectorXi(), true);
    Eigen::SparseMatrix<double> Mx = directional::face_vectors_mass_matrix_2D<double>(&mesh, true);
    Eigen::SparseMatrix<double> iMx = directional::face_vectors_mass_matrix_2D<double>(&mesh, true, true);
    Eigen::SparseMatrix<double> Mc = directional::edge_diamond_mass_matrix_2D<double>(&mesh, true);
    IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(&mesh);

    int bettiNumber;
    directional::cohomology_basis(G, C, Mx,  harmBasis, bettiNumber);

    std::cout<<"divergence of harmonic basis: "<<(G.adjoint()*Mx*harmBasis).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"curl of harmonic basis: "<<(C*harmBasis).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"bettiNumber: "<<bettiNumber<<std::endl;

    viewer.init();
    viewer.set_callback(&callbackFunc);
    viewer.set_mesh(mesh);
    //std::cout<<"before set extrinsic field "<<std::endl;
    harmField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    harmField.set_extrinsic_field(IE*harmBasis.col(currHarmBasis));
    //std::cout<<"after set extrinsic field "<<std::endl;
    viewer.set_field(harmField);
    viewer.launch();
}
