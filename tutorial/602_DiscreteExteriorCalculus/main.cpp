#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/dec.h>
#include <directional/poly_to_constant_subdivision.h>
#include <directional/CochainComplex.h>
#include <directional/mass_matrices.h>


directional::TriMesh mesh;
Eigen::VectorXd z1, z1Diag, z2Exact;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/Eight.off",mesh);

    Eigen::VectorXd z2(mesh.dcel.faces.size());
    /*for (int i=0;i<mesh.dcel.vertices.size();i++)
        z0[i] = sin(mesh.V(i,0)/20.0)*cos(mesh.V(i,1)/20.0);*/

    //giving a 2-form of constant (pointwise curl).
    for (int i=0;i<mesh.dcel.faces.size();i++){
        //Eigen::RowVector3d midFacePoint = (mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1)));
        z2[i] = mesh.faceAreas(i)*sin(mesh.barycenters(i,0)/20.0)*cos(mesh.barycenters(i,1)/20.0); //sin(midEdgePoint(1)/20.0)*cos(midEdgePoint(2)/20.0);
    }
    //making it physical curl
    z2.array()-=z2.mean();
    //z2 = mesh.faceAreas;

    Eigen::SparseMatrix<double> d0 = directional::d0_matrix<double>(&mesh);
    Eigen::SparseMatrix<double> d1 = directional::d1_matrix<double>(&mesh);
    Eigen::SparseMatrix<double> hodgeStar, invHodgeStar;
    directional::hodge_star_1_matrix(&mesh, hodgeStar, invHodgeStar, true);
    Eigen::SparseMatrix<double> M2 = directional::face_scalar_mass_matrix_2D<double>(&mesh);
    Eigen::SparseMatrix<double> M1;
    directional::linear_whitney_mass_matrix(&mesh, M1);

    Eigen::SparseMatrix<double> L1 = d0.adjoint()*M1*d0;
    Eigen::SparseMatrix<double> L2 = d0.adjoint()*hodgeStar*d0;
    Eigen::SparseMatrix<double> diff = M1-hodgeStar;
    //Eigen::VectorXd linFunc = mesh.V.col(0);
    //Eigen::VectorXd lapLin = L1*linFunc;

    //std::cout<<"lapLin.cwiseAbs().maxCoeff(): "<<lapLin.cwiseAbs().maxCoeff()<<std::endl;

    double maxAbsValue = 0.0;
    for (int k = 0; k < diff.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(diff, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));

    std::cout<<"exact laplacian identity (should be exactly zero): "<<maxAbsValue<<std::endl;
                                                                      //z1Exact = d0*z0;
    directional::project_exact(d1, hodgeStar, z2, z1Diag, z2Exact, true);

    std::cout<<"Reproducing the original curl (z2Diag): "<<(z2-z2Exact).cwiseAbs().maxCoeff()<<std::endl;

    directional::project_exact(d1, M1, z2, z1, z2Exact, true);

    std::cout<<"Reproducing the original curl (z2): "<<(z2-z2Exact).cwiseAbs().maxCoeff()<<std::endl;

    std::cout<<"Difference between M1 and hodge star : "<<(z1-z1Diag).cwiseAbs().maxCoeff()<<std::endl;

    //triangle mesh setup
    viewer.init();
    viewer.set_mesh(mesh);
    /*viewer.set_vertex_data(z0, z0.minCoeff(), z0.maxCoeff(),"Primal 0-form", 0);
    viewer.set_1form(&mesh, z1Exact,"Exact field", 0, 0, 2, 0.2, 2.0);*/
    viewer.set_face_data(z2, z2.minCoeff(), z2.maxCoeff(),"dual 0-form", 0);
    Eigen::MatrixXd formField = viewer.set_1form(&mesh, z1,"Coexact field", 0, 0, 2, 0.3, 10.0);
    Eigen::MatrixXd formFieldDiag = viewer.set_1form(&mesh, z1Diag,"Coexact field diag", 0, 1, 2, 0.3, 10.0);

    viewer.launch();
}
