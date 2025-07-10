#include <iostream>
#include <Eigen/Core>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/discrete_exterior_calculus.h>
//#include <directional/poly_to_constant_subdivision.h>
#include <directional/cochain_complex.h>
#include <directional/mass_matrices.h>


directional::TriMesh mesh;
Eigen::VectorXd z1, z1Diag, z2Exact;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOBJ(TUTORIAL_DATA_PATH "/697224__sf.obj",mesh);
    
    Eigen::VectorXd z2(mesh.dcel.faces.size());
    //giving a 2-form of constant (pointwise curl).
    for (int i=0;i<mesh.dcel.faces.size();i++)
        z2[i] = mesh.faceAreas(i)*sin(mesh.barycenters(i,0)/40.0)*cos(mesh.barycenters(i,1)/40.0);
    
    //making it a physical curl quantity
    z2.array()-=z2.mean();
    
    Eigen::SparseMatrix<double> d0 = directional::d0_matrix<double>(mesh);
    Eigen::SparseMatrix<double> d1 = directional::d1_matrix<double>(mesh);
    Eigen::SparseMatrix<double> hodgeStar, invHodgeStar;
    directional::hodge_star_1_matrix(mesh, hodgeStar, invHodgeStar);
    //Eigen::SparseMatrix<double> M2 = directional::face_mass_matrix_2D<double>(mesh, true);  //of inverse face areas, since -forms are integrated quantities
    Eigen::SparseMatrix<double> M1;
    directional::linear_whitney_mass_matrix(mesh, M1);
    
    Eigen::SparseMatrix<double> L1 = d0.adjoint()*M1*d0;
    Eigen::SparseMatrix<double> L2 = d0.adjoint()*hodgeStar*d0;
    Eigen::SparseMatrix<double> diff = L1-L2;
    
    double maxAbsValue = 0.0;
    for (int k = 0; k < diff.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(diff, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));
    
    std::cout<<"Exact laplacian identity (should be numerically zero): "<<maxAbsValue<<std::endl;
    directional::project_exact(d1, hodgeStar, z2, z1Diag, z2Exact, true);
    std::cout<<"Reproducing the original curl (z2Diag, should be numerically zero): "<<(z2-z2Exact).cwiseAbs().maxCoeff()<<std::endl;
    directional::project_exact(d1, M1, z2, z1, z2Exact, true);
    std::cout<<"Reproducing the original curl (z2, should be numerically zero): "<<(z2-z2Exact).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Difference between z1 and z1Diag (small, not zero): "<<(z1-z1Diag).cwiseAbs().maxCoeff()<<std::endl;
    
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_surface_face_data(z2, "Integrated face curl", 0);
    Eigen::MatrixXd formField = viewer.set_1form(z1,"Coexact field", 0, 0, 1.0, 2, 0.2);
    Eigen::MatrixXd formFieldDiag = viewer.set_1form(z1Diag,"Coexact field diag", 0, 1, 1.0, 2, 0.2);
    
    viewer.launch();
}
