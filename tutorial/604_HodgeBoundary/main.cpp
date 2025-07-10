#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/discrete_exterior_calculus.h>
#include <directional/cochain_complex.h>
#include <directional/mass_matrices.h>

directional::TriMesh mesh;
Eigen::VectorXd z0, z1, z2, z1Exact, z1Coexact, z1Harmonic;
directional::DirectionalViewer viewer;


int main()
{
    bool dirichletBoundary = false;
    directional::readOFF(TUTORIAL_DATA_PATH "/FlatMesh.off",mesh);
    
    Eigen::SparseMatrix<double> d0 = directional::d0_matrix<double>(mesh, false);
    Eigen::SparseMatrix<double> d1 = directional::d1_matrix<double>(mesh, false);
    Eigen::SparseMatrix<double> h1, invh1, h2, invh2;
    directional::hodge_star_1_matrix(mesh, h1, invh1);
    directional::hodge_star_2_matrix<double>(mesh, h2, invh2);
    
    //demonstrating the exact sequences even with boundary conditions
    Eigen::SparseMatrix<double> d1d0 = d1*d0;
    double maxAbsValue = 0.0;
    for (int k = 0; k < d1d0.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(d1d0, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));
    
    std::cout<<"exact sequence identity (should be exactly zero): "<<maxAbsValue<<std::endl;
    
    Eigen::MatrixXd harmBasis;
    int bettiNumber = mesh.EV.rows() - (mesh.V.rows()-1) - mesh.F.rows();
    directional::cohomology_basis(d0, d1, h1,  bettiNumber, harmBasis);
    std::cout<<"Betti number: "<<bettiNumber<<std::endl;
    
    Eigen::RowVector3d COM = mesh.V.colwise().mean();
    Eigen::VectorXd z0GT(mesh.dcel.vertices.size()), curlGT(mesh.dcel.faces.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        z0GT[i] = 3.0*cos((mesh.V(i,0)-COM(1))/4.0)*cos((mesh.V(i,2)-COM(2))/4.0);
    
    //vertexVec.array()/=vertexVec.mean();
    
    for (int i=0;i<mesh.dcel.faces.size();i++){
        Eigen::RowVector3d midFacePointCOM = mesh.barycenters.row(i) - COM;
        curlGT[i] = mesh.faceAreas[i]*(sin(midFacePointCOM(0)/2.0)+sin(midFacePointCOM(1)/2.0)+sin(midFacePointCOM(2)/2.0));
    }
    curlGT.array()-=curlGT.mean();  //Due to Neumann boundary conditions (tangent coexact field), vector potential function adds up to zero
    
    //Generating an artificial composition and then reproducing it through decomposition.
    Eigen::VectorXd harmGT = harmBasis.col(0);
    Eigen::VectorXd z1ExactGT = d0*z0GT;
    Eigen::VectorXd z1CoexactGT, curlFiltered;
    directional::project_exact(d1, h1, curlGT, z1CoexactGT, curlFiltered, true);
    
    //creating balanced GT results for exact, coexact, and harmonic
    z1ExactGT.array()/=sqrt(((z1ExactGT.transpose()*h1*z1ExactGT).coeff(0,0)));
    z1CoexactGT.array()/=sqrt(((z1CoexactGT.transpose()*h1*z1CoexactGT).coeff(0,0)));
    
    z1 = harmGT + z1ExactGT + z1CoexactGT;
    
    Eigen::SparseMatrix<double> I(mesh.F.rows(), mesh.F.rows());
    I.setIdentity();
    directional::hodge_decomposition<double>(d0, d1, h1, h2, z1, 1,  z1Exact, z1Coexact, z1Harmonic, z0,  z2);
    std::cout<<"Exact reproduction: "<<(z1Exact - z1ExactGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Coexact reproduction: "<<(z1Coexact - z1CoexactGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Harmonic reproduction: "<<(z1Harmonic - harmGT).cwiseAbs().maxCoeff()<<std::endl;
    
    //triangle mesh setup
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_surface_vertex_data(z0, "0-form potential", 0);
    viewer.set_1form(z1,"Original field", 0, 0, 5.0);
    viewer.set_1form(z1Exact,"Exact field", 0, 1, 5.0);
    Eigen::VectorXd dualz0 = z2.array()/mesh.faceAreas.array();
    viewer.set_surface_face_data(dualz0, "dual 0-form potential", 0);
    viewer.set_1form(z1Coexact,"Coexact field", 0, 2, 5.0);
    viewer.set_1form(z1Harmonic,"Harmonic field", 0, 3, 5.0);
    viewer.launch();
}
