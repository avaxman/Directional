#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/dec.h>
#include <directional/CochainComplex.h>
#include <directional/mass_matrices.h>

directional::TriMesh mesh;
Eigen::VectorXd z0, z1, z2, z1Exact, z1Coexact, z1Harmonic;
directional::DirectionalViewer viewer;


int main()
{
    bool dirichletBoundary = false;
    directional::readOFF(TUTORIAL_DATA_PATH "/FlatMesh.off",mesh);
    //ftb.init(mesh);

    Eigen::SparseMatrix<double> d0 = directional::d0_matrix<double>(mesh, false);
    Eigen::SparseMatrix<double> d1 = directional::d1_matrix<double>(mesh, false);
    Eigen::SparseMatrix<double> hodgeStar, invHodgeStar;
    directional::hodge_star_1_matrix(mesh, hodgeStar, invHodgeStar, false);
    Eigen::SparseMatrix<double> M2 = directional::face_scalar_mass_matrix_2D<double>(mesh);

    //demonstrating the exact sequences
    Eigen::SparseMatrix<double> d1d0 = d1*d0;
    double maxAbsValue = 0.0;
    for (int k = 0; k < d1d0.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(d1d0, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));


    std::cout<<"exact sequence identity (should be exactly zero): "<<maxAbsValue<<std::endl;

    Eigen::MatrixXd harmBasis;
    int bettiNumber;
    directional::cohomology_basis(d0, d1, hodgeStar,  harmBasis, bettiNumber);
    std::cout<<"Euler characteristic: "<<mesh.V.rows()-mesh.EV.rows()+mesh.F.rows()<<std::endl;
    std::cout<<"Betti number: "<<bettiNumber<<std::endl;
    std::cout<<"divergence of harmonic basis: "<<(d0.adjoint()*hodgeStar*harmBasis).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"curl of harmonic basis: "<<(d1*harmBasis).cwiseAbs().maxCoeff()<<std::endl;

    Eigen::VectorXd harmVec = harmBasis.col(0);
    Eigen::RowVector3d COM = mesh.V.colwise().mean();
    //std::cout<<"COM: "<<COM<<std::endl;
    Eigen::VectorXd vertexVec(mesh.dcel.vertices.size()), midFaceVec(mesh.dcel.faces.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        vertexVec[i] = 3.0*cos((mesh.V(i,0)-COM(1))/4.0)*cos((mesh.V(i,2)-COM(2))/4.0);

    vertexVec.array()/=vertexVec.mean();

    for (int i=0;i<mesh.dcel.faces.size();i++){
        Eigen::RowVector3d midFacePointCOM = mesh.barycenters.row(i) - COM;
        midFaceVec[i] = mesh.faceAreas[i]*(sin(midFacePointCOM(0)/2.0)+sin(midFacePointCOM(1)/2.0)+sin(midFacePointCOM(2)/2.0));
    }
    midFaceVec.array()/=midFaceVec.mean();

    Eigen::VectorXd z1CoexactPres, z2Filtered;
    directional::project_exact(d1, hodgeStar, midFaceVec, z1CoexactPres, z2Filtered, true);
    
    //creating balanced GT results for exact, coexact, and harmonic
    Eigen::VectorXd z1ExactGT = d0*vertexVec;
    z1ExactGT.array()/=sqrt(((z1ExactGT.transpose()*hodgeStar*z1ExactGT).coeff(0,0)));
    Eigen::VectorXd z1CoexactGT = z1CoexactPres.array()/sqrt(((z1CoexactPres.transpose()*hodgeStar*z1CoexactPres).coeff(0,0)));
    

    z1 = harmVec +z1ExactGT +  z1CoexactGT;  //harmvec is already normalized
   
    std::cout<<"before hodge decomposition"<<std::endl;
    directional::hodge_decomposition<double>(d0, d1, hodgeStar, M2, z1, z0, z1Exact, z1Coexact, z2, z1Harmonic);
    std::cout<<"after hodge decomposition"<<std::endl;
    std::cout<<"Exact reproduction: "<<(z1Exact - z1ExactGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Coexact reproduction: "<<(z1Coexact - z1CoexactGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Harmonic reproduction: "<<(z1Harmonic - harmVec).cwiseAbs().maxCoeff()<<std::endl;

    //triangle mesh setup
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_vertex_data(z0, z0.minCoeff(), z0.maxCoeff(),"0-form", 0);
    viewer.set_1form(z1,"Original field", 0, 0, 2, 0.2, 5.0);
    viewer.set_1form(z1Exact,"Exact field", 0, 1, 2, 0.2, 5.0);
    viewer.set_face_data(midFaceVec, midFaceVec.minCoeff(), midFaceVec.maxCoeff(),"2-form", 0);
    viewer.set_1form(z1Coexact,"Coexact field", 0, 2, 2, 0.2, 5.0);
    viewer.set_1form(z1Harmonic,"Harmonic field", 0, 3, 2, 0.2, 5.0);
    viewer.launch();
}
