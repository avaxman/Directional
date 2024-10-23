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
Eigen::VectorXd z1Exact, z1Coexact, z1CoexactDiag;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/Eight.off",mesh);

    Eigen::VectorXd z0(mesh.dcel.vertices.size()), z2(mesh.dcel.faces.size());
    /*for (int i=0;i<mesh.dcel.vertices.size();i++)
        z0[i] = sin(mesh.V(i,0)/20.0)*cos(mesh.V(i,1)/20.0);*/
    for (int i=0;i<mesh.dcel.faces.size();i++){
        Eigen::RowVector3d midFacePoint = (mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1)));
        z2[i] = mesh.faceAreas(i)*sin(mesh.barycenters(i,0)/20.0)*cos(mesh.barycenters(i,1)/20.0); //sin(midEdgePoint(1)/20.0)*cos(midEdgePoint(2)/20.0);
    }

    Eigen::SparseMatrix<double> d0 = directional::d0_matrix<double>(&mesh);
    Eigen::SparseMatrix<double> d1 = directional::d1_matrix<double>(&mesh);
    Eigen::SparseMatrix<double> hodgeStar, invHodgeStar;
    directional::hodge_star_1_matrix(&mesh, hodgeStar, invHodgeStar, true);
    Eigen::SparseMatrix<double> M2 = directional::face_scalar_mass_matrix_2D<double>(&mesh);
    Eigen::SparseMatrix<double> M1;
    directional::linear_whitney_mass_matrix(&mesh, M1);

    //checking M1:

    Eigen::SparseMatrix<double> L1 = d0.adjoint()*M1*d0;
    Eigen::SparseMatrix<double> L2 = d0.adjoint()*hodgeStar*d0;
    Eigen::SparseMatrix<double> diff = L1-L2;
    Eigen::VectorXd linFunc = mesh.V.col(0);
    Eigen::VectorXd lapLin = L1*linFunc;

    std::cout<<"lapLin.cwiseAbs().maxCoeff(): "<<lapLin.cwiseAbs().maxCoeff()<<std::endl;

    double maxAbsValue = 0.0;
    for (int k = 0; k < diff.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(diff, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));

    std::cout<<"exact laplacian identity (should be exactly zero): "<<maxAbsValue<<std::endl;
                                                                      //z1Exact = d0*z0;
    z1CoexactDiag = directional::project_coexact(d1, hodgeStar, M2, z2,
    Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
    Eigen::Vector<NumberType, Eigen::Dynamic>& nextCochain
            invHodgeStar * d1.transpose() * z2dual;   //There is a minus issue.
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(M1);
    if(solver.info()!=Eigen::Success) {
        // decomposition failed
        std::cout<<"factorization failed!!!"<<std::endl;
        exit(0);
    }
    z1Coexact = solver.solve(-d1.transpose() * z2dual);
    if(solver.info()!=Eigen::Success) {
        // solving failed
        std::cout<<"solving failed!!!"<<std::endl;
        exit(0);
    }
    //demonstrating the exact sequences
    /*Eigen::SparseMatrix<double> d1d0 = d1*d0;
    double maxAbsValue = 0.0;
    for (int k = 0; k < d1d0.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(d1d0, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));*/


   // std::cout<<"exact sequence identity (should be exactly zero): "<<maxAbsValue<<std::endl;

    //Eigen::VectorXd curl = d1*z1Coexact;

    //triangle mesh setup
    viewer.init();
    viewer.set_mesh(mesh);
    /*viewer.set_vertex_data(z0, z0.minCoeff(), z0.maxCoeff(),"Primal 0-form", 0);
    viewer.set_1form(&mesh, z1Exact,"Exact field", 0, 0, 2, 0.2, 2.0);*/
    viewer.set_face_data(z2dual, z2dual.minCoeff(), z2dual.maxCoeff(),"dual 0-form", 0);
    Eigen::MatrixXd formField = viewer.set_1form(&mesh, z1Coexact,"Coexact field", 0, 0, 2, 0.2, 20.0);
    Eigen::MatrixXd formFieldDiag = viewer.set_1form(&mesh, z1CoexactDiag,"Coexact field diag", 0, 1, 2, 0.2, 20.0);

    directional::TriMesh subdMesh;
    Eigen::VectorXi FFK;
    Eigen::MatrixXd rawFieldK;
    directional::poly_to_constant_subdivision(mesh, formField, 1, 1, subdMesh, FFK, rawFieldK);

    //std::cout<<"rawFieldK : "<<rawFieldK<<std::endl;
    directional::IntrinsicFaceTangentBundle subdftb;
    subdftb.init(subdMesh);
    directional::CartesianField rawFieldClass;
    rawFieldClass.init(subdftb, directional::fieldTypeEnum::RAW_FIELD,1);
    rawFieldClass.set_extrinsic_field(rawFieldK);

    directional::poly_to_constant_subdivision(mesh, formFieldDiag, 1, 1, subdMesh, FFK, rawFieldK);

    //std::cout<<"rawFieldK : "<<rawFieldK<<std::endl;
    directional::CartesianField rawFieldClassDiag;
    rawFieldClassDiag.init(subdftb, directional::fieldTypeEnum::RAW_FIELD,1);
    rawFieldClassDiag.set_extrinsic_field(rawFieldK);


    viewer.set_mesh(subdMesh, 1);
    //viewer.set_vertex_data(lapLin, lapLin.minCoeff(), lapLin.maxCoeff(), "lapLin",0);
    viewer.set_field(rawFieldClass, "subd-field", 1, 2);
    viewer.set_field(rawFieldClassDiag, "subd-field-diag", 1, 3);
    /*viewer.init_streamlines(1, Eigen::VectorXi(), 0.3);
    for (int i=0;i<100;i++)
        viewer.advance_streamlines(0.5, 1);  //to get the initial step*/

    viewer.launch();
}
