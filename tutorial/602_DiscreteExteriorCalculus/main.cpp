#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/dec.h>

directional::TriMesh mesh;
Eigen::VectorXd z1Exact, z1Coexact;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/Pegasus.off",mesh);

    Eigen::VectorXd z0(mesh.dcel.vertices.size()), z2dual(mesh.dcel.faces.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        z0[i] = sin(mesh.V(i,0)/20.0)*cos(mesh.V(i,1)/20.0);
    for (int i=0;i<mesh.dcel.faces.size();i++){
        Eigen::RowVector3d midFacePoint = (mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1)));
        z2dual[i] = sin(mesh.barycenters(i,0)/20.0)*cos(mesh.barycenters(i,1)/20.0); //sin(midEdgePoint(1)/20.0)*cos(midEdgePoint(2)/20.0);
    }

    Eigen::SparseMatrix<double> d0 = directional::d0_matrix<double>(&mesh);
    Eigen::SparseMatrix<double> d1 = directional::d1_matrix<double>(&mesh);
    Eigen::SparseMatrix<double> hodgeStar, invHodgeStar;
    directional::hodge_star_1_matrix(&mesh, hodgeStar, invHodgeStar, false);

    z1Exact = d0*z0;
    z1Coexact = invHodgeStar * d1.transpose() * z2dual;   //There is a minus issue.

    //demonstrating the exact sequences
    Eigen::SparseMatrix<double> d1d0 = d1*d0;
    double maxAbsValue = 0.0;
    for (int k = 0; k < d1d0.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(d1d0, k); it; ++it)
            maxAbsValue = std::max(maxAbsValue, std::abs(it.value()));


    std::cout<<"exact sequence identity (should be exactly zero): "<<maxAbsValue<<std::endl;

    Eigen::VectorXd curl = d1*z1Coexact;

    //triangle mesh setup
    viewer.init();
    viewer.set_mesh(mesh);
    viewer.set_vertex_data(z0, z0.minCoeff(), z0.maxCoeff(),"Primal 0-form", 0);
    viewer.set_1form(&mesh, z1Exact,"Exact field", 0, 0, 2, 0.2, 2.0);
    viewer.set_face_data(curl, curl.minCoeff(), curl.maxCoeff(),"dual 0-form", 0);
    viewer.set_1form(&mesh, z1Coexact,"Coexact field", 0, 1, 2, 0.2, 20.0);
    viewer.launch();
}
