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
#include <directional/PCFaceTangentBundle.h>

directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
Eigen::VectorXd gradFieldVec, rotCogradFieldVec, harmFieldVec, vertexFunction, diamondForm;
directional::CartesianField origField, gradField, rotCogradField, harmField;
directional::DirectionalViewer viewer;


int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/1146164.off",mesh);
    //directional::readOBJ(TUTORIAL_DATA_PATH "/torus.obj",mesh);
    
    ftb.init(mesh);

    //TODO: create the actual field
    //Must use intrinsic since otherwise the harmonic field will have spurious normal components
    Eigen::SparseMatrix<double> G = directional::conf_gradient_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(mesh, Eigen::VectorXi(), true);
    Eigen::SparseMatrix<double> Mx = directional::face_vectors_mass_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> iMx = directional::face_vectors_mass_matrix_2D<double>(mesh, true, true);
    Eigen::SparseMatrix<double> Mc = directional::edge_diamond_mass_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(mesh);

    Eigen::MatrixXd harmBasis;
    int bettiNumber;
    std::cout<<"Extracting the basis of harmonic fields..."<<std::endl;
    directional::cohomology_basis(G, C, Mx,  harmBasis, bettiNumber);
    //std::cout<<"Euler characteristic: "<<mesh.V.rows()-mesh.EV.rows()+mesh.F.rows()<<std::endl;
    //std::cout<<"Betti number: "<<bettiNumber<<std::endl;
    //std::cout<<"harmBasis: "<<harmBasis.block(0,0,20,harmBasis.cols())<<std::endl;

    std::cout<<"divergence of harmonic basis (numerically zero): "<<(G.adjoint()*Mx*harmBasis).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"curl of harmonic basis (numerically zero): "<<(C*harmBasis).cwiseAbs().maxCoeff()<<std::endl;

    Eigen::VectorXd harmVec = 10.0*harmBasis.col(0);
    Eigen::RowVector3d COM = mesh.V.colwise().mean();
    Eigen::VectorXd vertexVec(mesh.dcel.vertices.size()), midEdgeVec(mesh.dcel.edges.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        vertexVec[i] = 3.0*cos((mesh.V(i,0)-COM(1))/4.0)*cos((mesh.V(i,2)-COM(2))/4.0);

    vertexVec.array()/=vertexVec.mean();

    for (int i=0;i<mesh.dcel.edges.size();i++){
        Eigen::RowVector3d midEdgePointCOM = 0.5*(mesh.V.row(mesh.EV(i,0))+mesh.V.row(mesh.EV(i,1))) - COM;
        midEdgeVec[i] = sin(midEdgePointCOM(0)/4.0)+sin(midEdgePointCOM(1)/4.0)+sin(midEdgePointCOM(2)/4.0);
    }
    midEdgeVec.array()/=midEdgeVec.mean();

    origField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    Eigen::VectorXd gradientGT = G*vertexVec;
    gradientGT.array()/=sqrt((gradientGT.transpose()*Mx*gradientGT).coeff(0,0));
    Eigen::VectorXd rotCogradientGT = iMx*C.adjoint()*midEdgeVec;
    rotCogradientGT.array()/=sqrt((rotCogradientGT.transpose()*Mx*rotCogradientGT).coeff(0,0));
    origField.set_extrinsic_field(IE*(harmVec + gradientGT +  rotCogradientGT));

    std::cout<<"before hodge decomposition"<<std::endl;
    directional::hodge_decomposition<double>(G, C, Mx, Mc, origField.flatten(true), vertexFunction, gradFieldVec, rotCogradFieldVec, diamondForm, harmFieldVec);
    std::cout<<"after hodge decomposition"<<std::endl;
    std::cout<<"Exact reproduction (numerically zero): "<<(gradFieldVec -gradientGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Coexact reproduction (numerically zero): "<<(rotCogradFieldVec - rotCogradientGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Harmonic reproduction (numerically zero): "<<(harmFieldVec - harmVec).cwiseAbs().maxCoeff()<<std::endl;

    gradField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    gradField.set_extrinsic_field((const Eigen::VectorXd)(IE*gradFieldVec));
    //std::cout<<"exactField.extField: "<<exactField.extField<<std::endl;
    rotCogradField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    rotCogradField.set_extrinsic_field((const Eigen::VectorXd)(IE*rotCogradFieldVec));
    harmField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    harmField.set_extrinsic_field((const Eigen::VectorXd)(IE*harmVec));
   
    //triangle mesh setup
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_cartesian_field(origField,"Original Field", 0, 0, 2.0);
    viewer.set_vertex_data(vertexFunction, vertexFunction.minCoeff(), vertexFunction.maxCoeff(),"Vertex Function", 0);
    viewer.set_cartesian_field(gradField,"Gradient Field", 1, 0, 2.0);
    viewer.set_edge_data(diamondForm, diamondForm.minCoeff(), diamondForm.maxCoeff(),"Diamond Function", 0);
    viewer.set_cartesian_field(rotCogradField,"Rotated Cogradient Field", 2, 0, 2.0);
    viewer.set_cartesian_field(harmField,"Harmonic Field", 3, 0, 2.0);
    viewer.launch();
}
