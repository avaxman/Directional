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
Eigen::VectorXd gradFieldVec, rotCogradFieldVec, harmFieldVec, vertexFunction, edgeFunction;
directional::DirectionalViewer viewer;

int main()
{
    directional::readOFF(TUTORIAL_DATA_PATH "/1146164.off",mesh);
    ftb.init(mesh);
    
    //Must use intrinsic since otherwise the harmonic field will have spurious normal components
    Eigen::SparseMatrix<double> G = directional::conf_gradient_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(mesh, Eigen::VectorXi(), true);
    Eigen::SparseMatrix<double> Mx = directional::face_vectors_mass_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> iMx = directional::face_vectors_mass_matrix_2D<double>(mesh, true, true);
    //Eigen::SparseMatrix<double> Mc = directional::edge_diamond_mass_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(mesh);
    
    Eigen::MatrixXd harmBasis;
    int bettiNumber = mesh.EV.rows() - (mesh.V.rows()-1) - (mesh.F.rows()-1);
    std::cout<<"Extracting the basis of harmonic fields..."<<std::endl;
    directional::cohomology_basis(G, C, Mx,  bettiNumber, harmBasis);
    
    Eigen::VectorXd harmGT = harmBasis.col(0);
    Eigen::RowVector3d COM = mesh.V.colwise().mean();
    Eigen::VectorXd vertexVec(mesh.dcel.vertices.size()), midEdgeVec(mesh.dcel.edges.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        vertexVec[i] = (mesh.V.row(i)-mesh.V.row(296)).norm();
    
    for (int i=0;i<mesh.dcel.edges.size();i++)
        midEdgeVec[i] = sin(mesh.midEdges(i,0)/2.0)+sin(mesh.midEdges(i,1)/2.0)+sin(mesh.midEdges(i, 2)/2.0);
    
    //Doing an artifical composition and decomposition. For results to be meaningful, the fields must be of more or less the same magnitude, so normalizing them.
    Eigen::VectorXd gradientGT = G*vertexVec;
    gradientGT.array()/=sqrt((gradientGT.transpose()*Mx*gradientGT).coeff(0,0));
    Eigen::VectorXd rotCogradientGT = -iMx*C.adjoint()*midEdgeVec;    //This is equivalent to J*Ge*midEdgeVec
    rotCogradientGT.array()/=sqrt((rotCogradientGT.transpose()*Mx*rotCogradientGT).coeff(0,0));
    Eigen::VectorXd origFieldVec = harmGT + gradientGT +  rotCogradientGT;
    //Using identity matrix since the convention is that Hodge decomposition in face based fields provide the non-conforming mid-edge based function
    Eigen::SparseMatrix<double> I(mesh.EV.rows(), mesh.EV.rows());
    I.setIdentity();
    directional::hodge_decomposition<double>(G, C, Mx, I, origFieldVec, vertexFunction, gradFieldVec, rotCogradFieldVec, edgeFunction, harmFieldVec);
    std::cout<<"Exact reproduction (numerically zero): "<<(gradFieldVec -gradientGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Coexact reproduction (numerically zero): "<<(rotCogradFieldVec - rotCogradientGT).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"Harmonic reproduction (numerically zero): "<<(harmFieldVec - harmGT).cwiseAbs().maxCoeff()<<std::endl;
    
    //Visualization of the extrinsic fields
    Eigen::MatrixXd gradField(mesh.F.rows(),3), origField(mesh.F.rows(),3), rotCogradField(mesh.F.rows(),3), harmField(mesh.F.rows(),3);
    gradFieldVec = IE*gradFieldVec;
    rotCogradFieldVec = IE*rotCogradFieldVec;
    harmFieldVec = IE*harmFieldVec;
    origFieldVec = IE*origFieldVec;
    for (int i=0;i<mesh.F.rows();i++){
        gradField.row(i) = gradFieldVec.segment(3*i,3).transpose();
        rotCogradField.row(i) = rotCogradFieldVec.segment(3*i,3).transpose();
        harmField.row(i) = harmFieldVec.segment(3*i,3).transpose();
        origField.row(i) = origFieldVec.segment(3*i,3).transpose();
    }
    
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_surface_vertex_data(vertexFunction, "Vertex Function", 0)->setEnabled(false);
    viewer.set_surface_edge_data(edgeFunction, "Edge Function", 0)->setEnabled(false);
    viewer.set_raw_field(mesh.barycenters, origField, "Original Field", 0, 20.0*mesh.avgEdgeLength);
    viewer.set_raw_field(mesh.barycenters, gradField, "Gradient Field", 1, 20.0*mesh.avgEdgeLength);
    viewer.toggle_raw_field(false, 1);
    viewer.set_raw_field(mesh.barycenters, rotCogradField, "Rotated Cogradient Field", 2, 20.0*mesh.avgEdgeLength);
    viewer.toggle_raw_field(false, 2);
    viewer.set_raw_field(mesh.barycenters, harmField, "Harmonic Field", 3, 20.0*mesh.avgEdgeLength);
    viewer.toggle_raw_field(false, 3);
    viewer.launch();
}
