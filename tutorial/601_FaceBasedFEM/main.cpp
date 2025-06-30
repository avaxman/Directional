#include <iostream>
#include <Eigen/Core>
#include <directional/readOBJ.h>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>
#include <directional/directional_viewer.h>
#include <directional/gradient_matrices.h>
#include <directional/mass_matrices.h>
#include <directional/curl_matrices.h>
#include <directional/div_matrices.h>

directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField exactField, coexactField;
directional::DirectionalViewer viewer;


double maxAbsValue(const Eigen::SparseMatrix<double>& A) {
    double maxVal = 0.0;
    for (int k = 0; k < A.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(A, k); it; ++it) {
            maxVal = std::max(maxVal, std::abs(it.value()));
        }
    }
    return maxVal;
}


int main()
{
    
    directional::readOBJ(TUTORIAL_DATA_PATH "/60745__sf.obj",mesh);
    //directional::readOFF(TUTORIAL_DATA_PATH "/Flat.off",mesh);
    ftb.init(mesh);
    
    Eigen::VectorXd confVec(mesh.dcel.vertices.size()), nonConfVec(mesh.dcel.edges.size());
    for (int i=0;i<mesh.dcel.vertices.size();i++)
        confVec[i] = 10.0*sin(mesh.V(i,2)/20.0);//*cos(mesh.V(i,1)/20.0);
    for (int i=0;i<mesh.dcel.edges.size();i++)
        nonConfVec[i] = 5.0*mesh.midEdges(i,0);
    
    Eigen::SparseMatrix<double> Gv = directional::conf_gradient_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> Ge = directional::non_conf_gradient_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> J =  directional::face_vector_rotation_matrix_2D<double>(mesh, true);
    
    Eigen::SparseMatrix<double> C = directional::curl_matrix_2D<double>(mesh, Eigen::VectorXi(), true);
    Eigen::SparseMatrix<double> D = directional::div_matrix_2D<double>(mesh, true);
    
    //Debugging
    Eigen::SparseMatrix<double> Mx = directional::face_vectors_mass_matrix_2D<double>(mesh, true);
    Eigen::SparseMatrix<double> iMx = directional::face_vectors_mass_matrix_2D<double>(mesh, true, true);
    
    
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(mesh);

    std::cout<<"C: "<<std::endl;
    int numMembers = 10;
    int currMember = 0;
    for (int k = 0; k < C.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(C, k); it; ++it){
            std::cout<<it.row()<<","<<it.col()<<","<<it.value()<<std::endl;
            if (currMember>numMembers)
                break;
            currMember++;
        }
        if (currMember>numMembers)
            break;
        
    }
    
    currMember=0;
    
    
    std::cout<<"C2: "<<std::endl;
    Eigen::SparseMatrix<double> C2 = -(J*Ge).adjoint()*Mx;
    for (int k = 0; k < C2.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(C2, k); it; ++it){
            std::cout<<it.row()<<","<<it.col()<<","<<it.value()<<std::endl;
            if (currMember>numMembers)
                break;
            currMember++;
        }
        if (currMember>numMembers)
            break;
    }
    
    Eigen::SparseMatrix<double> exactRelation = Gv.adjoint()*Mx*J*Ge;
    Eigen::SparseMatrix<double> coexactRelation = (J*Ge).adjoint()*Mx*Gv;
    std::cout<<"max(abs(exactRelation): "<<maxAbsValue(exactRelation)<<std::endl;
    std::cout<<"max(abs(coexactRelation): "<<maxAbsValue(coexactRelation)<<std::endl;
    
    //std::cout<<"difference between C and (JG^T)M: "<<maxAbsValue(C-C2)<<std::endl;
    //std::cout<<"difference between iMx*C^T and J*Ge: "<<maxAbsValue(iMx*C.transpose()-J*Ge)<<std::endl;
    
    //demonstrating the exact sequences
    std::cout<<"max abs curl of exact field (should be numerically zero): "<<(C*Gv*confVec).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"max abs curl2 of exact field (should be numerically zero): "<<(C2*Gv*confVec).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"max abs divergence of coexact field (should be numerically zero): "<<(D*J*Ge*nonConfVec).cwiseAbs().maxCoeff()<<std::endl;
    std::cout<<"max abs divergence2 of coexact field (should be numerically zero): "<<(D*iMx*C.transpose()*nonConfVec).cwiseAbs().maxCoeff()<<std::endl;
    
    Eigen::MatrixXd gradRawField(mesh.F.rows(), 3), rotCogradRawField(mesh.F.rows(), 3);
    Eigen::VectorXd gradFieldVec = IE*Gv*confVec;
    Eigen::VectorXd rotCogradFieldVec = IE*(-iMx*C.adjoint())*nonConfVec;
    for (int i=0;i<mesh.F.rows();i++){
        gradRawField.row(i)<<gradFieldVec.segment(3*i,3).transpose();
        rotCogradRawField.row(i)<<rotCogradFieldVec.segment(3*i,3).transpose();
    }
    
    coexactField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, 1);
    coexactField.set_extrinsic_field(rotCogradFieldVec);
    
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_edge_data(nonConfVec, nonConfVec.minCoeff(), nonConfVec.maxCoeff(),"Non-Conforming Function");
    //viewer.set_raw_field(mesh.barycenters, rotCogradRawField, 20.0*mesh.avgEdgeLength, "Rot. Cogradient field", 0);
    viewer.set_cartesian_field(coexactField, "Rot. Cogradient field", 0, 0, 10.0);
    viewer.set_vertex_data(confVec, confVec.minCoeff(), confVec.maxCoeff(),"Conforming Function");
    viewer.set_raw_field(mesh.barycenters, gradRawField, 2.0*mesh.avgEdgeLength, "Gradient field", 1);
    viewer.toggle_raw_field(false, 1);
    viewer.launch();
}
