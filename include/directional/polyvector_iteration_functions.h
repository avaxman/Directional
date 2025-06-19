// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_ITERATION_FUNCTIONS_H
#define DIRECTIONAL_POLYVECTOR_ITERATION_FUNCTIONS_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/project_curl.h>
#include <directional/principal_matching.h>
#include <directional/sparse_diagonal.h>
#include <directional/PolyVectorData.h>

namespace directional{

typedef std::function<CartesianField(const CartesianField&, const PolyVectorData&)> PvIterationFunction;


//Normalizes the field to have unit length on all its vectors
CartesianField hard_normalization(const CartesianField& pvField, const PolyVectorData& pvData){
    CartesianField rawField, normPvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
    directional::raw_to_polyvector(rawField, normPvField);
    return normPvField;
}


CartesianField hard_rosy(const CartesianField& pvField, const PolyVectorData& pvData){
    Eigen::MatrixXcd rosyField = pvField.get_complex_intrinsic_field();
    rosyField.block(0,1,rosyField.rows(), rosyField.cols()-1).setZero();
    rosyField.col(0) = rosyField.col(0).array() / rosyField.col(0).array().abs();
    //std::cout<<"rosyField: "<<rosyField<<std::endl;
    CartesianField hardRosyField = pvField;
    hardRosyField.set_intrinsic_field(rosyField);
    return hardRosyField;
}



CartesianField soft_rosy(const CartesianField& pvField, const PolyVectorData& pvData){
    Eigen::MatrixXcd rosyField = pvField.get_complex_intrinsic_field();
    rosyField.block(0,1,rosyField.rows(), rosyField.cols()-1).setZero();
    rosyField.col(0) = rosyField.col(0).array() / rosyField.col(0).array().abs();
    //doing closed-form implicit step
    CartesianField softRosyField = pvField;
    Eigen::MatrixXcd interpField = (pvField.get_complex_intrinsic_field().array() + 2.0*pvData.currImplicitCoeff*rosyField.array())/(1+2.0*pvData.currImplicitCoeff);
    softRosyField.set_intrinsic_field(interpField);
    return softRosyField;
}

CartesianField curl_projection(const CartesianField& pvField, const PolyVectorData& pvData){
    assert(pvData.tb->discTangType()==directional::discTangTypeEnum::FACE_SPACES && "Projecting curl only works for face-based fields for now!");
    CartesianField rawField, curlFreeFieldRaw, curlFreeFieldPv;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, false);
    directional::principal_matching(rawField);
    project_curl(rawField, Eigen::VectorXi(), Eigen::MatrixXd(), curlFreeFieldRaw);
    directional::raw_to_polyvector(curlFreeFieldRaw,  curlFreeFieldPv);
    return curlFreeFieldPv;
}


Eigen::RowVectorXd project_on_quadric(const Eigen::RowVectorXd& y0, const Eigen::MatrixXd& H){
    // Step 1: Perform eigen-decomposition of H (since it's symmetric)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(H);
    assert("directional:project_on_quadric(): Eigen decomposition failed!" && eigensolver.info() == Eigen::Success);
    
    Eigen::MatrixXd U = eigensolver.eigenvectors();
    Eigen::VectorXd sigma = eigensolver.eigenvalues();
    Eigen::VectorXd z = U.transpose() * y0.transpose();
    double lambda = 0.0; // initial guess
    int maxIters = 100;
    double tol = 10e-9;
    for (int iter = 0; iter < maxIters; ++iter) {
        double f = 0.0;
        double df = 0.0;
        for (int i = 0; i < z.size(); ++i) {
            double zi = z[i], si = sigma[i];
            double denom = 1.0 + lambda * si;
            if (std::abs(denom) < 1e-14) continue; // skip unstable terms
            
            double denom2 = denom * denom;
            double denom3 = denom2 * denom;
            
            f  += zi * zi * si / denom2;
            df += -2.0 * zi * zi * si * si / denom3;
        }
        std::cout<<"f: "<<f<<std::endl;
        if (std::abs(f) < tol) break;
        
        double step = f / df;
        lambda -= step;
        
        if (std::abs(step) < tol) break;
        
        if (std::isnan(lambda) || std::isinf(lambda))
            throw std::runtime_error("Newton-Raphson diverged");
    }
    
    //checking lambda is correct
    double lambdasum = 0.0;
    for (int i=0;i<y0.size();i++)
        lambdasum+=z(i)*z(i)*sigma(i)/((1+sigma(i)*lambda)*(1+sigma(i)*lambda));
    std::cout<<"lambdasum: "<<lambdasum<<std::endl;
    Eigen::ArrayXd denom = (1.0 + lambda * sigma.array());
    Eigen::VectorXd z_scaled = (z.array() / denom).matrix();
    Eigen::VectorXd y = U * z_scaled;                                 // back to original coordinates
    
    //checking
    std::cout<<"y0.transpose() * H * y0: "<<y0 * H * y0.transpose()<<std::endl;
    std::cout<<"y.transpose() * H * y: "<<y.transpose() * H * y<<std::endl;
    
    return y.transpose(); // Return RowVectorXd
}


//Projecting a 2^2 Rosy field to a conjugate field
CartesianField conjugate(const CartesianField& pvField, const PolyVectorData& pvData){
    assert("directional::conjugate(): This method only works on symmetric 2^2 fields on faces!" && pvField.N==4 && pvData.signSymmetry && pvField.tb->discTangType()==discTangTypeEnum::FACE_SPACES);
    CartesianField rawField, conjugatePvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
    Eigen::MatrixXd extField = rawField.extField;
    PCFaceTangentBundle* tb =(PCFaceTangentBundle*)pvField.tb;
    for (int i=0;i<tb->mesh->F.rows();i++){
        Eigen::Matrix3d G1 =tb->mesh->facePrincipalCurvatures(i,0)*tb->mesh->minFacePrincipalDirectionals.row(i).transpose()*tb->mesh->minFacePrincipalDirectionals.row(i);
        Eigen::Matrix3d G2 =tb->mesh->facePrincipalCurvatures(i,1)*tb->mesh->maxFacePrincipalDirectionals.row(i).transpose()*tb->mesh->maxFacePrincipalDirectionals.row(i);
        Eigen::Matrix<double, 6,6> H; H.setZero();
        H.block(0,3,3,3) = G1;
        H.block(3,0,3,3) = G2;
        Eigen::RowVectorXd y0(6); y0<<rawField.extField.row(i).head(6);
        extField.row(i).head(6)<<project_on_quadric(y0, H);
    }
    rawField.set_extrinsic_field(extField);
    directional::raw_to_polyvector(rawField, conjugatePvField);
    return conjugatePvField;
}

}


#endif
