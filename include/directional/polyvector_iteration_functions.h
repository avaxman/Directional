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
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/project_curl.h>
#include <directional/principal_matching.h>
#include <directional/sparse_diagonal.h>
#include <directional/PolyVectorData.h>

namespace directional{


//This file contains a collection of functions that can be used in the "local" or "projection" steps of the extended polyvector algorithm.

typedef std::function<CartesianField(const CartesianField&, const PolyVectorData&)> PvIterationFunction;
typedef std::function<bool(const CartesianField&, const PolyVectorData&)>           PvTerminationFunction;


//Normalizes the field to have unit length on all its vectors
CartesianField hard_normalization(const CartesianField& pvField, const PolyVectorData& pvData){
    CartesianField rawField, normPvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
    directional::raw_to_polyvector(rawField, normPvField);
    return normPvField;
}

//Projects the field into the nearest RoSy field
CartesianField hard_rosy(const CartesianField& pvField, const PolyVectorData& pvData){
    Eigen::MatrixXcd rosyField = pvField.get_complex_intrinsic_field();
    rosyField.block(0,1,rosyField.rows(), rosyField.cols()-1).setZero();
    rosyField.col(0) = rosyField.col(0).array() / rosyField.col(0).array().abs();
    CartesianField hardRosyField = pvField;
    hardRosyField.set_intrinsic_field(rosyField);
    return hardRosyField;
}


//A single implicit step (with the pvData state coefficients) that makes the current field more RoSy.
CartesianField soft_rosy(const CartesianField& pvField, const PolyVectorData& pvData){
    Eigen::MatrixXcd rosyField = pvField.get_complex_intrinsic_field();
    rosyField.block(0,1,rosyField.rows(), rosyField.cols()-1).setZero();
    rosyField.col(0) = rosyField.col(0).array() / rosyField.col(0).array().abs();
    //doing closed-form implicit step
    CartesianField softRosyField = pvField;
    Eigen::MatrixXcd interpField = (pvField.get_complex_intrinsic_field().array() + 2.0*pvData.currImplicitFactor*rosyField.array())/(1+2.0*pvData.currImplicitFactor);
    softRosyField.set_intrinsic_field(interpField);
    return softRosyField;
}

//Projects the field onto the nearest curl-free field.
CartesianField curl_projection(const CartesianField& pvField, const PolyVectorData& pvData){
    assert(pvData.tb->discTangType()==directional::discTangTypeEnum::FACE_SPACES && "Projecting curl only works for face-based fields for now!");
    CartesianField rawField, curlFreeFieldRaw, curlFreeFieldPv;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0);
    directional::principal_matching(rawField);
    project_curl(rawField, Eigen::VectorXi(), Eigen::MatrixXd(), curlFreeFieldRaw);
    directional::raw_to_polyvector(curlFreeFieldRaw,  curlFreeFieldPv);
    return curlFreeFieldPv;
}

//Projects a vector on a quadric (used for the conjugacy projection)
Eigen::RowVectorXd project_on_quadric(const Eigen::RowVectorXd& y0, const Eigen::MatrixXd& H){
    // Step 1: Perform eigen-decomposition of H (since it's symmetric)
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(0.5*(H+H.transpose()));
    assert(eigensolver.info() == Eigen::Success && "directional:project_on_quadric(): Eigen decomposition failed!");
    
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
        if (std::abs(f) < tol) break;
        
        double step = f / df;
        lambda -= step;
        
        if (std::abs(step) < tol) break;
        
        if (std::isnan(lambda) || std::isinf(lambda))
            throw std::runtime_error("Newton-Raphson diverged");
    }
    
    Eigen::ArrayXd denom = (1.0 + lambda * sigma.array());
    Eigen::VectorXd z_scaled = (z.array() / denom).matrix();
    Eigen::VectorXd y = U * z_scaled;                                 // back to original coordinates
    
    return y.transpose(); // Return RowVectorXd
}


//Projecting a 2^2 Rosy field to a conjugate field
CartesianField conjugate(const CartesianField& pvField, const PolyVectorData& pvData){
    assert(pvField.N==4 && pvData.signSymmetry && pvField.tb->discTangType()==discTangTypeEnum::VERTEX_SPACES&& "directional::conjugate(): This method only works on symmetric 2^2 fields on vertices!");
    
    CartesianField rawField, conjugatePvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0);
    Eigen::MatrixXd extField = rawField.extField;
    IntrinsicVertexTangentBundle* tb =(IntrinsicVertexTangentBundle*)pvField.tb;
    Eigen::VectorXd conjugacy(tb->mesh->V.rows());
    Eigen::VectorXd conjugacyBefore(tb->mesh->V.rows());
    for (int i=0;i<tb->mesh->V.rows();i++){
        Eigen::Matrix3d G1 =tb->mesh->vertexPrincipalCurvatures(i,0)*tb->mesh->minVertexPrincipalDirections.row(i).transpose()*tb->mesh->minVertexPrincipalDirections.row(i);
        Eigen::Matrix3d G2 =tb->mesh->vertexPrincipalCurvatures(i,1)*tb->mesh->maxVertexPrincipalDirections.row(i).transpose()*tb->mesh->maxVertexPrincipalDirections.row(i);
        Eigen::Matrix<double, 6,6> H; H.setZero();
        H.block(0,3,3,3) = G1;
        H.block(3,0,3,3) = G2;
        Eigen::RowVectorXd y0(6); y0<<rawField.extField.row(i).head(6);
        conjugacyBefore(i) = extField.row(i).head(3)*tb->mesh->Sv[i]*extField.row(i).segment(3,3).transpose();
        if ((std::abs(tb->mesh->vertexPrincipalCurvatures(i,0))>10e-8)||(std::abs(tb->mesh->vertexPrincipalCurvatures(i,1))>10e-8))
            extField.row(i).head(6)<<project_on_quadric(y0, H);
        else
            extField.row(i).head(6) = y0;
        extField.row(i).tail(6) = - extField.row(i).head(6);
        //checking conjugacy
        conjugacy(i) = (extField.row(i).head(3)*tb->mesh->Sv[i]*extField.row(i).segment(3,3).transpose()).coeff(0,0);
    }
    //std::cout<<"conjugacy: "<<conjugacy.head(20)<<std::endl;
    if (pvData.verbose){
        std:: cout<<"conjugacy before conjugate(): "<<conjugacyBefore.cwiseAbs().maxCoeff()<<std::endl;
        std:: cout<<"conjugacy after  conjugate(): "<<conjugacy.cwiseAbs().maxCoeff()<<std::endl;
    }
    rawField.set_extrinsic_field(extField);
    directional::raw_to_polyvector(rawField, conjugatePvField);
    //std::cout<<"rawField.row(0)"<<rawField.extField.row(0)<<std::endl;
    //std::cout<<"conjugatePvField.row(0)"<<conjugatePvField.intField.row(0)<<std::endl;
    CartesianField rawField2;
    directional::polyvector_to_raw(conjugatePvField, rawField2, pvData.N%2==0);
    //std::cout<<"rawField2.row(0)"<<rawField2.extField.row(0)<<std::endl;
    //std::cout<<"rawField2 - rawField"<<(rawField.extField-rawField2.extField).maxCoeff()<<std::endl;
    return conjugatePvField;
}

bool conjugate_termination(const CartesianField& pvField, const PolyVectorData& pvData){
    assert(pvField.N==4 && pvData.signSymmetry && pvField.tb->discTangType()==discTangTypeEnum::VERTEX_SPACES&& "directional::conjugate_termination(): This method only works on symmetric 2^2 fields on vertices!");
    static double tolerance = 1e-2;
    CartesianField rawField, conjugatePvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0);
    //std::cout<<"rawField.row(0)"<<rawField.extField.row(0)<<std::endl;
    Eigen::MatrixXd extField = rawField.extField;
    IntrinsicVertexTangentBundle* tb =(IntrinsicVertexTangentBundle*)pvField.tb;
    Eigen::VectorXd conjugacy(tb->mesh->V.rows());
    for (int i=0;i<tb->mesh->V.rows();i++){
        Eigen::Matrix3d G1 =tb->mesh->vertexPrincipalCurvatures(i,0)*tb->mesh->minVertexPrincipalDirections.row(i).transpose()*tb->mesh->minVertexPrincipalDirections.row(i);
        Eigen::Matrix3d G2 =tb->mesh->vertexPrincipalCurvatures(i,1)*tb->mesh->maxVertexPrincipalDirections.row(i).transpose()*tb->mesh->maxVertexPrincipalDirections.row(i);
        Eigen::Matrix<double, 6,6> H; H.setZero();
        H.block(0,3,3,3) = G1;
        H.block(3,0,3,3) = G2;
        //conjugacy(i) = (extField.row(i).head(3)*tb->mesh->Sv[i]*extField.row(i).segment(3,3).transpose()).cwiseAbs().coeff(0,0);
        conjugacy(i) = (extField.row(i).head(6)*H*extField.row(i).head(6).transpose()).cwiseAbs().coeff(0,0);
        
    }
    if (conjugacy.maxCoeff()>tolerance){
        if (pvData.verbose)
            std::cout<<"Maximum conjugacy is "<<conjugacy.maxCoeff()<<" and above termination threshold "<<tolerance<<std::endl;
        return false;
    } else {
        if (pvData.verbose)
            std:: cout<<"conjugacy at termination: "<<conjugacy.cwiseAbs().maxCoeff()<<std::endl;
        return true;
    }
    
}

}


#endif
