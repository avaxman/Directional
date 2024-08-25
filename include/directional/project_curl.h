// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PROJECT_CURL
#define DIRECTIONAL_PROJECT_CURL

#include <Eigen/Core>
#include <vector>
#include <set>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/principal_matching.h>

namespace directional {

    inline void gradient_descent(const Eigen::SparseMatrix<double>& A,
                          const Eigen::VectorXd& gradMasses,
                          const Eigen::VectorXd& b,
                          const Eigen::VectorXd& initx,
                          const double tol,
                          const int maxIterations,
                          Eigen::VectorXd& resultx)
    {

        using namespace Eigen;
        using namespace std;

        resultx = initx;
        double denValue = b.norm();
        //cout<<"(A*initx-b).lpNorm<Infinity>(): "<<(A*initx-b).lpNorm<2>();
        for (int i=0;i<maxIterations;i++){
            VectorXd Ax=A*resultx;
            VectorXd Axb = Ax-b;
            double currValue = Axb.norm()/denValue;
            // cout<<"currValue: "<<currValue<<endl;
            if (currValue<tol)
                break;

            VectorXd gradient = gradMasses.array()*Axb.array();
            VectorXd Ag = A*gradient;

            double alpha = (Ag.dot(resultx))/(Ag.dot(gradient));
            resultx = resultx - alpha*gradient;
            //cout<<"(A*resultx-b).lpNorm<Infinity>(): "<<(A*resultx-b).lpNorm<Infinity>();
        }

    }


    inline void project_curl(const TriMesh& mesh,
                                 const Eigen::VectorXi& bc,
                                 const Eigen::MatrixXd& b,
                                 const PVData& pvData,
                                 const double smoothCoeff,
                                 Eigen::MatrixXd& rawField){

        using namespace Eigen;
        using namespace std;


        VectorXi matching;//, singVertices,singIndices, combedMatching;

        VectorXd curl, effort;
        directional::principal_matching(rawField);

        MatrixXd prevRawField=rawField;
        VectorXd rawFieldReducedVec(2*pvglData.d*F.rows());
        VectorXd rawFieldVec(2*pvglData.N*F.rows());
        VectorXi constinFace(bc.size());
        for (int i=0;i<F.rows();i++){
            for (int j=0;j<pvglData.N;j++){
                RowVector3d currVector3D = rawField.block(i, 3*j,1,3);
                rawFieldVec.segment(2*pvglData.N*i+2*j,2)<<B1.row(i).dot(currVector3D), B2.row(i).dot(currVector3D);
            }
        }

        for (int i=0;i<bc.size();i++){
            double maxDot=-32767.0;
            //cout<<"b.row(i): "<<b.row(i)<<endl;
            for (int j=0;j<pvglData.N;j++){
                RowVector3d currVector3D = rawField.block(bc(i), 3*j,1,3).normalized();
                //cout<<"currVector3D: "<<currVector3D<<endl;
                double currDot = currVector3D.dot(b.row(i));
                if (currDot>maxDot){
                    maxDot=currDot;
                    constinFace(i)=j;
                }
            }
            //cout<<"maxDot: "<<maxDot<<endl;
        }

        rawFieldReducedVec=pvglData.bigInvSymmMat*rawFieldVec;

        SparseMatrix<double> C,CNorm;
        VectorXd MeArray, dualWeightsArray;
        curl_matrix(V, F, EV, EF,  pvglData.N, matching, B1, B2,C, CNorm, MeArray,dualWeightsArray);

        SparseMatrix<double> CReduced=C*pvglData.bigSymmMat;
        SparseMatrix<double> CRedTrans=CReduced.transpose();

        SparseMatrix<double> CNormReduced=CNorm*pvglData.bigSymmMat;
        SparseMatrix<double> CNormRedTrans=CNorm.transpose();

        //cout<<"(C*rawFieldVec).lpNorm<Infinity>() before: "<<(C*rawFieldVec).lpNorm<Infinity>()<<endl;
        VectorXd x0 = rawFieldReducedVec;
        VectorXd x = x0;
        SparseMatrix<double> eyeMat,dualWeightsMat;
        igl::speye(x0.size(),x0.size(), eyeMat);
        SparseMatrix<double> MeMat, MfMat, MfConstMat;
        igl::diag(MeArray/MeArray.sum(), MeMat);
        igl::diag(dualWeightsArray/dualWeightsArray.sum(), dualWeightsMat);
        igl::diag(pvglData.MfArray/pvglData.MfArray.sum(), MfMat);
        VectorXd MfNConstArray=VectorXd::Zero(2*pvglData.N*F.rows());
        for (int i=0;i<bc.rows();i++)
            for (int j=0;j<2;j++)
                MfNConstArray(2*pvglData.N*bc(i)+2*constinFace(i)+j)=pvglData.MfNArray(2*pvglData.N*bc(i)+2*constinFace(i)+j);

        //creating orthogonal energy
        SparseMatrix<double> orthAlignMat;
        vector<Triplet<double> > orthAlignTris;
        for (int i=0;i<bc.size();i++){
            RowVector2d currVector2D; currVector2D<<-B2.row(bc(i)).dot(b.row(i)), B1.row(bc(i)).dot(b.row(i));  //dot product with the orthogonal vector
            for (int j=0;j<2;j++)
                orthAlignTris.push_back(Triplet<double>(i, 2*pvglData.N*bc(i)+2*constinFace(i)+j,currVector2D(j)));
        }

        orthAlignMat.resize(bc.rows(), 2*pvglData.N*F.rows());
        orthAlignMat.setFromTriplets(orthAlignTris.begin(), orthAlignTris.end());

        //cout<<"orthAlignMat*rawFieldVec: "<<orthAlignMat*rawFieldVec<<endl;

        VectorXd MfOrthAlignArray(bc.rows());
        for (int i=0;i<bc.size();i++)
            MfOrthAlignArray(i)=pvglData.MfArray(2*pvglData.d*bc(i));

        SparseMatrix<double> MfOrthAlignMat;
        igl::diag(MfOrthAlignArray/MfOrthAlignArray.sum(), MfOrthAlignMat);

        SparseMatrix<double> EAlign = pvglData.bigSymmMat.transpose()*orthAlignMat.transpose()*MfOrthAlignMat*orthAlignMat*pvglData.bigSymmMat;

        cout<<"Max Curl before: "<<(CNormReduced*x).lpNorm<Infinity>() <<endl;
        //cout<<"Align energy before: "<<x.transpose()*EAlign*x<<endl;

        SparseMatrix<double> CEnergyMat=CReduced.transpose()*MeMat*CReduced;

        for (double rho = 1.0; rho>1e-5;rho/=1.5){
            SparseMatrix<double> currMat = CEnergyMat + (MfMat+pvglData.alignCoeff*EAlign)*rho;
            VectorXd rhs = rho*(MfMat*x);
            VectorXd prevx=x;
            gradient_descent(currMat, pvglData.MfInvArray, rhs, prevx, 10e-4, 50, x);

            //cout<<"(C*rawFieldVec).lpNorm<Infinity>() inside iteration: "<<(C*pvglData.bigSymmMat*x).lpNorm<Infinity>()<<endl;
            //ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
            //cg.compute(currMat);
            //x=cg.solveWithGuess(rhs,prevx);
        }

        cout<<"Max Curl after: "<<(CNormReduced*x).lpNorm<Infinity>() <<endl;
        //cout<<"Align energy after: "<<x.transpose()*EAlign*x<<endl;

        //normalizing result
        rawFieldVec=pvglData.bigSymmMat*x;
        rawField.resize(rawField.rows(), rawField.cols());
        double totalLengthSum=0.0;
        for (int i=0;i<F.rows();i++)
            for (int j=0;j<pvglData.N;j++){
                rawField.block(i, 3*j,1,3)=rawFieldVec(2*pvglData.N*i+2*j)*B1.row(i)+rawFieldVec(2*pvglData.N*i+2*j+1)*B2.row(i);
                totalLengthSum+=rawField.block(i, 3*j,1,3).norm();
            }

        rawField.array()/=totalLengthSum/(double)(F.rows()*pvglData.N);
        rawFieldVec.array()/=totalLengthSum/(double)(F.rows()*pvglData.N);
        cout<<"(C*rawFieldVec).lpNorm<Infinity>() after: "<<(C*rawFieldVec).lpNorm<Infinity>()<<endl;
    }
};


#endif
