// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_FIELD_H
#define DIRECTIONAL_POLYVECTOR_FIELD_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>
#include <iostream>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/principal_matching.h>
#include <directional/sparse_diagonal.h>
#include <directional/polyvector_iteration_functions.h>
#include <directional/PolyVectorData.h>

namespace directional
{


// Precalculate the operators needed for PolyVector computation according to the user-prescribed parameters. Must be called whenever any of them changes
// Input:
//  pvField: POLYVECTOR_FIELD cartesian field initalized with the tangent bundle
//  pvData: the details of the algorithm
//
// Output:
//  pvField: the computer field (returned in the input parameter)
//  pvData:  Updated structure with all operators
inline void polyvector_precompute(directional::CartesianField& pvField,
                                  PolyVectorData& pvData)
{
    
    using namespace std;
    using namespace Eigen;
    
    pvField.init(*(pvData.tb), fieldTypeEnum::POLYVECTOR_FIELD, pvData.N);
    
    //Building the smoothness matrices, with an energy term for each inner edge and degree
    int rowCounter=0;
    std::vector< Triplet<complex<double> > > dTriplets, WTriplets;
    if (pvData.N%2!=0) pvData.signSymmetry=false;  //it has to be for odd N
    
    pvData.totalSmoothWeight = pvField.tb->connectionMass.sum();
    
    //vector<Triplet<complex<double>>> WSmoothTriplets, MTriplets;
    for (int n = 0; n < pvData.N; n++)
    {
        for (int i=0;i<pvField.tb->adjSpaces.rows();i++)
        {
            if ((pvField.tb->adjSpaces(i,0)==-1)||(pvField.tb->adjSpaces(i,1)==-1))
                continue;  //boundary edge
            
            // differential matrix between two tangent spaces
            dTriplets.push_back(Triplet<complex<double> >(rowCounter, n*pvField.intField.rows()+pvField.tb->adjSpaces(i,0), pow(pvField.tb->connection(i),pvData.N-n)));
            dTriplets.push_back(Triplet<complex<double> >(rowCounter, n*pvField.intField.rows()+pvField.tb->adjSpaces(i,1), -1.0));
            
            //stiffness weights
            //WSmoothTriplets.push_back(Triplet<complex<double> >(rowCounter, rowCounter, pvField.tb->connectionMass(i)));
            rowCounter++;
        }
        
        //for (int i=0;i<pvField.intField.rows();i++)
        //    MTriplets.push_back(Triplet<complex<double>>(n*pvField.intField.rows()+i, n*pvField.intField.rows()+i, pvField.tb->tangentSpaceMass(i)));
    }
    
    pvData.smoothMat.resize(rowCounter, pvData.N*pvField.intField.rows());
    pvData.smoothMat.setFromTriplets(dTriplets.begin(), dTriplets.end());
    
    std::vector<Eigen::SparseMatrix<double>> WVector, MVector;
    MatrixXi blkDiagIndices(pvData.N,2);
    for (int n=0; n<pvData.N;n++) {
        WVector.push_back(pvField.tb->connectionMass);
        MVector.push_back(pvField.tb->tangentSpaceMass);
        blkDiagIndices.row(n) << n, n;
    }
    SparseMatrix<double> WSmooth, M;
    directional::sparse_diagonal(WVector, WSmooth);
    directional::sparse_diagonal(MVector, M);
    pvData.WSmooth = WSmooth.cast<std::complex<double>>();
    pvData.M = M.cast<std::complex<double>>();
    
    //pvData.WSmooth.resize(rowCounter, rowCounter);
    //pvData.WSmooth.setFromTriplets(WSmoothTriplets.begin(), WSmoothTriplets.end());
    
    //pvData.M.resize(pvData.N*pvField.intField.rows(), pvData.N*pvField.intField.rows());
    //pvData.M.setFromTriplets(MTriplets.begin(), MTriplets.end());
    
    //creating reduction transformation
    VectorXi numSpaceConstraints = Eigen::VectorXi::Zero(pvField.intField.rows());
    int realN = (pvData.signSymmetry ? pvData.N/2 : pvData.N);
    realN = (pvData.wRoSy < 0.0 ? 1 : realN);
    //MatrixXcd faceConstraints(F.rows(),realN);
    std::vector<MatrixXcd> localSpaceReducMats; localSpaceReducMats.resize(pvField.intField.rows());
    std::vector<VectorXcd> localSpaceReducRhs;  localSpaceReducRhs.resize(pvField.intField.rows());
    
    for (int i=0;i<pvField.intField.rows();i++){
        localSpaceReducMats[i]=MatrixXcd::Identity(realN,realN);
        localSpaceReducRhs[i]=VectorXcd::Zero(realN);
    }
    
    /*************Hard-constraint reduction matrices******************/
    MatrixXd constVectorsIntrinsic=pvField.tb->project_to_intrinsic(pvData.constSpaces,pvData.constVectors);
    //cout<<"constVectorsIntrinsic: "<<constVectorsIntrinsic<<endl;
    for (int i=0;i<pvData.constSpaces.size();i++){
        if (pvData.wAlignment(i)>=0.0)
            continue;  //here we only handle the reduction caused by a fixed dof
        if (numSpaceConstraints(pvData.constSpaces(i))==realN)
            continue; //overconstrained; we ignore any further constraints on that face
        
        int currSpaceNumDof = numSpaceConstraints(pvData.constSpaces(i));
        
        complex<double> constVectorComplexRaw = complex<double>(constVectorsIntrinsic(i,0),constVectorsIntrinsic(i,1));
        complex<double> constVectorComplex = (pvData.signSymmetry ? constVectorComplexRaw*constVectorComplexRaw : constVectorComplexRaw);
        constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorComplexRaw, pvData.N) :constVectorComplex);
        
        MatrixXcd singleReducMat=MatrixXcd::Zero(realN-currSpaceNumDof,realN-currSpaceNumDof-1);
        VectorXcd singleReducRhs=VectorXcd::Zero(realN-currSpaceNumDof);
        for (int j=0;j<realN-currSpaceNumDof-1;j++){
            singleReducMat(j,j)=-constVectorComplex;
            singleReducMat(j+1,j)=complex<double>(1.0,0.0);
        }
        singleReducRhs(realN-currSpaceNumDof-1) = -constVectorComplex;
        
        localSpaceReducRhs[pvData.constSpaces(i)] = localSpaceReducMats[pvData.constSpaces(i)]*singleReducRhs + localSpaceReducRhs[pvData.constSpaces(i)];
        localSpaceReducMats[pvData.constSpaces(i)] = localSpaceReducMats[pvData.constSpaces(i)]*singleReducMat;
        
        numSpaceConstraints(pvData.constSpaces(i))++;
        
        //faceConstraints(pvData.constSpaces(i), numSpaceConstraints(pvData.constSpaces(i))++) = constVectorComplex;
    }
    
    //creating the global reduction matrices
    double colCounter=0;
    pvData.reducRhs=VectorXcd::Zero(pvData.N*pvField.intField.rows());
    vector<Triplet<complex<double>>> reducMatTriplets;
    int jump = (pvData.signSymmetry ? 2 : 1);
    jump = (pvData.wRoSy < 0.0 ? pvData.N : jump);
    for (int i=0;i<pvField.intField.rows();i++){
        for (int j=0;j<pvData.N;j+=jump){
            for (int k=0;k<localSpaceReducMats[i].cols();k++)
                reducMatTriplets.push_back(Triplet<complex<double>>(j*pvField.intField.rows()+i, colCounter+k, localSpaceReducMats[i](j/jump,k)));
            
            pvData.reducRhs(j*pvField.intField.rows()+i) = localSpaceReducRhs[i](j/jump);
        }
        
        colCounter+=localSpaceReducMats[i].cols();
    }
    
    pvData.reducMat.resize(pvData.N*pvField.intField.rows(), colCounter);
    pvData.reducMat.setFromTriplets(reducMatTriplets.begin(), reducMatTriplets.end());
    
    
    /****************rotational-symmetry matrices********************/
    //TODO: use new massweights
    if (pvData.wRoSy >= 0.0){ //this is anyhow enforced, this matrix is unnecessary)
        vector<Triplet<complex<double>>> roSyTriplets;
        for (int i=pvField.intField.rows();i<pvData.N*pvField.intField.rows();i++){
            roSyTriplets.push_back(Triplet<complex<double>>(i,i,1.0));
            //WRoSyTriplets.push_back(Triplet<complex<double>>(i,i,pvField.tb->tangentSpaceMass(i%pvField.intField.rows())));
        }
        
        std::vector<Eigen::SparseMatrix<double>> WRosyVec;
        WRosyVec.push_back(Eigen::SparseMatrix<double>(pvField.intField.rows(), pvField.intField.rows()));  //for the free coefficient
        Eigen::MatrixXi WRosyIndices(pvData.N,2);
        WRosyIndices.row(0)<<0,0;
        for (int i=1;i<pvData.N;i++){
            WRosyVec.push_back(pvField.tb->tangentSpaceMass);
            WRosyIndices.row(i)<<i,i;
        }
        
        pvData.roSyMat.resize(pvData.N*pvField.intField.rows(), pvData.N*pvField.intField.rows());
        pvData.roSyMat.setFromTriplets(roSyTriplets.begin(), roSyTriplets.end());
        
        //pvData.WRoSy.resize(N*pvField.intField.rows(), N*pvField.intField.rows());
        //pvData.WRoSy.setFromTriplets(WRoSyTriplets.begin(), WRoSyTriplets.end());
        Eigen::SparseMatrix<double> WRosy;
        directional::sparse_diagonal(WRosyVec, WRosy);
        pvData.WRoSy = WRosy.cast<std::complex<double>>();
        
        pvData.totalRoSyWeight=((double)pvData.N)*pvField.tb->tangentSpaceMass.sum();
    } else {
        pvData.roSyMat.resize(0, pvData.N*pvField.intField.rows());
        pvData.WRoSy.resize(0,0);  //Even necessary?
        pvData.totalRoSyWeight=1.0;
    }
    
    /*****************Soft alignment matrices*******************/
    rowCounter=0;
    vector<Triplet<complex<double>>> alignTriplets;
    vector<VectorXcd> alignRhsList;
    vector<Triplet<complex<double>>> WAlignTriplets;
    pvData.totalConstrainedWeight=0.0;
    bool noSoftAlignment = true;
    for (int i=0;i<pvData.constSpaces.size();i++){
        if (pvData.wAlignment(i)<0.0)
            continue;  //here we only handle soft alignments
        
        noSoftAlignment=false;
        
        //complex<double> constVectorComplexSingle=std::complex<double>(pvData.constVectors.row(i).dot(mesh.Bx.row(pvData.constSpaces(i))), pvData.constVectors.row(i).dot(mesh.By.row(pvData.constSpaces(i))));
        complex<double> constVectorComplexRaw = complex<double>(constVectorsIntrinsic(i,0),constVectorsIntrinsic(i,1));
        complex<double> constVectorComplex = (pvData.signSymmetry ? constVectorComplexRaw*constVectorComplexRaw : constVectorComplexRaw);
        constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorComplexRaw, pvData.N) : constVectorComplex);
        numSpaceConstraints(pvData.constSpaces(i))++;
        //faceConstraints(pvData.constSpaces(i), numSpaceConstraints(pvData.constSpaces(i))++) = constVectorComplex;
        
        MatrixXcd singleReducMat=MatrixXcd::Zero(realN,realN-1);
        VectorXcd singleReducRhs=VectorXcd::Zero(realN);
        for (int j=0;j<realN-1;j++){
            singleReducMat(j,j)=-constVectorComplex;
            singleReducMat(j+1,j)=complex<double>(1.0,0.0);
        }
        singleReducRhs(realN-1) = -constVectorComplex;
        
        MatrixXcd IAiA;
        if (realN>1){
            MatrixXcd invSingleReducMat = singleReducMat.completeOrthogonalDecomposition().pseudoInverse();
            IAiA = MatrixXcd::Identity(realN,realN) - singleReducMat*invSingleReducMat;
        } else IAiA = MatrixXcd::Ones(realN,realN);
        singleReducRhs = IAiA*singleReducRhs;
        for (int j=0;j<IAiA.rows();j++)
            for (int k=0;k<IAiA.cols();k++)
                alignTriplets.push_back(Triplet<complex<double>>(rowCounter+j, k*jump*pvField.intField.rows()+pvData.constSpaces(i), IAiA(j,k)));
        
        //TODO: not a great summing... also using random access which is quite bad...
        alignRhsList.push_back(singleReducRhs);
        for (int j=0;j<singleReducRhs.size();j++){
            WAlignTriplets.push_back(Triplet<complex<double>>(rowCounter+j, rowCounter+j, pvData.wAlignment(i)*pvField.tb->tangentSpaceMass.coeff(pvData.constSpaces(i),pvData.constSpaces(i))));
            pvData.totalConstrainedWeight+=pvField.tb->tangentSpaceMass.coeff(pvData.constSpaces(i),pvData.constSpaces(i));
        }
        rowCounter+=realN;
    }
    
    if (noSoftAlignment)
        pvData.totalConstrainedWeight=1.0;  //it wouldn't be used, except just to avoid a division by zero in the energy formulation
    
    pvData.alignRhs.resize(rowCounter);
    for (int i=0;i<alignRhsList.size();i++)
        pvData.alignRhs.segment(i*realN,realN)=alignRhsList[i];
    
    pvData.alignMat.resize(rowCounter, pvData.N*pvField.intField.rows());
    pvData.alignMat.setFromTriplets(alignTriplets.begin(), alignTriplets.end());
    
    pvData.WAlign.resize(rowCounter,rowCounter);
    pvData.WAlign.setFromTriplets(WAlignTriplets.begin(), WAlignTriplets.end());
}


// Computes a polyvector field on the entire mesh, where precomputation has taken place.
// Inputs:
//  PolyVectorData: The data structure which should have been initialized with polyvector_precompute()
// Outputs:
//  pvField: a POLYVECTOR_FIELD type cartesian field object
inline void polyvector_field(PolyVectorData& pvData,
                             directional::CartesianField& pvField)
{
    using namespace std;
    using namespace Eigen;
    
    //Using a temporary pvData so it could be updated and given to the iteration functions if needed
    polyvector_precompute(pvField,pvData);
    
    //forming total energy matrix;
    SparseMatrix<complex<double>> totalUnreducedLhs =(pvData.smoothMat.adjoint()*pvData.WSmooth*pvData.smoothMat) * (pvData.wSmooth / pvData.totalSmoothWeight);
    if (pvData.roSyMat.rows()!=0)
        totalUnreducedLhs=totalUnreducedLhs+(pvData.wRoSy*pvData.roSyMat.adjoint()*pvData.WRoSy*pvData.roSyMat)/pvData.totalRoSyWeight;
    if (pvData.alignMat.rows()!=0)
        totalUnreducedLhs=totalUnreducedLhs+(pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignMat)/pvData.totalConstrainedWeight;
    VectorXcd totalUnreducedRhs= (pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignRhs)/pvData.totalConstrainedWeight;
    
    pvData.totalLhs = pvData.reducMat.adjoint()*totalUnreducedLhs*pvData.reducMat;
    pvData.totalRhs = pvData.reducMat.adjoint()*(totalUnreducedRhs - totalUnreducedLhs*pvData.reducRhs);
    
    //Initial solution (or only solution if numIterations = 0)
    SimplicialLDLT<SparseMatrix<complex<double>>> solver;
    solver.compute(pvData.totalLhs);
    pvData.reducedDofs = solver.solve(pvData.totalRhs);
    assert(solver.info() == Success && "PolyVector solver failed!");
    VectorXcd fullDofs = pvData.reducMat*pvData.reducedDofs+pvData.reducRhs;
    
    MatrixXcd intField(pvData.tb->numSpaces, pvData.N);
    for (int i=0;i<pvData.N;i++)
        intField.col(i) = fullDofs.segment(i*pvData.tb->numSpaces,pvData.tb->numSpaces);
    
    pvField.fieldType = fieldTypeEnum::POLYVECTOR_FIELD;
    pvField.set_intrinsic_field(intField);
    
    double totalMass;
    if (pvData.iterationMode){
        if (pvData.verbose)
            std::cout<<"Iteration Mode"<<std::endl;
        //approximating the smallest non-zero eigenvalues
        complex<double> energy = (0.5 * pvData.reducedDofs.adjoint() * pvData.totalLhs* pvData.reducedDofs- pvData.reducedDofs.adjoint() * pvData.totalRhs).coeff(0,0);
        complex<double> mass = (pvData.reducedDofs.adjoint() * pvData.reducMat.adjoint()*pvData.M*pvData.reducMat * pvData.reducedDofs).coeff(0,0);
        double approxEig = std::abs(energy/mass);
        pvData.currImplicitCoeff = pvData.initImplicitFactor/approxEig;
        pvData.currIteration = 0;
        
        if (pvData.verbose)
            cout<<"initial implicit coefficient: "<<pvData.currImplicitCoeff<<endl;
        pvData.reducProjSolver.compute(pvData.reducMat.adjoint()*pvData.M*pvData.reducMat);
        assert(pvData.reducProjSolver.info() == Success && "Reduction Projection solver failed!");
        
        pvData.implicitLhs = pvData.reducMat.adjoint()*pvData.M*pvData.reducMat +  pvData.currImplicitCoeff*pvData.totalLhs;
        pvData.implicitSolver.compute(pvData.implicitLhs);
        assert(pvData.implicitSolver.info() == Success && "Implicit factorization failed!");
    }
}

// Computes iterations of the extended PolyVector algorithm, which includes a single implicit step and running the set of projection iteration functions by order once.
// Inputs:
//  pvData: The data structure which should have been initialized with polyvector_precompute()
//  pvField: The current POLYVECTOR_FIELD type cartesian field object
//  iterationFunctions: the iteration functions, as objects of type PvIterationFunction (see examples in polyvector_iteration_functions.h)
//  numIterations: how many iterations (of everyting) to run
// Outputs:
//  pvField: a POLYVECTOR_FIELD type cartesian field object
inline void polyvector_iterate(PolyVectorData& pvData,
                               directional::CartesianField& pvField,
                               const std::vector<directional::PvIterationFunction> iterationFunctions,
                               const int numIterations=1){
    
    assert(pvData.iterationMode && "polyvector_iterate(): Iteration mode has not been set");
    for (int i=0;i<numIterations;i++){
        if (pvData.verbose)
            std::cout<<"Iteration no. "<<pvData.currIteration<<std::endl;
        
        //An implicit step to reduce the energy
        if (pvData.verbose)
            std::cout<<"Energy before implicit step: "<<(pvData.totalLhs*pvData.reducedDofs-pvData.totalRhs).cwiseAbs().maxCoeff()<<std::endl;
        pvData.implicitRhs = pvData.currImplicitCoeff*pvData.totalRhs + pvData.reducMat.adjoint()*pvData.M*pvData.reducMat * pvData.reducedDofs;
        pvData.reducedDofs = pvData.implicitSolver.solve(pvData.implicitRhs);
        if (pvData.verbose)
            std::cout<<"Energy after implicit step: "<<(pvData.totalLhs*pvData.reducedDofs-pvData.totalRhs).cwiseAbs().maxCoeff()<<std::endl;
        
        Eigen::VectorXcd fullDofs = pvData.reducMat*pvData.reducedDofs+pvData.reducRhs;
        
        Eigen::MatrixXcd intField(pvData.tb->numSpaces, pvData.N);
        for (int i=0;i<pvData.N;i++)
            intField.col(i) = fullDofs.segment(i*pvData.tb->numSpaces,pvData.tb->numSpaces);
        
        pvField.set_intrinsic_field(intField);
        
        //running the iteration over the prescribed functions
        for (int i=0;i<iterationFunctions.size();i++)
            pvField = iterationFunctions[i](pvField, pvData);
        
        pvData.currImplicitCoeff*=pvData.implicitScheduler;
        if (std::abs(pvData.implicitScheduler-1.0)>10e-6){
            if (pvData.verbose)
                std::cout<<"Current implicitCoeff: "<<pvData.currImplicitCoeff<<std::endl;
            pvData.implicitLhs = pvData.reducMat.adjoint()*pvData.M*pvData.reducMat +  pvData.currImplicitCoeff*pvData.totalLhs;
            pvData.implicitSolver.compute(pvData.implicitLhs);
            assert(pvData.implicitSolver.info() == Eigen::Success && "Implicit factorization failed!");
        }
        pvData.currIteration++;
        
        //recreating prevSolution with reducedDof
        for (int i = 0; i < pvField.intField.cols() / 2; ++i) {
            fullDofs.segment(pvField.intField.rows() * i, pvField.intField.rows()).real() = pvField.intField.col(2 * i);
            fullDofs.segment(pvField.intField.rows() * i, pvField.intField.rows()).imag() = pvField.intField.col(2 * i + 1);
        }
        pvData.reducedDofs = pvData.reducProjSolver.solve(pvData.reducMat.adjoint()*pvData.M*(fullDofs-pvData.reducRhs));
    }
}

}

#endif
