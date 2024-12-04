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
#include <directional/project_curl.h>
#include <directional/principal_matching.h>

namespace directional
{

    //Data for the precomputation of the PolyVector algorithm
    struct PolyVectorData{
    public:

        //User parameters
        Eigen::VectorXi constSpaces;    // List of tangent spaces where there are (partial) constraints. The faces can repeat to constrain more vectors
        Eigen::MatrixXd constVectors;   // Corresponding to constSpaces.

        int N;                          // Degree of field
        int sizeT;                      // #tangent spaces
        bool signSymmetry;              // Whether field enforces a ssign symmetry (only when N is even, otherwise by default set to false)
        bool perfectRoSy;               // Whether the field must be perfect rotationally-symmetric (but not unit).
        double wSmooth;                 // Weight of smoothness
        double wRoSy;                   // Weight of rotational-symmetry. "-1" means a perfect RoSy field (power field)
        Eigen::VectorXd wAlignment;     // Weight of alignment per each of the constfaces. "-1" means a fixed vector
        bool projectCurl;               // Project out the curl of the field
        bool normalizeField;            // Normalize the field (per vector)
        int numIterations;              //  Iterate energy reduction->(possibly)normalize->(possibly)project curl

        Eigen::SparseMatrix<std::complex<double>> smoothMat;    //Smoothness energy
        Eigen::SparseMatrix<std::complex<double>> roSyMat;      //Rotational-symmetry energy
        Eigen::SparseMatrix<std::complex<double>> alignMat;     //(soft) alignment energy.
        Eigen::SparseMatrix<std::complex<double>> reducMat;     //reducing the fixed dofs (for instance with sign symmetry or fixed partial constraints)
        Eigen::VectorXcd reducRhs;                              //The uncompressed PV coeffs are reducMat*true_dofs+reducRhs
        Eigen::VectorXcd alignRhs;                              //encoding the soft constraints

        //Mass and stiffness matrices
        Eigen::SparseMatrix<std::complex<double>> WSmooth, WAlign, WRoSy, M;
        double totalRoSyWeight, totalConstrainedWeight, totalSmoothWeight;    //for co-scaling energies

        PolyVectorData():signSymmetry(true),  wSmooth(1.0), wRoSy(0.0), numIterations(0), projectCurl(false), normalizeField(false) {wAlignment.resize(0); constSpaces.resize(0); constVectors.resize(0,3);}
        ~PolyVectorData(){}
    };


    // Precalculate the operators needed for PolyVector computation according to the user-prescribed parameters. Must be called whenever any of them changes
    // Input:
    //  tb:     underlying tangent bundle
    //  N:      degree of the field
    //
    // Output:
    //  pvField: POLYVECTOR_FIELD cartesian field initalized with the tangent bundle
    //  pvData:  Updated structure with all operators
    inline void polyvector_precompute(const directional::TangentBundle& tb,
                                          const int N,
                                          directional::CartesianField& pvField,
                                          PolyVectorData& pvData)
    {

        using namespace std;
        using namespace Eigen;

        pvField.init(tb, fieldTypeEnum::POLYVECTOR_FIELD, N);

        assert(pvData.projectCurl && tb.discTangType()==directional::FACE_SPACES && "Projecting curl only works for face-based fields for now!");

        //Building the smoothness matrices, with an energy term for each inner edge and degree
        int rowCounter=0;
        std::vector< Triplet<complex<double> > > dTriplets, WTriplets;
        pvData.N = N;
        pvData.sizeT = pvField.intField.rows();
        if (pvData.N%2!=0) pvData.signSymmetry=false;  //it has to be for odd N

        pvData.totalSmoothWeight = pvField.tb->connectionMass.sum();

        vector<Triplet<complex<double>>> WSmoothTriplets, MTriplets;
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
                WSmoothTriplets.push_back(Triplet<complex<double> >(rowCounter, rowCounter, pvField.tb->connectionMass(i)));
                rowCounter++;
            }

            for (int i=0;i<pvField.intField.rows();i++)
                MTriplets.push_back(Triplet<complex<double>>(n*pvField.intField.rows()+i, n*pvField.intField.rows()+i, pvField.tb->tangentSpaceMass(i)));
        }

        pvData.smoothMat.resize(rowCounter, pvData.N*pvField.intField.rows());
        pvData.smoothMat.setFromTriplets(dTriplets.begin(), dTriplets.end());

        pvData.WSmooth.resize(rowCounter, rowCounter);
        pvData.WSmooth.setFromTriplets(WSmoothTriplets.begin(), WSmoothTriplets.end());

        pvData.M.resize(pvData.N*pvField.intField.rows(), pvData.N*pvField.intField.rows());
        pvData.M.setFromTriplets(MTriplets.begin(), MTriplets.end());

        //creating reduction transformation
        VectorXi numSpaceConstraints = Eigen::VectorXi::Zero(pvField.intField.rows());
        int realN = (pvData.signSymmetry ? N/2 : N);
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
            constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorComplexRaw, N) :constVectorComplex);

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
            vector<Triplet<complex<double>>> roSyTriplets, WRoSyTriplets;
            for (int i=pvField.intField.rows();i<pvData.N*pvField.intField.rows();i++){
                roSyTriplets.push_back(Triplet<complex<double>>(i,i,1.0));
                WRoSyTriplets.push_back(Triplet<complex<double>>(i,i,pvField.tb->tangentSpaceMass(i%pvField.intField.rows())));
            }

            pvData.roSyMat.resize(N*pvField.intField.rows(), N*pvField.intField.rows());
            pvData.roSyMat.setFromTriplets(roSyTriplets.begin(), roSyTriplets.end());

            pvData.WRoSy.resize(N*pvField.intField.rows(), N*pvField.intField.rows());
            pvData.WRoSy.setFromTriplets(WRoSyTriplets.begin(), WRoSyTriplets.end());

            pvData.totalRoSyWeight=((double)pvData.N)*pvField.tb->tangentSpaceMass.sum();
        } else {
            pvData.roSyMat.resize(0, N*pvField.intField.rows());
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
            constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorComplexRaw, N) : constVectorComplex);
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

            alignRhsList.push_back(singleReducRhs);
            for (int j=0;j<singleReducRhs.size();j++){
                WAlignTriplets.push_back(Triplet<complex<double>>(rowCounter+j, rowCounter+j, pvData.wAlignment(i)*pvField.tb->tangentSpaceMass(pvData.constSpaces(i))));
                pvData.totalConstrainedWeight+=pvField.tb->tangentSpaceMass(pvData.constSpaces(i));
            }
            rowCounter+=realN;
        }

        if (noSoftAlignment)
            pvData.totalConstrainedWeight=1.0;  //it wouldn't be used, except just to avoid a division by zero in the energy formulation

        pvData.alignRhs.resize(rowCounter);
        for (int i=0;i<alignRhsList.size();i++)
            pvData.alignRhs.segment(i*realN,realN)=alignRhsList[i];

        pvData.alignMat.resize(rowCounter, N*pvField.intField.rows());
        pvData.alignMat.setFromTriplets(alignTriplets.begin(), alignTriplets.end());

        pvData.WAlign.resize(rowCounter,rowCounter);
        pvData.WAlign.setFromTriplets(WAlignTriplets.begin(), WAlignTriplets.end());
    }


    // Computes a polyvector field on the entire mesh, where precomputation has taken place.
    // Inputs:
    //  PolyVectorData: The data structure which should have been initialized with polyvector_precompute()
    // Outputs:
    //  pvField: a POLYVECTOR_FIELD type cartesian field object

    inline void polyvector_field(const PolyVectorData& pvData,
                                 directional::CartesianField& pvField)
    {
        using namespace std;
        using namespace Eigen;

        //forming total energy matrix;
        SparseMatrix<complex<double>> totalUnreducedLhs =(pvData.smoothMat.adjoint()*pvData.WSmooth*pvData.smoothMat) * (pvData.wSmooth / pvData.totalSmoothWeight);
        if (pvData.roSyMat.rows()!=0)
            totalUnreducedLhs=totalUnreducedLhs+(pvData.wRoSy*pvData.roSyMat.adjoint()*pvData.WRoSy*pvData.roSyMat)/pvData.totalRoSyWeight;
        if (pvData.alignMat.rows()!=0)
            totalUnreducedLhs=totalUnreducedLhs+(pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignMat)/pvData.totalConstrainedWeight;
        VectorXcd totalUnreducedRhs= (pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignRhs)/pvData.totalConstrainedWeight;

        SparseMatrix<complex<double>> totalLhs = pvData.reducMat.adjoint()*totalUnreducedLhs*pvData.reducMat;
        VectorXcd totalRhs = pvData.reducMat.adjoint()*(totalUnreducedRhs - totalUnreducedLhs*pvData.reducRhs);

        //Initial solution (or only solution if numIterations = 0)
        SimplicialLDLT<SparseMatrix<complex<double>>> solver;
        solver.compute(totalLhs);
        VectorXcd reducedDofs = solver.solve(totalRhs);
        assert(solver.info() == Success & "PolyVector solver failed!");
        VectorXcd fullDofs = pvData.reducMat*reducedDofs+pvData.reducRhs;
        MatrixXcd intField(pvData.sizeT, pvData.N);
        for (int i=0;i<pvData.N;i++)
            intField.col(i) = fullDofs.segment(i*pvData.sizeT,pvData.sizeT);

        pvField.fieldType = fieldTypeEnum::POLYVECTOR_FIELD;
        pvField.set_intrinsic_field(intField);

        if (pvData.normalizeField){
           CartesianField rawField;
           polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
           raw_to_polyvector(rawField,  pvField);
        }

        //testing raw_to_polyvector
        CartesianField rawField, pvField2;
        polyvector_to_raw(pvField, rawField);
        directional::raw_to_polyvector(rawField, pvField2);

        std::cout<<"raw_to_polyvector(polyvector_to_raw()) test: "<<(pvField.intField-pvField2.intField).cwiseAbs().maxCoeff()<<std::endl;

        //Doing reduce energy-renormalize-project curl iterations
        for (int i=0;i<pvData.numIterations;i++){

            //TODO: iteration of implicit Euler step



            if (pvData.normalizeField){
                CartesianField rawField;
                polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
                directional::raw_to_polyvector(rawField, pvField);
            }

            if (pvData.projectCurl){
                CartesianField rawField, curlFreeField;
                polyvector_to_raw(pvField, rawField, pvData.N%2==0, false);
                directional::principal_matching(rawField);
                project_curl(rawField, Eigen::VectorXi(), Eigen::MatrixXd(), curlFreeField);
                directional::raw_to_polyvector(rawField,  pvField);


        }



    }


    // minimal version without auxiliary data
    inline void polyvector_field(const TangentBundle& tb,
                                 const Eigen::VectorXi& constSpaces,
                                 const Eigen::MatrixXd& constVectors,
                                 const double smoothWeight,
                                 const double roSyWeight,
                                 const Eigen::VectorXd& alignWeights,
                                 const int N,
                                 directional::CartesianField& pvField,
                                 const bool projectCurl = false,
                                 const bool normalizeField = false,
                                 const int numIterations = false)
    {
        PolyVectorData pvData;
        if (constSpaces.size()!=0) {
            pvData.constSpaces = constSpaces;
            pvData.constVectors = constVectors;
            pvData.wAlignment = alignWeights;
            pvData.projectCurl = projectCurl;
            pvData.normalizeField = normalizeField;
            pvData.numIterations = numIterations;
        }else{
            pvData.constSpaces.resize(1); pvData.constSpaces(0)=0;
            Eigen::RowVector2d intConstVector; intConstVector<<1.0,0.0;
            pvData.constVectors = tb.project_to_extrinsic(pvData.constSpaces, intConstVector);
            pvData.wAlignment = Eigen::VectorXd::Constant(pvData.constSpaces.size(),-1.0);
        }
        pvData.wSmooth = smoothWeight;
        pvData.wRoSy = roSyWeight;
        pvField.init(tb,fieldTypeEnum::POLYVECTOR_FIELD,N);
        polyvector_precompute(tb,N,pvField,pvData);
        polyvector_field(pvData, pvField);
    }


//A version with default parameters (in which alignment is hard by default).
    inline void polyvector_field(const TangentBundle& tb,
                                 const Eigen::VectorXi& constSpaces,
                                 const Eigen::MatrixXd& constVectors,
                                 const int N,
                                 directional::CartesianField& pvField)
    {

        PolyVectorData pvData;
        //in case the const spaces is zero, using a single const space which is default, just to offset the smoothest field rotation null-space
        if (constSpaces.size()!=0) {
            pvData.constSpaces = constSpaces;
            pvData.constVectors = constVectors;
        }else{
            pvData.constSpaces.resize(1); pvData.constSpaces(0)=0;
            Eigen::RowVector2d intConstVector; intConstVector<<1.0,0.0;
            pvData.constVectors = tb.project_to_extrinsic(pvData.constSpaces, intConstVector);
            pvData.wAlignment = Eigen::VectorXd::Constant(pvData.constSpaces.size(),-1.0);
        }
        pvData.wAlignment = Eigen::VectorXd::Constant(constSpaces.size(),-1.0);
        pvData.wSmooth = 1.0;
        pvData.wRoSy = 0.0;
        pvData.normalizeField = pvData.projectCurl = false;
        pvData.numIterations = 0;
        polyvector_precompute(tb, N, pvField,pvData);
        polyvector_field(pvData, pvField);
    }

}

#endif
