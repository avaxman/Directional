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
#include <igl/doublearea.h>
#include <igl/speye.h>
#include <igl/eigs.h>
#include <iostream>
#include <directional/circumcircle.h>
#include <directional/complex_eigs.h>
#include <directional/CartesianField.h>

namespace directional
{

#define UNIFORM_WEIGHTS 0
#define BARYCENTRIC_WEIGHTS 1
#define INV_COT_WEIGHTS 2

  struct PolyVectorData{
  public:
    
    //User parameters
    Eigen::VectorXi constSpaces;   //list of tangent spaces where there are (partial) constraints. The faces can repeat to constrain more vectors
    Eigen::MatrixXd constVectors; //corresponding to constSpaces.
    
    int N;                        //Degree of field
    int sizeT;                    //#tangent spaces
    bool signSymmetry;            //Whether field enforces a ssign symmetry (only when N is even, otherwise by default set to false)
    bool perfectRoSy;             //Whether the field must be perfect rotationally-symmetric (but not unit).
    double wSmooth;               //Weight of smoothness
    double wRoSy;                 //Weight of rotational-symmetry. "-1" means a perfect RoSy field (power field)
    Eigen::VectorXd wAlignment;   //Weight of alignment per each of the constfaces. "-1" means a fixed vector
    
    int lapType;                  //Choice of weights (from UNIFORM_WEIGHTS,BARYCENTRIC_WEIGHTS or INV_COT_WEIGHTS)
    
    Eigen::SparseMatrix<std::complex<double>> smoothMat;    //Smoothness energy
    Eigen::SparseMatrix<std::complex<double>> roSyMat;      //Rotational-symmetry energy
    Eigen::SparseMatrix<std::complex<double>> alignMat;     //(soft) alignment energy.
    Eigen::SparseMatrix<std::complex<double>> reducMat;     //reducing the fixed dofs (for instance with sign symmetry or fixed partial constraints)
    Eigen::VectorXcd reducRhs;                                   //The uncompressed PV coeffs are reducMat*true_dofs+reducRhs
    Eigen::VectorXcd alignRhs;                                   //encoding the soft constraints
    
    //Mass and stiffness matrices
    Eigen::SparseMatrix<std::complex<double>> WSmooth, WAlign, WRoSy, M;
    double totalRoSyWeight, totalConstrainedWeight, totalSmoothWeight;    //for co-scaling energies
    
    PolyVectorData():signSymmetry(true), lapType(BARYCENTRIC_WEIGHTS), wSmooth(1.0), wRoSy(0.0) {wAlignment.resize(0); constSpaces.resize(0); constVectors.resize(0,3);}
    ~PolyVectorData(){}
  };
  
  
  // Precalculate the operators according to the user-prescribed parameters. Must be called whenever any of them changes
  // Inputs:
  //  V:      #V by 3 vertex coordinates.
  //  F:      #F by 3 face vertex indices.
  //  EV:     #E by 2 matrix of edges (vertex indices)
  //  EF:     #E by 2 matrix of oriented adjacent faces
  //  B1, B2: #F by 3 matrices representing the local base of each face.
  //  PolyVectorData (must fill non-default values in advance)
  // Outputs:
  //  PolyVectorData:       Updated structure with all operators
  IGL_INLINE void polyvector_precompute(const directional::CartesianField& pvField,
                                        const int N,
                                        PolyVectorData& pvData)
  {
    
    using namespace std;
    using namespace Eigen;
    
    //assert(pvField.connection.cols()==N && "The PolyVector field should be initialized to be of order N");
    
    //Building the smoothness matrices, with an energy term for each inner edge and degree
    int rowCounter=0;
    std::vector< Triplet<complex<double> > > dTriplets, WTriplets;
    pvData.N = N;
    pvData.sizeT = pvField.intField.rows();
    if (pvData.N%2!=0) pvData.signSymmetry=false;  //it has to be for odd N
    
   
    
    pvData.totalSmoothWeight = pvField.stiffnessWeights.sum();
    
    vector<Triplet<complex<double>>> WSmoothTriplets, MTriplets;
    //VectorXd doubleAreas;
    //igl::doublearea(mesh.V,mesh.F,doubleAreas);
    for (int n = 0; n < pvData.N; n++)
    {
      for (int i=0;i<pvField.adjSpaces.rows();i++)
      {
        if ((pvField.adjSpaces(i,0)==-1)||(pvField.adjSpaces(i,1)==-1))
          continue;  //boundary edge
        
        // Compute the complex representation of the common edge
        /*RowVector3d e = mesh.V.row(mesh.EV(i,1)) - mesh.V.row(mesh.EV(i,0));
        RowVector2d vef = Vector2d(e.dot(mesh.Bx.row(mesh.EF(i,0))), e.dot(mesh.By.row(mesh.EF(i,0)))).normalized();
        complex<double> ef(vef(0), vef(1));
        Vector2d veg = Vector2d(e.dot(mesh.Bx.row(mesh.EF(i,1))), e.dot(mesh.By.row(mesh.EF(i,1)))).normalized();
        complex<double> eg(veg(0), veg(1));*/
        
        // differential matrix between two tangent spaces
        dTriplets.push_back(Triplet<complex<double> >(rowCounter, n*pvField.intField.rows()+pvField.adjSpaces(i,0), pow(pvField.connection(i),pvData.N-n)));
        dTriplets.push_back(Triplet<complex<double> >(rowCounter, n*pvField.intField.rows()+pvField.adjSpaces(i,1), -1.0));
        
        //stiffness weights
        WSmoothTriplets.push_back(Triplet<complex<double> >(rowCounter, rowCounter, pvField.stiffnessWeights(i)));
        rowCounter++;
      }
      
      for (int i=0;i<pvField.intField.rows();i++)
        MTriplets.push_back(Triplet<complex<double>>(n*pvField.intField.rows()+i, n*pvField.intField.rows()+i, pvField.massWeights(i)));
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
    std::vector<MatrixXcd> localFaceReducMats; localFaceReducMats.resize(pvField.intField.rows());
    std::vector<VectorXcd> localSpaceReducRhs;  localSpaceReducRhs.resize(pvField.intField.rows());
    
    for (int i=0;i<pvField.intField.rows();i++){
      localFaceReducMats[i]=MatrixXcd::Identity(realN,realN);
      localSpaceReducRhs[i]=VectorXcd::Zero(realN);
    }
    
    /*************Hard-constraint reduction matrices******************/
    VectorXcd constVectorsIntrinsic=pvField.project_to_intrinsic(pvData.constSpaces,pvData.constVectors);
    for (int i=0;i<pvData.constSpaces.size();i++){
      if (pvData.wAlignment(i)>=0.0)
        continue;  //here we only handle the reduction caused by a fixed dof
      if (numSpaceConstraints(pvData.constSpaces(i))==realN)
        continue; //overconstrained; we ignore any further constraints on that face
      
      int currSpaceNumDof = numSpaceConstraints(pvData.constSpaces(i));
      
      
      //std::complex<double>(pvData.constVectors.row(i).dot(mesh.Bx.row(pvData.constSpaces(i))), pvData.constVectors.row(i).dot(mesh.By.row(pvData.constSpaces(i))));
      complex<double> constVectorComplex = (pvData.signSymmetry ? constVectorsIntrinsic(i)*constVectorsIntrinsic(i) : constVectorsIntrinsic(i));
      constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorsIntrinsic(i), N) : constVectorsIntrinsic(i));
      
      MatrixXcd singleReducMat=MatrixXcd::Zero(realN-currSpaceNumDof,realN-currSpaceNumDof-1);
      VectorXcd singleReducRhs=VectorXcd::Zero(realN-currSpaceNumDof);
      for (int j=0;j<realN-currSpaceNumDof-1;j++){
        singleReducMat(j,j)=-constVectorComplex;
        singleReducMat(j+1,j)=complex<double>(1.0,0.0);
      }
      singleReducRhs(realN-currSpaceNumDof-1) = -constVectorComplex;
      
      localSpaceReducRhs[pvData.constSpaces(i)] = localFaceReducMats[pvData.constSpaces(i)]*singleReducRhs + localSpaceReducRhs[pvData.constSpaces(i)];
      localFaceReducMats[pvData.constSpaces(i)] = localFaceReducMats[pvData.constSpaces(i)]*singleReducMat;
      
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
      //std::cout<<"localFaceReducMats[i]: "<<localFaceReducMats[i]<<std::endl;
      for (int j=0;j<pvData.N;j+=jump){
        for (int k=0;k<localFaceReducMats[i].cols();k++)
          reducMatTriplets.push_back(Triplet<complex<double>>(j*pvField.intField.rows()+i, colCounter+k, localFaceReducMats[i](j/jump,k)));
        
        pvData.reducRhs(j*pvField.intField.rows()+i) = localSpaceReducRhs[i](j/jump);
      }
     
      colCounter+=localFaceReducMats[i].cols();
    }
    
    pvData.reducMat.resize(pvData.N*pvField.intField.rows(), colCounter);
    pvData.reducMat.setFromTriplets(reducMatTriplets.begin(), reducMatTriplets.end());
  
    
    /****************rotational-symmetry matrices********************/
    
    if (pvData.wRoSy >= 0.0){ //this is anyhow enforced, this matrix is unnecessary)
      vector<Triplet<complex<double>>> roSyTriplets, WRoSyTriplets;
      for (int i=pvField.intField.rows();i<pvData.N*pvField.intField.rows();i++){
        roSyTriplets.push_back(Triplet<complex<double>>(i,i,1.0));
        WRoSyTriplets.push_back(Triplet<complex<double>>(i,i,pvField.massWeights(i%pvField.intField.rows())));
      }
      
      pvData.roSyMat.resize(N*pvField.intField.rows(), N*pvField.intField.rows());
      pvData.roSyMat.setFromTriplets(roSyTriplets.begin(), roSyTriplets.end());
      
      pvData.WRoSy.resize(N*pvField.intField.rows(), N*pvField.intField.rows());
      pvData.WRoSy.setFromTriplets(WRoSyTriplets.begin(), WRoSyTriplets.end());
      
      pvData.totalRoSyWeight=((double)pvData.N)*pvField.massWeights.sum();
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
      complex<double> constVectorComplex = (pvData.signSymmetry ? constVectorsIntrinsic(i)*constVectorsIntrinsic(i) : constVectorsIntrinsic(i));
      constVectorComplex = (pvData.wRoSy < 0.0 ? pow(constVectorsIntrinsic(i), N) : constVectorsIntrinsic(i));
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
        WAlignTriplets.push_back(Triplet<complex<double>>(rowCounter+j, rowCounter+j, pvData.wAlignment(i)*pvField.massWeights(pvData.constSpaces(i))));
        pvData.totalConstrainedWeight+=pvField.massWeights(pvData.constSpaces(i));
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


  // Computes a polyvector on the entire mesh
  // Inputs:
  //  PolyVectorData: The data structure which should have been initialized with polyvector_precompute()
  // Outputs:
  //  polyVectorField: #F by N The output interpolated field, in polyvector (complex polynomial) format.
  IGL_INLINE void polyvector_field(const PolyVectorData& pvData,
                                   directional::CartesianField& pvField)
  {
    using namespace std;
    using namespace Eigen;
    
    //forming total energy matrix;
    SparseMatrix<complex<double>> totalUnreducedLhs = pvData.wSmooth * (pvData.smoothMat.adjoint()*pvData.WSmooth*pvData.smoothMat)/pvData.totalSmoothWeight + (pvData.wRoSy*pvData.roSyMat.adjoint()*pvData.WRoSy*pvData.roSyMat)/pvData.totalRoSyWeight + (pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignMat)/pvData.totalConstrainedWeight;
    VectorXcd totalUnreducedRhs= (pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignRhs)/pvData.totalConstrainedWeight;
        
    SparseMatrix<complex<double>> totalLhs = pvData.reducMat.adjoint()*totalUnreducedLhs*pvData.reducMat;
    VectorXcd totalRhs = pvData.reducMat.adjoint()*(totalUnreducedRhs - totalUnreducedLhs*pvData.reducRhs);
    if (pvData.constSpaces.size() == 0)  //alignmat should be empty and the reduction matrix should be only sign symmetry, if applicable
    {
     //using a matrix with only the first sizeT x sizeT block
      vector<Triplet<complex<double>>> X0LhsTriplets, X0MTriplets;
      SparseMatrix<complex<double>> X0Lhs, X0M;
      for (int k=0; k<totalUnreducedLhs.outerSize(); ++k)
        for (SparseMatrix<std::complex<double>>::InnerIterator it(totalUnreducedLhs,k); it; ++it)
          if ((it.row()<pvData.sizeT)&&(it.col()<pvData.sizeT))
            X0LhsTriplets.push_back(Triplet<complex<double>>(it.row(), it.col(), it.value()*pvData.totalSmoothWeight));  //to bypass the early convergence of igl eigenvalues...
      
      X0Lhs.resize(pvData.sizeT, pvData.sizeT);
      X0Lhs.setFromTriplets(X0LhsTriplets.begin(), X0LhsTriplets.end());
      
      for (int k=0; k<pvData.M.outerSize(); ++k)
        for (SparseMatrix<std::complex<double>>::InnerIterator it(pvData.M,k); it; ++it)
          if ((it.row()<pvData.sizeT)&&(it.col()<pvData.sizeT))
            X0MTriplets.push_back(Triplet<complex<double>>(it.row(), it.col(), it.value()*pvData.totalSmoothWeight));
      
      X0M.resize(pvData.sizeT, pvData.sizeT);
      X0M.setFromTriplets(X0MTriplets.begin(), X0MTriplets.end());
          
      //Extracting first eigenvector
      Eigen::MatrixXcd U;
      Eigen::VectorXcd S;
      complex_eigs(X0Lhs, X0M, 10, U, S);
      int smallestIndex; S.cwiseAbs().minCoeff(&smallestIndex);
      
      pvField.fieldType = POLYVECTOR_FIELD;
      MatrixXcd intField(U.col(0).rows(),pvData.N);
      intField.col(0)=U.col(smallestIndex);
      pvField.set_intrinsic_field(intField);
    } else { //just solving the system
      SimplicialLDLT<SparseMatrix<complex<double>>> solver;
      //solver.analyzePattern(totalLhs);   // for this step the numerical values of A are not used
      solver.compute(totalLhs);
      VectorXcd reducedDofs = solver.solve(totalRhs);
      assert(solver.info() == Success);
      VectorXcd fullDofs = pvData.reducMat*reducedDofs+pvData.reducRhs;
      MatrixXcd intField;
      for (int i=0;i<pvData.N;i++)
        intField.col(i) = fullDofs.segment(i*pvData.sizeT,pvData.sizeT);
      
      pvField.fieldType = POLYVECTOR_FIELD;
      pvField.set_intrinsic_field(intField);
      
      //std::cout<<"Smoothness energy: "<<pvData.wSmooth * (fullDofs.adjoint()*pvData.smoothMat.adjoint()*pvData.WSmooth*pvData.smoothMat*fullDofs)/pvData.totalSmoothWeight<<std::endl;
      //std::cout<<"RoSy Energy: "<<fullDofs.adjoint()* ((pvData.wRoSy*pvData.roSyMat.adjoint()*pvData.WRoSy*pvData.roSyMat)/pvData.totalRoSyWeight)*fullDofs<<std::endl;
      //TODO: measure energy in the correct way
      /*std::cout<<"Alignment energy: "<< fullDofs.adjoint()*((pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignMat*fullDofs - pvData.alignMat.adjoint()*pvData.WAlign*pvData.alignRhs)/pvData.totalConstrainedWeight)<<std::endl;*/
      
    }
    

  }

  
  // minimal version without auxiliary data
  IGL_INLINE void polyvector_field(directional::CartesianField& pvField,
                                   const Eigen::VectorXi& constSpaces,
                                   const Eigen::MatrixXd& constVectors,
                                   const double smoothWeight,
                                   const double roSyWeight,
                                   const Eigen::VectorXd& alignWeights,
                                   const int N)
  {
    PolyVectorData pvData;
    pvData.constSpaces=constSpaces;
    pvData.constVectors=constVectors;
    pvData.wAlignment = alignWeights;
    pvData.wSmooth = smoothWeight;
    pvData.wRoSy = roSyWeight;
    polyvector_precompute(pvField, N, pvData);
    polyvector_field(pvData, pvField);
  }

IGL_INLINE void polyvector_field(directional::CartesianField& pvField,
                                 const Eigen::VectorXi& constSpaces,
                                 const Eigen::MatrixXd& constVectors,
                                 const int N)
{
  
  PolyVectorData pvData;
  pvData.constSpaces=constSpaces;
  pvData.constVectors=constVectors;
  pvData.wAlignment = Eigen::VectorXd::Constant(constSpaces.size(),-1.0);
  pvData.wSmooth = 1.0;
  pvData.wRoSy = 0.0;
  polyvector_precompute(pvField, N, pvData);
  polyvector_field(pvData, pvField);
}
}

#endif
