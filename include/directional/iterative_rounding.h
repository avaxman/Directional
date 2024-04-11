// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef iterative_rounding_h
#define iterative_rounding_h

//#include <SaddlePoint/LMSolver.h>
//#include <SaddlePoint/EigenSolverWrapper.h>
//#include <SaddlePoint/check_traits.h>
//#include <SaddlePoint/DiagonalDamping.h>
//#include <directional/SIInitialSolutionTraits.h>
//#include <directional/IterativeRoundingTraits.h>
#include <iostream>
#include <Eigen/Core>
#include <iomanip>





namespace directional{


template<typename T> void printElement(T t, const int& width)
{
  std::cout << std::right << std::setw(width) << std::setfill(' ') << t;
}


bool iterative_rounding(const Eigen::SparseMatrix<double>& A,
                        const Eigen::MatrixXd& rawField,
                        const Eigen::VectorXi& fixedIndices,
                        const Eigen::VectorXd& fixedValues,
                        const Eigen::VectorXi& singularIndices,
                        const Eigen::VectorXi& integerIndices,
                        const double lengthRatio,
                        const Eigen::VectorXd& b,
                        const Eigen::SparseMatrix<double> C,
                        const Eigen::SparseMatrix<double> G,
                        const Eigen::MatrixXd& FN,
                        const int N,
                        const int n,
                        const Eigen::MatrixXd& V,
                        const Eigen::MatrixXi& F,
                        const Eigen::SparseMatrix<double>& x2CornerMat,
                        const bool fullySeamless,
                        const bool roundSeams,
                        const bool localInjectivity,
                        const bool verbose,
                        Eigen::VectorXd& fullx){
  
  using namespace Eigen;
  using namespace std;
  
  typedef SaddlePoint::EigenSolverWrapper<Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > > LinearSolver;
  
  SIInitialSolutionTraits<LinearSolver> slTraits;
  LinearSolver lSolver1,lSolver2;
  SaddlePoint::DiagonalDamping<SIInitialSolutionTraits<LinearSolver>> dISTraits(localInjectivity ? 0.01 : 0.0);
  SaddlePoint::LMSolver<LinearSolver,SIInitialSolutionTraits<LinearSolver>, SaddlePoint::DiagonalDamping<SIInitialSolutionTraits<LinearSolver> > > initialSolutionLMSolver;
  
  IterativeRoundingTraits<LinearSolver> irTraits;
  SaddlePoint::DiagonalDamping<IterativeRoundingTraits<LinearSolver>> dIRTraits(localInjectivity ? 0.01 : 0.0);
  SaddlePoint::LMSolver<LinearSolver,IterativeRoundingTraits<LinearSolver>, SaddlePoint::DiagonalDamping<IterativeRoundingTraits<LinearSolver> > > iterativeRoundingLMSolver;
  
  slTraits.A=A;
  slTraits.rawField=rawField;
  slTraits.fixedIndices=fixedIndices;
  slTraits.fixedValues=fixedValues;
  slTraits.singularIndices=singularIndices;
  slTraits.b=b;
  slTraits.C=C;
  slTraits.G=G;
  slTraits.FN=FN;
  slTraits.N=N;
  slTraits.lengthRatio=lengthRatio;
  slTraits.n=n;
  slTraits.V=V;
  slTraits.F=F;
  slTraits.x2CornerMat=x2CornerMat;
  slTraits.integerIndices=integerIndices;
  slTraits.localInjectivity=localInjectivity;
  
  //initial solution
  if (verbose)
    cout<<"Computing initial solution..."<<endl;
  slTraits.init(verbose);
  initialSolutionLMSolver.init(&lSolver1, &slTraits, &dISTraits, 100);
  //SaddlePoint::check_traits(slTraits, slTraits.initXandFieldSmall);
  initialSolutionLMSolver.solve(false);
  if (verbose){
    cout<<"Done!"<<endl;
    cout<<"Integrability error: "<<slTraits.integrability<<endl;
    cout<<"LM Solver energy: "<<initialSolutionLMSolver.energy<<endl;
    cout<<"LM 1st-order optimality: "<<initialSolutionLMSolver.fooOptimality<<endl;
  }
  
  irTraits.init(slTraits, initialSolutionLMSolver.x, roundSeams);
  
  if (!fullySeamless){
    fullx=irTraits.x0;
    return true;
  }
  
  //Iterative rounding
  int colWidth=20;
  if (verbose){
    cout<<"Starting Iterative rounding..."<<endl;
    cout << std::right << setw(colWidth) << setfill(' ') << "Index";
    cout << std::right << setw(colWidth) << setfill(' ') << "Orig. value";
    cout << std::right << setw(colWidth) << setfill(' ') << "Integer value";
    cout << std::right << setw(colWidth) << setfill(' ') << "Energy";
    cout << std::right << setw(colWidth) << setfill(' ') << "1st-ord. Optimality";
    cout << std::right << setw(colWidth) << setfill(' ') << "# Iterations";
    cout<<endl;
  }
  
  bool success=true;
  bool hasRounded=false;
  while (irTraits.leftIndices.size()!=0){
    //cout<<"i: "<<i++<<endl;
    if (!irTraits.initFixedIndices())
      continue;
    hasRounded=true;
    dIRTraits.currLambda=(localInjectivity ? 0.01 : 0.0);
    iterativeRoundingLMSolver.init(&lSolver2, &irTraits, &dIRTraits, 100, 1e-7, 1e-7);
    iterativeRoundingLMSolver.solve(false);
    if (verbose){
      printElement(irTraits.currRoundIndex, colWidth);
      printElement(irTraits.origValue, colWidth);
      printElement(irTraits.roundValue, colWidth);
      printElement(iterativeRoundingLMSolver.energy, colWidth);
      printElement(iterativeRoundingLMSolver.fooOptimality, colWidth);
      printElement(iterativeRoundingLMSolver.currIter, colWidth);
      cout<<endl;
    }
    if (!irTraits.post_checking(iterativeRoundingLMSolver.x)){
      success=false;
      if (verbose)
        cout<<"Failed to round!"<<endl;
      break;
    }
  }
  
  if (verbose)
    cout<<"Iterative rounding "<<(success ? "succeeded!" : "failed!")<<endl;
  
  if (hasRounded)
    fullx=irTraits.UFull*iterativeRoundingLMSolver.x;
  else
    fullx=irTraits.x0;  //in case nothing happens
  return success;
  
  
}
}



#endif /* iterative_rounding_h */
