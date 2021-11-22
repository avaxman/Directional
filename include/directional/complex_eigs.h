// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_COMPLEX_EIGS_H
#define DIRECTIONAL_COMPLEX_EIGS_H


namespace directional {

  IGL_INLINE void complex_eigs(const Eigen::SparseMatrix<std::complex<double>>& Q,
                               const Eigen::SparseMatrix<double> M,
                               const int numEigs,
                               Eigen::MatrixXcd& U,
                               Eigen::VectorXcd& S)
  {
    Eigen::SparseMatrix<double> QReal(2*Q.rows(),2*Q.rows());
    std::vector<Eigen::Triplet<double> > QRealTriplets, MRealTriplets;
    for (int k=0; k<Q.outerSize(); ++k)
      for (Eigen::SparseMatrix<std::complex<double>>::InnerIterator it(Q,k); it; ++it)
      {
        QRealTriplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value().real()));
        QRealTriplets.push_back(Eigen::Triplet<double>(it.row(), Q.cols()+it.col(), -it.value().imag()));
        QRealTriplets.push_back(Eigen::Triplet<double>(Q.rows()+it.row(), it.col(), it.value().imag()));
        QRealTriplets.push_back(Eigen::Triplet<double>(Q.rows()+it.row(), Q.cols()+it.col(), it.value().real()));
      }
    QReal.setFromTriplets(QRealTriplets.begin(), QRealTriplets.end());
    
    
    Eigen::SparseMatrix<double> MReal(2*M.rows(),2*M.rows());
    for (int k=0; k<M.outerSize(); ++k)
      for (Eigen::SparseMatrix<double>::InnerIterator it(M,k); it; ++it)
      {
        MRealTriplets.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
        MRealTriplets.push_back(Eigen::Triplet<double>(it.row()+M.rows(), it.col()+M.cols(), it.value()));
        
      }
    
    MReal.setFromTriplets(MRealTriplets.begin(), MRealTriplets.end());
    Eigen::MatrixXd UReal;
    Eigen::VectorXd SReal;
    
    Eigen::MatrixXd U;
    Eigen::VectorXd S;
    igl::eigs(L,M,numEigs,igl::EIGS_TYPE_SM,UReal,SReal);
    int smallestIndex; S.minCoeff(&smallestIndex);
    
    U.resize(UReal.rows()/2, N);
    S.resize(SReal.rows()/2);
    
    U = UReal.block(0,smallestIndex,UReal.rows()/2,N).cast<std::complex<double> >().array()*std::complex<double>(1,0)+
    UReal.block(UReal.rows()/2,smallestIndex,UReal.rows()/2,N).cast<std::complex<double> >().array()*std::complex<double>(0,1);
    
    S = SReal.segment(0,S.Real.rows()/2).array()*std::complex<double>(1,0)+SReal.segment(S.Real.rows()/2,S.Real.rows()).array()*std::complex<double>(0,1);
    
  }

}


#endif /* complex_eigs_h */
