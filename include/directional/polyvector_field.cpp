// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2017 Daniele Panozzo <daniele.panozzo@gmail.com>, Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <directional/polyvector_field.h>


#include <iterator>
#include <complex>
#include <cmath>
#include <stdexcept>

#include <Eigen/Sparse>
#include <Eigen/Geometry>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>

#include <igl/triangle_triangle_adjacency.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <igl/speye.h>
#include <igl/eigs.h>


namespace directional
{
  class PolyVectorComputer
  {

    // ----------------------- Types ------------------------------
  public:
    /// Sparse and complex valued matrix type
    typedef Eigen::SparseMatrix<std::complex<double> > SparseMatrix;
    /// An alias over a  complex valued triplet
    typedef Eigen::Triplet<std::complex<double> > Triplet;
    /// Container type storing triplets. It has to support a random access and back inserter iterators
    typedef std::vector<Triplet> TContainer;
    /// A random, const iterator over a container storing triplets
    typedef TContainer::const_iterator TCConstIterator;
    /// A back inserter over a container storing triplets
    typedef std::back_insert_iterator<TContainer> OutputIterator;


    // ----------------------- Standard services ------------------------------
  public:
    /**
     *
     * @param vertices #V by 3 list of the vertex positions
     * @param B1 each row represent a vector of the local coordinate basis at the respective face identified
     *        by the row number. Each such vector is orthogonal to a respective vector in B2.
     * @param B2 each row represent a vector of the local coordinate basis at the respective face identified
     *        by the row number. Each such vector is orthogonal to a respective vector in B1.
     * @param hardConstrIDs list of constrained faces' IDs. The number of rows must be equal
     *        to the number of rows of hardConstrDir.
     * @param hardConstrDir a list of constrained directions at the faces given in hardConstrIndices.
     *        Should be given in either as a matrix of size: no. of constrained faces by N, i.e.,
     *        X1, Y1, Z1, X2, Y2, Z2, Xn, Yn, Zn, or as a matrix of size: number of constrained face
     *        by 3, i.e., a single direction per row, implying N-RoSy. The number of rows must be equal
     *        to the number of rows in hardConstrIndices.
     * @param softConstrIDs list of constrained faces' IDs. The number of rows must be equal
     *        to the number of rows in softConstrWeights and softConstrDir.
     * @param softConstrWeights weights for the soft constraints. It can be given either as a matrix of size: no. of
     *        constrained face times N, or as a matrix of size number of constrained faces times 1. The number of
     *        rows must be equal to the number of rows in softConstrID and softConstrDir.
     * @param softConstrDir a list of constrained directions at the faces given in hardConstrIndices.
     *        Should be given in either as a matrix of size: no. of constrained faces by N, i.e.,
     *        X1, Y1, Z1, X2, Y2, Z2, Xn, Yn, Zn, or as a matrix of size: number of constrained face
     *        by 3, i.e., a single direction per row, implying N-RoSy. The number of
     *        rows must be equal to the number of rows in softConstrID and softConstrWeights.
     * @param N the degree of the field
     */
    PolyVectorComputer(const Eigen::MatrixXd & vertices,
                       const Eigen::MatrixXd & B1,
                       const Eigen::MatrixXd & B2,
                       const Eigen::VectorXi & hardConstrIDs,
                       const Eigen::MatrixXd & hardConstrDir,
                       const Eigen::VectorXi & softConstrIDs,
                       const Eigen::MatrixXd & softConstrWeights,
                       const Eigen::MatrixXd & softConstrDir,
                       unsigned int N);

    /**
     * It builds the system and precalculate the polyvector LDLt solvers.
     * @param EV a matrix of edges (vertex indices), its size has to be no. of edges times 2.
     *           Each row containts IDs of vertices which constitute the edge identified by the row.
     * @param EF a matrix of faces adjecent to edges, it size has to be no. of edges times 2.
     *           Each row contains IDs of faces adjacent to a given edge identified with a row.
     */
    void precompute(const Eigen::MatrixXi & EV, const Eigen::MatrixXi & EF);

    /**
     * It solves the system and returns the field. If no constraints (hard/soft) are given the Fielder eigenvector
     * field will be returned.
     * @param[out] polyVectorField  the output interpolated field, in polyvector (complex polynomial) format.
     *             The size of the output is the number of faces times N.
     */
    void eval(Eigen::MatrixXcd &polyVectorField);

    /**
     * Destructor.
     */
    ~PolyVectorComputer() = default;


    // ----------------------- Removed services ------------------------------
  public:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    PolyVectorComputer(const PolyVectorComputer & other) = delete;


    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    PolyVectorComputer & operator=(const PolyVectorComputer & other) = delete;


    // ------------------------- Hidden services ------------------------------
  private:

    /**
     * It is used to identify the variables i.e., terms that were not hard constrained.
     * A vector of size: the degree of the field times number of faces, is build and each hard constrained face
     * (identified by the row number) has value -1, and >= 0 otherwise.
     */
    void tagVariables();

    /**
     * Builds an internal representation of the hard constraints.
     */
    void treatHardConstraints();

    /**
     * Builds an internal representation of the soft constraints.
     */
    void treatSoftConstraints();

    /**
     *
     * @param EV a matrix of edges (vertex indices), its size has to be no. of edges times 2.
     *           Each row containts IDs of vertices which constitute the edge identified by the row.
     * @param EF a matrix of faces adjecent to edges, it size has to be no. of edges times 2.
     *           Each row contains IDs of faces adjacent to a given edge identified with a row.
     * @param result a back inserter iterator over a container storing Triplets.
     * @return number of rows of the energy matrix
     */
    unsigned int buildFullEnergyMatrix(const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF, OutputIterator result);

    /**
     * It extracts entries of the full energy matrix which corresponds to non-hard constrained faces.
     * @param begin an interator pointing at the first element of the container storing triples which
     *        correspond to the entries in the full energy matrix.
     * @param end an interator pointing at the end of the container storing triples which
     *        correspond to the entries in the full energy matrix.
     * @param result a writtable iterator over a container storing triples which correspond to non-hard
     *        constrained elements of the energy matrix.
     */
    void extractVarPartOfEnergy(TCConstIterator begin, TCConstIterator end, OutputIterator result);

    /**
     * It builds a diagonal matrix with diagonal entries corresponding to weights provided with the soft constraints.
     * @param result a writtable iterator over a container storing triples which correspond to a diagonal
     * matrix which has at the diagonal weights provided with the soft constraints.
     */
    void buildSoftConstraintsEnergyMatrix(OutputIterator result);

    /**
     * A wrap-up method which simply calls buildFullEnergyMatrix, extractVarPartOfEnergy
     * and buildSoftConstraintsEnergyMatrix.
     * @param EV a matrix of edges (vertex indices), its size has to be no. of edges times 2.
     *           Each row containts IDs of vertices which constitute the edge identified by the row.
     * @param EF a matrix of faces adjecent to edges, it size has to be no. of edges times 2.
     *           Each row contains IDs of faces adjacent to a given edge identified with a row.
     */
    void buildEnergyMatrics(const Eigen::MatrixXi & EV, const Eigen::MatrixXi & EF);

    /**
     * Computes the Fielder eigenvector field.
     * @param polyVectorField the output interpolated field, in polyvector (complex polynomial) format.
     *                        The size of the output is the number of faces times N.
     */
    void evalNoConstraints(Eigen::MatrixXcd & polyVectorField);


    // ------------------------- Private Data --------------------------------
  private:
    /// Internal representation of faces subjected to hard constraints.
    Eigen::VectorXi constIndices;
    /// Internal representation of directions at the faces subjected to hard constraints.
    Eigen::VectorXcd constValues;

    /// Internal representation of faces subjected to soft constraints.
    Eigen::VectorXi softIndices;
    /// Internal representation of directions at the faces subjected to soft constraints.
    Eigen::VectorXcd softValues;
    /// Internal representation of wieghts at the faces subjected to soft constraints.
    Eigen::VectorXcd softWeights;

    /// the size of this vector is equal to degree of the field times the no. of faces.
    /// Each element has value either -1 when the respective face (identified by the row number times 0 <= n < N)
    /// is subjected to hard constraints or value >= 0 which correspond to an ID of the respective
    /// variable in the system.
    Eigen::VectorXi full2var;

    /// the field degree
    unsigned int N;

    /// Internal reference to the mesh vertices
    const Eigen::MatrixXd & mVertices;
    //// Internal reference to the first vectors of the local basis at the faces
    const Eigen::MatrixXd & mB1;
    //// Internal reference to the second vectors of the local basis at the faces
    const Eigen::MatrixXd & mB2;
    /// Internal reference to the IDs of faces subjected to hard constraints
    const Eigen::VectorXi & mHardConstrIDs;
    /// Internal reference to the directions at faces subjected to hard constraints
    const Eigen::MatrixXd & mHardConstrDir;

    /// Internal reference to the IDs of faces subjected to soft constraints
    const Eigen::VectorXi & mSoftConstrIDs;
    /// Internal reference to the weights at faces subjected to soft constraints
    const Eigen::MatrixXd & mSoftConstrWeights;
    /// Internal reference to the directions at faces subjected to soft constraints
    const Eigen::MatrixXd & mSoftConstrDir;

    /// A matrix representing the full system
    SparseMatrix mAfull;
    /// A matrix representing only the part of the system which is not subjected to the hard constraints.
    SparseMatrix mAVar;
    /// A diagonal matrix of size: the field degree times no. of face times the field degree times no. of face, where
    /// each diagonal entry is equal to a weight at the respective face.
    SparseMatrix mASoft;

    /// Cholesky solver
    Eigen::SimplicialLDLT<SparseMatrix> mSolver;
  };


  PolyVectorComputer::PolyVectorComputer(const Eigen::MatrixXd & vertices,
                                         const Eigen::MatrixXd & B1,
                                         const Eigen::MatrixXd & B2,
                                         const Eigen::VectorXi & hardConstrIDs,
                                         const Eigen::MatrixXd & hardConstrDir,
                                         const Eigen::VectorXi & softConstrIDs,
                                         const Eigen::MatrixXd & softConstrWeights,
                                         const Eigen::MatrixXd & softConstrDir,
                                         unsigned int N) : mVertices(vertices), mB1(B1), mB2(B2),
                                         mHardConstrIDs(hardConstrIDs), mHardConstrDir(hardConstrDir),
                                         mSoftConstrDir(softConstrDir), mSoftConstrWeights(softConstrWeights),
                                         mSoftConstrIDs(softConstrIDs)
  {
    if (mHardConstrIDs.size() != mHardConstrDir.rows())
      throw std::runtime_error("directional::PolyVectorComputer: The hard constraines data are inconsistant!");

    if ((mSoftConstrDir.rows() + mSoftConstrIDs.rows() + mSoftConstrWeights.rows()) / 3. != mSoftConstrDir.rows())
      throw std::runtime_error("directional::PolyVectorComputer: The soft constraines data are inconsistant!");
    if(N == 0)
      throw std::runtime_error("directional::PolyVectorComputer: The field degree cannot be 0!");
    this->N = N;
  }

  void PolyVectorComputer::tagVariables()
  {
    // removing columns pertaining to constant indices
    Eigen::VectorXi varMask = Eigen::VectorXi::Constant(N * mB1.rows(), 1);
    for (unsigned int i = 0; i < constIndices.size(); i++)
      varMask(constIndices(i)) = 0;

    full2var = Eigen::VectorXi::Constant(N * mB1.rows(), -1);
    unsigned int varCounter = 0;
    for (unsigned int i = 0; i < N * mB1.rows(); i++)
      if (varMask(i))
        full2var(i) = varCounter++;
    assert(varCounter == N * (mB1.rows() - mHardConstrIDs.size()));
  }

  void PolyVectorComputer::treatHardConstraints()
  {
    Eigen::MatrixXcd constValuesMat(mHardConstrDir.rows(), N);
    assert((mHardConstrDir.cols() == 3 * N) || (mHardConstrDir.cols() == 3));
    if (mHardConstrDir.cols() == 3)  //N-RoSy constraint
    {
      constValuesMat.setZero();
      for (unsigned int i = 0; i < mHardConstrDir.rows(); i++) {
        std::complex<double> bComplex(mHardConstrDir.row(i).dot(mB1.row(mHardConstrIDs(i))),
                                      mHardConstrDir.row(i).dot(mB2.row(mHardConstrIDs(i))));
        constValuesMat(i, 0) = std::pow(bComplex, N);
      }
    }
    else
    {
      for (unsigned int i = 0; i < mHardConstrDir.rows(); i++)
      {
        Eigen::RowVectorXcd poly, roots(N);
        for (unsigned int n = 0; n < N; n++)
        {
          Eigen::RowVector3d vec = mHardConstrDir.block(i, 3 * n, 1, 3);
          roots(n) = std::complex<double>(vec.dot(mB1.row(mHardConstrIDs(i))), vec.dot(mB2.row(mHardConstrIDs(i))));
        }
        roots_to_monicPolynomial(roots, poly);
        constValuesMat.row(i) << poly.head(N);
      }
    }
    constValues.resize(N * mHardConstrDir.size());
    constValues.setZero();
    constIndices.resize(N * mHardConstrIDs.size());
    constIndices.setZero();
    for (unsigned int n = 0; n < N; n++)
    {
      constIndices.segment(mHardConstrIDs.rows() * n, mHardConstrIDs.rows()) = mHardConstrIDs.array() + n * mB1.rows();
      constValues.segment(mHardConstrDir.rows() * n, mHardConstrDir.rows()) = constValuesMat.col(n);
    }
  }

 void PolyVectorComputer::treatSoftConstraints()
 {
   Eigen::MatrixXcd softValuesMat(mSoftConstrDir.rows(), N);
   softValuesMat.setZero();
   Eigen::MatrixXcd softWeightsMat(mSoftConstrWeights.rows(), N);
   softWeightsMat.setZero();
   if (!((mSoftConstrDir.cols() == 3 * N) || (mSoftConstrDir.cols() == 3)) ||
       !((mSoftConstrWeights.cols() == N) || (mSoftConstrWeights.cols() == 1)))
     throw std::runtime_error("directional::PolyVectorComputer:::treatSoftConstraints: Missing information!");

   if (mSoftConstrDir.cols() == 3)  //N-RoSy constraint
   {
     for (unsigned int i = 0; i < mSoftConstrDir.rows(); i++)
     {
       std::complex<double> bComplex(mSoftConstrDir.row(i).dot(mB1.row(mSoftConstrIDs(i))),
                                     mSoftConstrDir.row(i).dot(mB2.row(mSoftConstrIDs(i))));
       softValuesMat(i, 0) = std::pow(mSoftConstrWeights(i, 0) * bComplex, N);
       softWeightsMat(i, 0) = std::complex<double>(mSoftConstrWeights(i, 0), 0);
     }
   }
   else
   {
     for (unsigned int i = 0; i < mSoftConstrDir.rows(); i++)
     {
       Eigen::RowVectorXcd poly, roots(N);
       for (unsigned int n = 0; n < N; n++)
       {
         Eigen::RowVector3d vec = mSoftConstrDir.block(i, 3 * n, 1, 3);
         roots(n) = mSoftConstrWeights(i, n) * std::complex<double>(vec.dot(mB1.row(mSoftConstrIDs(i))),
                                                                    vec.dot(mB2.row(mSoftConstrIDs(i))));
         softWeightsMat(i, n) = std::complex<double>(mSoftConstrWeights(i, n), 0);
       }
       roots_to_monicPolynomial(roots, poly);
       softValuesMat.row(i) << poly.head(N);
     }
   }

   softValues.resize(N * mSoftConstrDir.size());
   softValues.setZero();
   softWeights.resize(N * mSoftConstrWeights.rows());
   softWeights.setZero();
   softIndices.resize(N * mSoftConstrIDs.size());
   for (unsigned int n = 0; n < N; n++)
   {
     softIndices.segment(mSoftConstrIDs.rows() * n, mSoftConstrIDs.rows()) = mSoftConstrIDs.array() + n * mB1.rows();
     softWeights.segment(mSoftConstrWeights.rows() * n, mSoftConstrWeights.rows()) = softWeightsMat.col(n);
     softValues.segment(mSoftConstrDir.rows() * n, mSoftConstrDir.rows()) = softValuesMat.col(n);
   }
 }

 unsigned int PolyVectorComputer::buildFullEnergyMatrix(const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF,
                                                        PolyVectorComputer::OutputIterator result)

 {
   unsigned int rowCounter = 0;
   // Build the sparse matrix, with an energy term for each edge and degree
   for (unsigned int n = 0; n < N; n++)
   {
     for (unsigned int i = 0; i < EF.rows(); i++)
     {
       if ((EF(i, 0) == -1) || (EF(i, 1) == -1))
         continue;  //boundary edge

       // Compute the complex representation of the common edge
       Eigen::RowVector3d e = mVertices.row(EV(i, 1)) - mVertices.row(EV(i, 0));
       Eigen::RowVector2d vef = Eigen::Vector2d(e.dot(mB1.row(EF(i, 0))),
                                                e.dot(mB2.row(EF(i, 0)))).normalized();
       std::complex<double> ef(vef(0), vef(1));
       Eigen::Vector2d veg = Eigen::Vector2d(e.dot(mB1.row(EF(i, 1))),
                                             e.dot(mB2.row(EF(i, 1)))).normalized();
       std::complex<double> eg(veg(0), veg(1));

       // Add the term conj(f)^n*ui - conj(g)^n*uj to the energy matrix
       // ensure the sign consistance over the soft constraints
       *result++ = Triplet(rowCounter, n * mB1.rows() + EF(i, 0), std::pow(conj(ef), N - n));
       if ((softIndices.array() == EF(i, 0)).any() || (softIndices.array() == EF(i, 1)).any())
         *result++ = Triplet(rowCounter++, n * mB1.rows() + EF(i, 1), std::pow(conj(eg), N - n));
       else
         *result++ = Triplet(rowCounter++, n * mB1.rows() + EF(i, 1), -1. * std::pow(conj(eg), N - n));
     }
   }
   return rowCounter;
 }

 void PolyVectorComputer::extractVarPartOfEnergy(PolyVectorComputer::TCConstIterator begin,
                                                 PolyVectorComputer::TCConstIterator end,
                                                 PolyVectorComputer::OutputIterator result)
 {
   for (auto it = begin; it != end; ++it)
     if (full2var(it->col()) != -1)
       *result++ = Triplet(it->row(), full2var(it->col()), it->value());
 }

 void PolyVectorComputer::buildSoftConstraintsEnergyMatrix(
      std::back_insert_iterator<std::vector<Eigen::Triplet<std::complex<double> > > > result)
 {
   Eigen::VectorXcd soft(N * mB1.rows(), 1);
   soft.setZero();
   for (size_t i = 0; i < softIndices.size(); i++)
     soft(softIndices(i)) = softWeights(i);

   for (unsigned int i = 0; i < full2var.rows(); i++)
     if (full2var(i) != -1 && soft(i) != 0.)
       *result++ = Triplet(full2var(i), full2var(i), soft(i));
 }

 void PolyVectorComputer::buildEnergyMatrics(const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF)
 {
   TContainer AfullTriplets;
   unsigned int rowCounter = buildFullEnergyMatrix(EV, EF, OutputIterator(AfullTriplets));
   mAfull.resize(rowCounter, N * mB1.rows());
   mAfull.setFromTriplets(AfullTriplets.cbegin(), AfullTriplets.cend());

   // no need to compute this when no constraints
   if(softIndices.size() > 0 || constIndices.size() > 0)
   {
     TContainer AVarTriplets;
     extractVarPartOfEnergy(AfullTriplets.cbegin(), AfullTriplets.cend(), OutputIterator(AVarTriplets));
     mAVar.resize(rowCounter, N * (mB1.rows() - mHardConstrIDs.size()));
     mAVar.setFromTriplets(AVarTriplets.cbegin(), AVarTriplets.cend());
   }

   // no need to compute this when no constraints
   if(softIndices.size() > 0)
   {
     TContainer TSoft;
     buildSoftConstraintsEnergyMatrix(OutputIterator(TSoft));
     mASoft.resize(N * (mB1.rows() - mHardConstrIDs.size()), N * (mB1.rows() - mHardConstrIDs.size()));
     mASoft.setFromTriplets(TSoft.cbegin(), TSoft.cend());
   }
 }

 void PolyVectorComputer::evalNoConstraints(Eigen::MatrixXcd &polyVectorField)
 {
   //extracting first eigenvector into the field
   //Have to use reals because libigl does not currently support complex eigs.
   Eigen::SparseMatrix<double> M;
   igl::speye(2 * mB1.rows(), 2 * mB1.rows(), M);
   //creating a matrix of only the N-rosy interpolation
   SparseMatrix AfullNRosy(int((double)mAfull.rows() / (double) N), int((double) mAfull.cols() / (double) N));

   TContainer AfullNRosyTriplets;
   for (int k = 0; k < mAfull.outerSize(); ++k)
     for (SparseMatrix::InnerIterator it(mAfull, k); it; ++it)
       if ((it.row() < (double) mAfull.rows() / (double) N) &&
           (it.col() < (double) mAfull.cols() / (double) N))
         AfullNRosyTriplets.emplace_back(it.row(), it.col(), it.value());

   AfullNRosy.setFromTriplets(AfullNRosyTriplets.begin(), AfullNRosyTriplets.end());

   SparseMatrix LComplex = AfullNRosy.adjoint() * AfullNRosy;
   Eigen::SparseMatrix<double> L(2 * mB1.rows(), 2 * mB1.rows());
   std::vector<Eigen::Triplet<double> > LTriplets;
   for (unsigned int k = 0; k < LComplex.outerSize(); ++k)
     for (SparseMatrix::InnerIterator it(LComplex, k); it; ++it)
     {
       LTriplets.emplace_back(it.row(), it.col(), it.value().real());
       LTriplets.emplace_back(it.row(), LComplex.cols() + it.col(), -it.value().imag());
       LTriplets.emplace_back(LComplex.rows() + it.row(), it.col(), it.value().imag());
       LTriplets.emplace_back(LComplex.rows() + it.row(), LComplex.cols() + it.col(), it.value().real());
     }
   L.setFromTriplets(LTriplets.begin(), LTriplets.end());
   Eigen::MatrixXd U;
   Eigen::VectorXd S;
   igl::eigs(L, M, 5, igl::EIGS_TYPE_SM, U, S);

   polyVectorField = Eigen::MatrixXcd::Constant(mB1.rows(), N, std::complex<double>(0., 0.));
   polyVectorField.col(0) = U.block(0, 0, (long int) ((double) U.rows() / 2.), 1).cast<std::complex<double> >().array() * std::complex<double>(1., 0.) +
                            U.block(int((double) U.rows() / 2.), 0, int((double) U.rows() / 2.), 1).cast<std::complex<double> >().array() * std::complex<double>(0., 1.);
 }

 void PolyVectorComputer::precompute(const Eigen::MatrixXi &EV, const Eigen::MatrixXi &EF)
 {
   if (mHardConstrDir.rows() > 0)
     treatHardConstraints();
   if (mSoftConstrDir.rows() > 0)
     treatSoftConstraints();
   tagVariables();
   buildEnergyMatrics(EV, EF);

 }

 void PolyVectorComputer::eval(Eigen::MatrixXcd &polyVectorField)
 {
   if (mHardConstrIDs.size() == 0 && mSoftConstrWeights.size() == 0)
   {
     evalNoConstraints(polyVectorField);
     return;
   }

   assert(mSolver.rows() != 0);
   mSolver.compute(mAVar.adjoint() * mAVar + mASoft);
   Eigen::VectorXcd torhs(N * mB1.rows(), 1);
   torhs.setZero();
   for (size_t i = 0; i < constIndices.size(); i++)
     torhs(constIndices(i)) = constValues(i);

   Eigen::VectorXcd torhsSoft(N * mB1.rows(), 1);
   torhsSoft.setZero();
   for (size_t i = 0; i < softIndices.size(); i++)
     torhsSoft(softIndices(i)) = softValues(i);

   Eigen::VectorXcd rhs = -mAVar.adjoint() * mAfull * (torhs - torhsSoft);
   Eigen::VectorXcd varFieldVector = mSolver.solve(rhs);
   if (mSolver.info() != Eigen::Success)
     throw std::runtime_error("directional::PolyVectorComputer::eval: Solving the system finished with a failure!");

   // plug the hard constraints to the result
   Eigen::VectorXcd polyVectorFieldVector(N * mB1.rows());
   for (size_t i = 0; i < constIndices.size(); i++)
     polyVectorFieldVector(constIndices(i)) = constValues(i);

   // extract non-hard constrained results
   for (unsigned int i = 0; i < N * mB1.rows(); i++)
     if (full2var(i) != -1)
       polyVectorFieldVector(i) = varFieldVector(full2var(i));

   //converting to matrix form
   polyVectorField.conservativeResize(mB1.rows(), N);
   for (unsigned int n = 0; n < N; n++)
     polyVectorField.col(n) = polyVectorFieldVector.segment(n * mB1.rows(), mB1.rows());
 }


  IGL_INLINE void polyvector_field(const Eigen::MatrixXd & vertices,
                                   const Eigen::MatrixXi & faces,
                                   const Eigen::VectorXi & hardConstrIDs,
                                   const Eigen::MatrixXd & hardConstrDir,
                                   const Eigen::VectorXi & softConstrIDs,
                                   const Eigen::MatrixXd & softConstrWeights,
                                   const Eigen::MatrixXd & softConstrDir,
                                   unsigned int N,
                                   Eigen::MatrixXcd& polyVectorField)
  {
    Eigen::MatrixXi EV, EF;
    Eigen::MatrixXd B1, B2;
    // put this into a block to discard usless data
    {
      Eigen::MatrixXd B3;
      Eigen::MatrixXi FE;
      igl::edge_topology(vertices, faces, EV, FE, EF);
      igl::local_basis(vertices, faces, B1, B2, B3);
    }
    PolyVectorComputer pvComputer(vertices, B1, B2, hardConstrIDs, hardConstrDir, softConstrIDs,
                                  softConstrWeights, softConstrDir, N);
    pvComputer.precompute(EV, EF);
    pvComputer.eval(polyVectorField);
  }

  IGL_INLINE void polyvector_field(const Eigen::MatrixXd & vertices,
                                   const Eigen::MatrixXi & faces,
                                   Eigen::MatrixXd & B1,
                                   Eigen::MatrixXd & B2,
                                   Eigen::MatrixXi & EV,
                                   Eigen::MatrixXi & EF,
                                   const Eigen::VectorXi & hardConstrIDs,
                                   const Eigen::MatrixXd & hardConstrDir,
                                   const Eigen::VectorXi & softConstrIDs,
                                   const Eigen::MatrixXd & softConstrWeights,
                                   const Eigen::MatrixXd & softConstrDir,
                                   unsigned int N,
                                   Eigen::MatrixXcd & polyVectorField)
  {
    if(EV.rows() == 0 || EF.rows() == 0)
    {
      Eigen::MatrixXi FE;
      igl::edge_topology(vertices, faces, EV, FE, EF);
    }
    if(B1.rows() == 0 || B2.rows() == 0 || B1.rows() != B2.rows())
    {
      Eigen::MatrixXd B3;
      igl::local_basis(vertices, faces, B1, B2, B3);
    }

    PolyVectorComputer pvComputer(vertices, B1, B2, hardConstrIDs, hardConstrDir, softConstrIDs,
                                  softConstrWeights, softConstrDir, N);
    pvComputer.precompute(EV, EF);
    pvComputer.eval(polyVectorField);
  }
}
