// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SETUP_INTEGRATION_H
#define DIRECTIONAL_SETUP_INTEGRATION_H

#include <queue>
#include <vector>
#include <cmath>
#include <unordered_set>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/tree.h>
#include <directional/principal_matching.h>
#include <directional/dcel.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/combing.h>
#include <directional/sparse_identity.h>
#include <directional/sparse_block.h>

namespace directional
{
//The data structure for seamless integration
struct IntegrationData
{
    int N;                                              // # uncompressed parametric functions
    int n;                                              // # independent parameteric functions
    Eigen::MatrixXi linRed;                             // Linear Reduction tying the n dofs to the full N
    Eigen::MatrixXi periodMat;                          // Function spanning integers
    Eigen::SparseMatrix<double> vertexTrans2CutMat;     // Map between the whole mesh (vertex + translational jump) representation to the vertex-based representation on the cut mesh
    Eigen::SparseMatrix<double> constraintMat;          // Linear constraints (resulting from non-singular nodes and feature alignment)
    //Eigen::SparseMatrix<double> constraintProjectMat;   // A constraint projection matrix
    Eigen::SparseMatrix<double> linRedMat;              // Global uncompression of n->N
    Eigen::SparseMatrix<double> intSpanMat;             // Spanning the translational jump lattice
    Eigen::SparseMatrix<double> singIntSpanMat;         // Layer for the singularities
    Eigen::SparseMatrix<double> featureAlignMat;        // The feature layer that projects a given cut-mesh solution to a feature-aligned one
    Eigen::SparseMatrix<double> featureAlignConstMat;   // A constraint matrix for the feature alignment (directly on the cut mesh)
    Eigen::VectorXi constrainedVertices;                // Constrained vertices (fixed points in the parameterization)
    Eigen::VectorXi integerVars;                        // Variables that are to be rounded.
    Eigen::MatrixXi face2cut;                           // |F|x3 map of which edges of faces are seams
    Eigen::VectorXd nVertexFunction;                    // Final compressed result (used for meshing)
    
    Eigen::VectorXi fixedIndices;                       // Translation fixing indices
    Eigen::VectorXd fixedValues;                        // Translation fixed values
    Eigen::VectorXi singularIndices;                    // Singular-vertex indices
    Eigen::VectorXi featureIndices;                     // Indices of the feature edges
    Eigen::VectorXi featureIntegerIndices;              //Indices (into the whole mesh + transitions) of integer variables from the feature alignment
    
    //integer versions, for exact seamless parameterizations (good for error-free meshing)
    Eigen::SparseMatrix<int> vertexTrans2CutMatInteger;
    Eigen::SparseMatrix<int> constraintMatInteger;
    Eigen::SparseMatrix<int> linRedMatInteger;
    Eigen::SparseMatrix<int> intSpanMatInteger;
    Eigen::SparseMatrix<int> singIntSpanMatInteger;
    //Eigen::SparseMatrix<int> FeatureMatInteger;
    
    //Debug variables
    Eigen::VectorXi featureFaces;
    Eigen::VectorXi featureInFaceVectors;
    
    double lengthRatio;                                 // Global scaling of functions
    //Flags
    bool integralSeamless;                              // Whether to do full translational seamless.
    bool roundSeams;                                    // Whether to round seams or round singularities
    bool verbose;                                       // Output the integration log.
    bool featureAlignment;                              // Whether to align to feature edges
    bool autoFeatureFunc;                               // Whether to automatically derive the function index that should be aligned (the one most orthogonal to the matched vectors on each side)
    
    IntegrationData(int _N):lengthRatio(0.02), integralSeamless(false), roundSeams(true), verbose(false), autoFeatureFunc(true), featureAlignment(false){
        N=_N;
        n=(N%2==0 ? N/2 : N);
        if (N%2==0)
            set_sign_symmetry(N);
        else linRed=Eigen::MatrixXi::Identity(N,n);
        set_default_period_matrix(n);
    }
    ~IntegrationData(){}
    
    inline void set_linear_reduction(const Eigen::MatrixXi& _linRed, const Eigen::MatrixXi& _periodMat){linRed =_linRed; N=linRed.rows(); n=linRed.cols(); periodMat=_periodMat;}
    
    //the default symmetry, where for even N there are N/2 lines
    inline void set_sign_symmetry(int N){
        assert(N%2==0);
        linRed.resize(N,N/2);
        linRed<<Eigen::MatrixXi::Identity(N/2,N/2),-Eigen::MatrixXi::Identity(N/2,N/2);
        n=N/2;
        set_default_period_matrix(n);
    }
    
    //the entire first N/3 lines are symmetric w.r.t. to the next two (N/3) packets, and where if N is even we also add sign symmetry.
    inline void set_triangular_symmetry(int N){
        assert(N%3==0);
        if (N%2==0){
            linRed.resize(N,N/3);
            linRed.block(0,0,N/2,N/3)<<Eigen::MatrixXi::Identity(N/3,N/3),-Eigen::MatrixXi::Identity(N/6,N/6),Eigen::MatrixXi::Identity(N/6,N/6);
            linRed.block(N/2,0,N/2,N/3)=-linRed.block(0,0,N/2,N/3);
            n=N/3;
        } else {
            linRed.resize(N,2*N/3);
            linRed<<Eigen::MatrixXi::Identity(2*N/3,2*N/3),-Eigen::MatrixXi::Identity(N/3,N/3),-Eigen::MatrixXi::Identity(N/3,N/3);
            n=2*N/3;
        }
        set_default_period_matrix(n);
    }
    
    inline void set_default_period_matrix(int n){
        periodMat=Eigen::MatrixXi::Identity(n,n);
    }
};


/**
 * Computes connected components directly from the constraint matrix C * y = 0.
 * * @param C  A sparse matrix of size (num_constraints x num_y), where each row
 * contains exactly one 1.0 and one -1.0.
 * @return   An Eigen::VectorXi of size num_y, where the i-th entry is the
 * contiguous, zero-indexed component ID (the z-index).
 */
Eigen::VectorXi dofsFromConstraintMatrix(const Eigen::SparseMatrix<double>& C) {
    const int numy = C.cols();
    
    // 1. Construct the adjacency topology: A = C^T * C
    // This creates a symmetric matrix where A(i,j) is non-zero if y_i and y_j
    // share an equality constraint.
    Eigen::SparseMatrix<double> A = C.transpose() * C;
    A.makeCompressed(); // Ensure optimal traversal memory layout

    Eigen::VectorXi dofMap = Eigen::VectorXi::Constant(numy, -1);
    int currentDofId = 0;
    
    // Reuse a single vector as a queue for memory efficiency during BFS
    std::vector<int> bfsQueue;
    bfsQueue.reserve(numy);

    // 2. Linear-time standard BFS traversal across the sparse graph topology
    for (int i = 0; i < numy; ++i) {
        // If this node is already assigned to a z-component, skip it
        if (dofMap(i) != -1) continue;

        // Found a new connected component! Start BFS.
        dofMap(i) = currentDofId;
        bfsQueue.push_back(i);
        
        size_t head = 0;
        while (head < bfsQueue.size()) {
            int u = bfsQueue[head++];
            
            // Traverse all structural neighbors of node 'u'
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, u); it; ++it) {
                int v = it.index();
                
                if (dofMap(v) == -1) {
                    dofMap(v) = currentDofId;
                    bfsQueue.push_back(v);
                }
            }
        }
        
        bfsQueue.clear();
        currentDofId++; // Increment for the next isolated variable group
    }

    return dofMap;
}

/**
 * Computes the sparse transfer operator T = (E^T * E)^-1 * E^T directly
 * from the connected component mappings.
 *
 * @param z_mapping  The VectorXi output from computeZMappingFromConstraintMatrix.
 * @return           A sparse matrix T of size (num_z x num_y).
 */
Eigen::SparseMatrix<double> featureProjectionMat(const Eigen::VectorXi& dofMapping) {
    const int numy = dofMapping.size();
    const int numz = dofMapping.maxCoeff() + 1;

    // 1. Group the original y-indices into their respective z-clusters
    std::vector<std::vector<int>> clusters(numz);
    for (int i = 0; i < numy; ++i) {
        clusters[dofMapping(i)].push_back(i);
    }

    // 2. Assemble the triplets for the square projection matrix P
    // Each cluster forms an isolated, fully connected dense block
    // where every entry inside that block has a weight of 1.0 / cluster_size.
    std::vector<Eigen::Triplet<double>> triplets;
    
    // Estimate allocation size to prevent vector reallocations
    size_t estimated_elements = 0;
    for (const auto& cluster : clusters) {
        estimated_elements += cluster.size() * cluster.size();
    }
    triplets.reserve(estimated_elements);

    for (const auto& cluster : clusters) {
        if (cluster.empty()) continue;
        
        double weight = 1.0 / static_cast<double>(cluster.size());
        
        // Populate the cross-product block for this component
        for (int row_idx : cluster) {
            for (int col_idx : cluster) {
                triplets.push_back(Eigen::Triplet<double>(row_idx, col_idx, weight));
            }
        }
    }

    // 3. Construct the square (numy x numy) filter matrix
    Eigen::SparseMatrix<double> P(numy, numy);
    P.setFromTriplets(triplets.begin(), triplets.end());
    P.makeCompressed();

    return P;
}

/**
 * Computes a sparse matrix P such that:
 *   - P * v is the 2-norm minimizing projection of v onto ker(C)
 *   - C * P = 0  (by construction)
 *   - P = I - C^T (C C^T)^+ C  (pseudoinverse formulation)
 *
 * Uses sparse QR on C^T to handle rank-deficient C of any size without
 * ever forming a dense k x k matrix.
 *
 * @param C  Sparse constraint matrix (k x n), any rank
 * @return   Sparse projection matrix P (n x n)
 */
Eigen::SparseMatrix<double> nullSpaceProjection(
    const Eigen::SparseMatrix<double>& C,
    double rankTolerance = 1e-10)
{
    using namespace Eigen;

    const int n = C.cols();
    const int k = C.rows();

    // P = I - C^T (C C^T)^+ C
    // We compute this implicitly via sparse QR of C^T:
    //
    // C^T = Q * R * P_col   (SparseQR decomposition)
    //
    // The first `rank` columns of Q span the row space of C.
    // P_C (projection onto row space) = Q1 * Q1^T
    // P_ker(C) = I - Q1 * Q1^T
    //
    // We materialise this as a sparse matrix by applying it to
    // the identity: P = I - Q1 * Q1^T, where Q1 is dense (n x rank)
    // but rank << n typically.

    // Step 1: sparse QR of C^T
    SparseMatrix<double> Ct = C.transpose();
    Ct.makeCompressed();

    SparseQR<SparseMatrix<double>, COLAMDOrdering<int>> sqr;
    sqr.compute(Ct);

    if (sqr.info() != Success)
        throw std::runtime_error("nullSpaceProjection: SparseQR failed");

    // Step 2: determine numerical rank
    const int rank = sqr.rank();

    if (rank == 0)
    {
        // C is zero matrix — projection is identity
        SparseMatrix<double> I(n, n);
        I.setIdentity();
        return I;
    }

    // Step 3: extract Q1 = first `rank` columns of Q (n x rank, dense)
    // matrixQ() is an implicit operator; materialise only the needed columns
    // by multiplying against the first `rank` standard basis vectors.
    //
    // Eigen's SparseQR returns Q as an implicit Householder sequence.
    // The safest way to get Q1 is to apply Q to I_{n x rank}.
    MatrixXd I_rank = MatrixXd::Identity(n, rank);
    MatrixXd Q1 = sqr.matrixQ() * I_rank;   // n x rank, dense

    // Step 4: P = I - Q1 * Q1^T  (as sparse)
    // Q1 * Q1^T is dense n x n — too large to form explicitly for big n.
    // Instead build P column by column:
    //   P[:,j] = e_j - Q1 * (Q1^T * e_j)
    //          = e_j - Q1 * Q1.row(j)^T
    // and store only entries above a sparsity threshold.
    //
    // For typical meshes rank << n, so Q1 * Q1^T is low-rank and
    // P is well approximated as sparse after thresholding near-zeros.

    // Build P = I - Q1*Q1^T column by column into triplets
    std::vector<Triplet<double>> triplets;
    triplets.reserve(n * (rank + 1));  // rough estimate

    const double thresh = 1e-12;

    for (int j = 0; j < n; j++)
    {
        // Q1^T * e_j = j-th row of Q1, as a column vector
        VectorXd q1j = Q1.row(j).transpose();   // rank x 1

        // Q1 * q1j  = projection component  (n x 1)
        VectorXd proj = Q1 * q1j;               // n x 1

        // Column j of P: e_j - proj
        for (int i = 0; i < n; i++)
        {
            double val = (i == j ? 1.0 : 0.0) - proj(i);
            if (std::fabs(val) > thresh)
                triplets.emplace_back(i, j, val);
        }
    }

    SparseMatrix<double> P(n, n);
    P.setFromTriplets(triplets.begin(), triplets.end());
    P.makeCompressed();
    return P;
}


/*** Setting up the seamless integration algorithm. Seamless integration only works on IntrinsicFaceTangentBundle at the moment.
 Input:
 field:        a face-based raw field that is to be integrated.
 Output:
 intData:      updated integration data.
 meshCut:      a mesh which is face-corresponding with meshWhole, but is cut so that it has disc-topology.
 combedField:  The raw field combed so that all singularities are on the seams of the cut mesh
 ***/
inline void setup_integration(const directional::CartesianField& field,
                              IntegrationData& intData,
                              directional::TriMesh& meshCut,
                              directional::CartesianField& combedField)
{
    
    using namespace Eigen;
    using namespace std;
    
    assert(field.tb->discTangType()==discTangTypeEnum::FACE_SPACES && "setup_integration() only works with face-based fields");
    
    const directional::TriMesh& meshWhole = *((PCFaceTangentBundle*)(field.tb))->mesh;
    //cutting mesh and combing field.
    cut_mesh_with_singularities(meshWhole, field.singLocalCycles, intData.face2cut);
    combing(field, combedField, intData.face2cut);
    
    // it stores number of edges per face, for now only tirangular
    //VectorXi D = VectorXi::Constant(meshWhole.F.rows(), 3);
    
    // mark vertices as being a singularity vertex of the vector field
    VectorXi isSingular = VectorXi::Zero(meshWhole.V.rows());
    for (int i = 0; i < field.singLocalCycles.size(); i++)
        isSingular(field.singLocalCycles(i)) = 1;

    
    intData.constrainedVertices = VectorXi::Zero(meshWhole.V.rows());
    
    // Find boundary vertices and mark them
    VectorXi isBoundary = VectorXi::Zero(meshWhole.V.rows());
    for (int i = 0; i < meshWhole.dcel.halfedges.size(); i++)
        if (meshWhole.twinH(i) == -1){
            isBoundary(meshWhole.HV(i)) = 1;
            isSingular(meshWhole.HV(i)) = 0; //boundary vertices cannot be singular
        }
    
    // Compute a permutation matrix
    vector<MatrixXi> constParmMatrices(intData.N);
    MatrixXi unitPermMatrix = MatrixXi::Zero(intData.N, intData.N);
    for (int i = 0; i < intData.N; i++)
        unitPermMatrix((i + 1) % intData.N, i) = 1;
    
    // Generate all the members of the permutation group
    constParmMatrices[0] = MatrixXi::Identity(intData.N, intData.N);
    for (int i = 1; i < intData.N; i++)
        constParmMatrices[i] = unitPermMatrix * constParmMatrices[i - 1];
    
    // Each edge which is on the cut seam is marked by 1 and 0 otherwise
    VectorXi isSeam = VectorXi::Zero(meshWhole.EV.rows());
    for(int i = 0; i < meshWhole.FE.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
            if (intData.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
                isSeam(meshWhole.FE(i, j)) = 1;
    }
    
    // Do the same for the half-edges, mark edges which correspond to the cut seam
    VectorXi isHEcut = VectorXi::Zero(meshWhole.dcel.halfedges.size());
    for(int i = 0; i < meshWhole.F.rows(); i++)
    {
        int hebegin = meshWhole.FH(i);
        int heiterate = hebegin;
        while (meshWhole.HV(hebegin)!=meshWhole.F(i,0))
            hebegin = meshWhole.nextH(hebegin);
        for (int j = 0; j < 3; j++) {
            if (intData.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
                isHEcut(heiterate) = 1; // FH is face to half-edge mapping
            heiterate = meshWhole.nextH(heiterate);
        }
    }
    
    // Calculate valence of the vertices which lay on the seam
    VectorXi cutValence = VectorXi::Zero(meshWhole.V.rows());
    for(int i = 0; i < meshWhole.EV.rows(); i++)
    {
        if (isSeam(i))
        {
            cutValence(meshWhole.EV(i, 0))++;
            cutValence(meshWhole.EV(i, 1))++;
        }
    }
    
    // Establish transition variables by tracing cut curves
    VectorXi Halfedge2TransitionIndices = VectorXi::Constant(meshWhole.dcel.halfedges.size(), 32767);
    VectorXi Halfedge2Matching(meshWhole.dcel.halfedges.size());
    VectorXi isHEClaimed = VectorXi::Zero(meshWhole.dcel.halfedges.size());
    
    // Convert the matching that was calculated for the vector field over edges to half-edges
    for (int i = 0; i < meshWhole.dcel.halfedges.size(); i++)
    {
        // HE is a map between half-edges to edges, but it does not carry the direction
        // EH edge to half-edge mapping
        Halfedge2Matching(i) = (meshWhole.EH(meshWhole.HE(i), 0) == i ? -combedField.matching(meshWhole.HE(i)) : combedField.matching(meshWhole.HE(i)));
        if(Halfedge2Matching(i) < 0)
            Halfedge2Matching(i) = (intData.N + (Halfedge2Matching(i) % intData.N)) % intData.N;
    }
    
    int currTransition = 1;
    
    /*****Cutting mesh and creating map between wholeF and cutF********/
    
    //cutting the mesh
    vector<int> cut2whole;
    vector<RowVector3d> cutVlist;
    MatrixXi cutF;
    MatrixXd cutV;
    cutF.resize(meshWhole.F.rows(),3);
    for (int i = 0; i < meshWhole.dcel.vertices.size(); i++)
    {
        //creating corners whereever we have non-trivial matching
        int beginH = meshWhole.VH(i);
        int currH = beginH;
        
        //reseting to first cut or first boundary, if exists
        if (!isBoundary(i))
        {
            do
            {
                if (isHEcut(currH)!=0)
                    break;
                currH=meshWhole.nextH(meshWhole.twinH(currH));
            } while (beginH!=currH);
        }
        else
        {
            do
            {
                if (meshWhole.twinH(currH)==-1)
                    break;
                currH=meshWhole.nextH(meshWhole.twinH(currH));
            } while(meshWhole.twinH(currH)!=-1);
        }
        
        beginH = currH;
        
        do
        {
            if ((isHEcut(currH) != 0) || (beginH == currH))
            {
                cut2whole.push_back(i);
                cutVlist.push_back(meshWhole.V.row(i));
            }
            
            for (int j = 0; j < 3; j++)
                if (meshWhole.F(meshWhole.HF(currH), j) == i)
                    cutF(meshWhole.HF(currH), j) = cut2whole.size() - 1;
            currH = meshWhole.twinH(meshWhole.prevH(currH));
        } while((beginH != currH) && (currH != -1));
    }
    
    cutV.resize(cutVlist.size(), 3);
    for(int i = 0; i < cutVlist.size(); i++)
        cutV.row(i) = cutVlist[i];
    
    //starting from each cut-graph node, we trace cut curves
    for(int i = 0;  i < meshWhole.V.rows(); i++)
    {
        if (((cutValence(i) == 2) && (!isSingular(i))) || (cutValence(i) == 0))
            continue;  //either mid-cut curve or non at all
        
        //tracing curves until next node, if not already filled
        int beginH = meshWhole.VH(i);
        
        //reseting to first boundary
        int currH = beginH;
        
        if (isBoundary(i))
        {
            do
            {
                if (meshWhole.twinH(currH) == -1)
                    break;
                currH = meshWhole.nextH(meshWhole.twinH(currH));
            } while(meshWhole.twinH(currH) != -1);
        }
        
        beginH = currH;
        
        int nextHalfedgeInCut = -1;
        do
        {
            //unclaimed inner halfedge
            if ((isHEcut(currH) != 0) && (isHEClaimed(currH) == 0) && (meshWhole.twinH(currH) != -1))
            {
                nextHalfedgeInCut = currH;
                Halfedge2TransitionIndices(nextHalfedgeInCut) = currTransition;
                Halfedge2TransitionIndices(meshWhole.twinH(nextHalfedgeInCut)) = -currTransition;
                isHEClaimed(nextHalfedgeInCut) = 1;
                isHEClaimed(meshWhole.twinH(nextHalfedgeInCut)) = 1;
                int nextCutVertex=meshWhole.HV(meshWhole.nextH(nextHalfedgeInCut));
                //advancing on the cut until next node
                while ((cutValence(nextCutVertex) == 2) && (!isSingular(nextCutVertex)) && (!isBoundary(nextCutVertex)))
                {
                    int beginH = meshWhole.VH(nextCutVertex);
                    int currH = beginH;
                    int nextHalfedgeInCut = -1;
                    do
                    {
                        //unclaimed cut halfedge
                        if ((isHEcut(currH) != 0) && (isHEClaimed(currH) == 0))
                        {
                            nextHalfedgeInCut = currH;
                            break;
                        }
                        currH=meshWhole.twinH(meshWhole.prevH(currH));
                    } while (beginH != currH);
                    Halfedge2TransitionIndices(nextHalfedgeInCut) = currTransition;
                    Halfedge2TransitionIndices(meshWhole.twinH(nextHalfedgeInCut)) = -currTransition;
                    isHEClaimed(nextHalfedgeInCut) = 1;
                    isHEClaimed(meshWhole.twinH(nextHalfedgeInCut)) = 1;
                    nextCutVertex = meshWhole.HV(meshWhole.nextH(nextHalfedgeInCut));
                }
                currTransition++;
            }
            currH = meshWhole.twinH(meshWhole.prevH(currH));
        } while((beginH != currH) && (currH != -1));
    }
    // end of cutting
    
    int numTransitions = currTransition - 1;
    vector<Triplet<double> > vertexTrans2CutTriplets, constTriplets;
    vector<Triplet<int> > vertexTrans2CutTripletsInteger, constTripletsInteger;
    // Forming the constraints and the singularity positions
    int currConst = 0;
    // Set up the transition functions (vector field matching) across the cuts
    for (int i = 0; i < meshWhole.V.rows(); i++)
    {
        std::vector<MatrixXi> permMatrices;
        std::vector<int> permIndices;  //in the space #V + #transitions
        //The initial corner gets the identity without any transition
        permMatrices.push_back(MatrixXi::Identity(intData.N, intData.N));
        permIndices.push_back(i);
        
        int beginH = meshWhole.VH(i);
        int currH = beginH;
        
        //reseting to first cut or boundary, if exists
        if (!isBoundary(i))
        {
            // travel throu the start of the vertex and stop once the edge on the cut is found
            do
            {
                if (isHEcut(currH) != 0)
                    break;
                currH = meshWhole.nextH(meshWhole.twinH(currH));
            } while(beginH != currH);
        }
        else
        {
            do
            {
                // travel until an edge without a twin is found, i.e., boundary
                if (meshWhole.twinH(currH) == -1)
                    break;
                currH = meshWhole.nextH(meshWhole.twinH(currH));
            } while(meshWhole.twinH(currH) != -1);
        }
        
        // set the beginning to the edge on the cut or on the boundary
        beginH = currH;
        
        int currCutVertex = -1;
        do
        {
            int currFace = meshWhole.HF(currH); // face containing the half-edge
            int newCutVertex = -1;
            //find position of the vertex i in the face of the initial mesh
            for (int j = 0; j < 3; j++)
            {
                if (meshWhole.F(currFace, j) == i)
                    newCutVertex = cutF(currFace, j);
            }
            
            //currCorner gets the permutations so far
            if (newCutVertex != currCutVertex)
            {
                currCutVertex = newCutVertex;
                for(int i = 0; i < permIndices.size(); i++)
                {
                    // place the perumtation matrix in a bigger matrix, we need to know how things are connected along the cut, no?
                    for(int j = 0; j < intData.N; j++)
                        for(int k = 0; k < intData.N; k++){
                            vertexTrans2CutTriplets.emplace_back(intData.N * currCutVertex + j, intData.N * permIndices[i] + k, (double) permMatrices[i](j, k));
                            vertexTrans2CutTripletsInteger.emplace_back(intData.N * currCutVertex + j, intData.N * permIndices[i] + k, permMatrices[i](j, k));
                        }
                }
            }
            
            //updating the matrices for the next corner
            int nextHalfedge = meshWhole.twinH(meshWhole.prevH(currH));
            //reached a boundary
            if(nextHalfedge == -1)
            {
                currH = nextHalfedge;
                continue;
            }
            
            // constParmMatrices contains all the members of the permutation group
            MatrixXi nextPermMatrix = constParmMatrices[Halfedge2Matching(nextHalfedge) % intData.N];
            //no update needed
            if(isHEcut(nextHalfedge) == 0)
            {
                currH = nextHalfedge;
                continue;
            }
            
            //otherwise, updating matrices with transition
            int nextTransition = Halfedge2TransitionIndices(nextHalfedge);
            //Pe*f + Je
            if(nextTransition > 0)
            {
                for(int j = 0; j < permMatrices.size(); j++)
                    permMatrices[j] = nextPermMatrix * permMatrices[j];
                
                //and identity on the fresh transition
                permMatrices.push_back(MatrixXi::Identity(intData.N, intData.N));
                permIndices.push_back(meshWhole.V.rows() + nextTransition - 1);
            }
            // (Pe*(f-Je))  matrix is already inverse since halfedge matching is minused
            else
            {
                //reverse order
                permMatrices.push_back(-MatrixXi::Identity(intData.N, intData.N));
                permIndices.push_back(meshWhole.V.rows() - nextTransition - 1);
                
                for(int j = 0; j < permMatrices.size(); j++)
                    permMatrices[j] = nextPermMatrix * permMatrices[j];
            }
            currH = nextHalfedge;
        } while((currH != beginH) && (currH != -1));
        
        //cleaning parmMatrices and permIndices to see if there is a constraint or reveal singularity-from-transition
        std::set<int> cleanPermIndicesSet(permIndices.begin(), permIndices.end());
        std::vector<int> cleanPermIndices(cleanPermIndicesSet.begin(), cleanPermIndicesSet.end());
        std::vector<MatrixXi> cleanPermMatrices(cleanPermIndices.size());
        
        for (int j = 0; j < cleanPermIndices.size(); j++)
        {
            cleanPermMatrices[j] = MatrixXi::Zero(intData.N, intData.N);
            for(int k = 0;k < permIndices.size(); k++)
                if(cleanPermIndices[j] == permIndices[k])
                    cleanPermMatrices[j] += permMatrices[k];
            if(cleanPermIndices[j] == i)
                cleanPermMatrices[j] -= MatrixXi::Identity(intData.N, intData.N);
        }
        
        //if not all matrices are zero, there is a constraint
        bool isConstraint = false;
        for(int j = 0; j < cleanPermMatrices.size(); j++)
            if (cleanPermMatrices[j].cwiseAbs().maxCoeff() != 0)
                isConstraint = true;
        
        if((isConstraint) && (!isBoundary(i)))
        {
            for(int j = 0; j < cleanPermMatrices.size(); j++)
            {
                for(int k = 0; k < intData.N; k++)
                    for(int l = 0; l < intData.N; l++){
                        constTriplets.emplace_back(intData.N * currConst + k, intData.N * cleanPermIndices[j] + l, (double) cleanPermMatrices[j](k, l));
                        constTripletsInteger.emplace_back(intData.N * currConst + k, intData.N * cleanPermIndices[j] + l, cleanPermMatrices[j](k, l));
                    }
            }
            currConst++;
            intData.constrainedVertices(i) = 1;
        }
    }
    
    vector< Triplet< double > > cleanTriplets;
    vector< Triplet< int > > cleanTripletsInteger;
    
    intData.vertexTrans2CutMat.resize(intData.N * cutV.rows(), intData.N * (meshWhole.V.rows() + numTransitions));
    intData.vertexTrans2CutMatInteger.resize(intData.N * cutV.rows(), intData.N * (meshWhole.V.rows() + numTransitions));
    cleanTriplets.clear();
    cleanTripletsInteger.clear();
    for(int i = 0; i < vertexTrans2CutTriplets.size(); i++){
        if(vertexTrans2CutTripletsInteger[i].value() != 0){
            cleanTripletsInteger.push_back(vertexTrans2CutTripletsInteger[i]);
            cleanTriplets.push_back(vertexTrans2CutTriplets[i]);
        }
        // if(std::abs((float)vertexTrans2CutTriplets[i].value())>10e-7)
    }
    intData.vertexTrans2CutMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    intData.vertexTrans2CutMatInteger.setFromTriplets(cleanTripletsInteger.begin(), cleanTripletsInteger.end());
    

    intData.constraintMat.resize(intData.N * currConst, intData.N * (meshWhole.V.rows() + numTransitions));
    intData.constraintMatInteger.resize(intData.N * currConst, intData.N * (meshWhole.V.rows() + numTransitions));
    cleanTriplets.clear();
    cleanTripletsInteger.clear();
    for(int i = 0; i < constTriplets.size(); i++){
        if(constTripletsInteger[i].value() != 0){
            cleanTripletsInteger.push_back(constTripletsInteger[i]);
            cleanTriplets.push_back(constTriplets[i]);
        }
    }
    intData.constraintMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    intData.constraintMatInteger.setFromTriplets(cleanTripletsInteger.begin(), cleanTripletsInteger.end());
    
    //doing the integer spanning matrix
    intData.intSpanMat.resize(intData.n * (meshWhole.V.rows() + numTransitions), intData.n * (meshWhole.V.rows() + numTransitions));
    intData.intSpanMatInteger.resize(intData.n * (meshWhole.V.rows() + numTransitions), intData.n * (meshWhole.V.rows() + numTransitions));
    vector<Triplet<double> > intSpanMatTriplets;
    vector<Triplet<int> > intSpanMatTripletsInteger;
    for (int i=0;i<intData.n*numTransitions;i+=intData.n){
        for(int k = 0; k < intData.n; k++)
            for(int l = 0; l < intData.n; l++){
                if (intData.periodMat(k,l)!=0){
                    intSpanMatTriplets.emplace_back(intData.n * meshWhole.V.rows()+i+k, intData.n * meshWhole.V.rows()+i+l, (double)intData.periodMat(k,l));
                    intSpanMatTripletsInteger.emplace_back(intData.n * meshWhole.V.rows()+i+k, intData.n * meshWhole.V.rows()+i+l, intData.periodMat(k,l));
                }
            }
    }
    for (int i=0;i<intData.n * meshWhole.V.rows();i++){
        intSpanMatTriplets.emplace_back(i,i,1.0);
        intSpanMatTripletsInteger.emplace_back(i,i,1);
    }
    
    intData.intSpanMat.setFromTriplets(intSpanMatTriplets.begin(), intSpanMatTriplets.end());
    intData.intSpanMatInteger.setFromTriplets(intSpanMatTripletsInteger.begin(), intSpanMatTripletsInteger.end());
    
    //filtering out barycentric symmetry, including sign symmetry. The parameterization should always only include n dof for the surface
    //TODO: this assumes n divides N!
    intData.linRedMat.resize(intData.N * (meshWhole.V.rows() + numTransitions), intData.n * (meshWhole.V.rows() + numTransitions));
    intData.linRedMatInteger.resize(intData.N * (meshWhole.V.rows() + numTransitions), intData.n * (meshWhole.V.rows() + numTransitions));
    vector<Triplet<double> > linRedMatTriplets;
    vector<Triplet<int> > linRedMatTripletsInteger;
    for(int i = 0; i < intData.N*(meshWhole.V.rows() + numTransitions); i +=intData.N)
        for(int k = 0; k < intData.N; k++)
            for(int l = 0; l < intData.n; l++){
                if (intData.linRed(k,l)!=0){
                    linRedMatTriplets.emplace_back(i + k, i*intData.n/intData.N + l, (double)intData.linRed(k,l));
                    linRedMatTripletsInteger.emplace_back(i + k, i*intData.n/intData.N + l, intData.linRed(k,l));
                }
            }
    
    intData.linRedMat.setFromTriplets(linRedMatTriplets.begin(), linRedMatTriplets.end());
    intData.linRedMatInteger.setFromTriplets(linRedMatTripletsInteger.begin(), linRedMatTripletsInteger.end());
    
    //integer variables are per single "d" packet, and the rounding is done for the N functions with projection over linRed
    intData.integerVars.resize(numTransitions);
    intData.integerVars.setZero();
    for(int i = 0; i < numTransitions; i++)
        intData.integerVars(i) = meshWhole.V.rows() + i;
    
    //fixed values
    intData.fixedIndices.resize(intData.n);
    if (isSingular.sum()==0){  //no inner singular vertices; vertex 0 is set to (0....0)
        for (int j=0;j<intData.n;j++)
            intData.fixedIndices(j)=j;
    }else {  //fixing first singularity to (0.5,....0.5)
        int firstSing;
        for (firstSing=0;firstSing<isSingular.size();firstSing++)
            if (isSingular(firstSing))
                break;
        //firstSing=0; //like before
        for (int j=0;j<intData.n;j++)
            intData.fixedIndices(j)=intData.n*firstSing+j;
    }
    
    //creating list of singular corners and singular integer matrix
    VectorXi singularIndices(intData.n * isSingular.sum());
    int counter=0;
    for (int i=0;i<isSingular.size();i++){
        if (isSingular(i))
            for (int j=0;j<intData.n;j++)
                singularIndices(counter++)=intData.n*i+j;
    }
    
    //doing the integer spanning matrix
    intData.singIntSpanMat.resize(intData.n * (meshWhole.V.rows() + numTransitions), intData.n * (meshWhole.V.rows() + numTransitions));
    intData.singIntSpanMatInteger.resize(intData.n * (meshWhole.V.rows() + numTransitions), intData.n * (meshWhole.V.rows() + numTransitions));
    vector<Triplet<double> > singIntSpanMatTriplets;
    vector<Triplet<int> > singIntSpanMatTripletsInteger;
    for (int i=0;i<isSingular.size();i++){
        if (!isSingular(i)){
            for (int j=0;j<intData.n;j++){
                singIntSpanMatTriplets.emplace_back(intData.n * i+j,intData.n * i+j,1.0);
                singIntSpanMatTripletsInteger.emplace_back(intData.n * i+j,intData.n * i+j,1);
            }
        } else {
            for(int k = 0; k < intData.n; k++)
                for(int l = 0; l < intData.n; l++){
                    if (intData.periodMat(k,l)!=0){
                        singIntSpanMatTriplets.emplace_back(intData.n*i+k, intData.n*i+l, (double)intData.periodMat(k,l));
                        singIntSpanMatTripletsInteger.emplace_back(intData.n*i+k, intData.n*i+l, intData.periodMat(k,l));
                    }
                }
        }
    }
    
    for (int i=intData.n * meshWhole.V.rows() ; i<intData.n*(meshWhole.V.rows()+numTransitions);i++){
        singIntSpanMatTriplets.emplace_back(i,i,1.0);
        singIntSpanMatTripletsInteger.emplace_back(i,i,1);
    }
    
    intData.singIntSpanMat.setFromTriplets(singIntSpanMatTriplets.begin(), singIntSpanMatTriplets.end());
    intData.singIntSpanMatInteger.setFromTriplets(singIntSpanMatTripletsInteger.begin(), singIntSpanMatTripletsInteger.end());
    
    intData.singularIndices=singularIndices;
    intData.fixedValues.resize(intData.n);
    intData.fixedValues.setConstant(0);
    
    meshCut.set_mesh(cutV, cutF);
    
    //Doing feature lines
    //This is a constraint matrix done directly on the cut mesh
    vector<int> featureCutIntegerIndices;
    if (intData.featureAlignment){
        vector<int> featureFaceList, featureInFaceVectorList;
        //featureAlignMat.resize(featureIndices.size(), intData.N*meshCut.V.rows());
        vector<Eigen::Triplet<double>> featureAlignMatTriplets;
        int featureIndex=0;
        for (int e=0;e<intData.featureIndices.size();e++){
            int faceLeft = meshWhole.EF(intData.featureIndices(e), 0);
            int faceRight = meshWhole.EF(intData.featureIndices(e), 1);
            int inFaceLeft = meshWhole.EFi(intData.featureIndices(e), 0);
            int inFaceRight = meshWhole.EFi(intData.featureIndices(e), 1);
            //std::cout<<"meshWhole.EV.row(intData.featureIndices(e)): "<<meshWhole.EV.row(intData.featureIndices(e))<<std::endl;
            //std::cout<<"meshWhole.F.row(faceLeft): "<<meshWhole.F.row(faceLeft)<<std::endl;
            //std::cout<<"meshWhole.F.row(faceRight): "<<meshWhole.F.row(faceRight)<<std::endl;
            //std::cout<<"inFaceLeft, inFaceRight: "<<inFaceLeft<<","<<inFaceRight<<std::endl;
            if (intData.autoFeatureFunc){
                //figuring out the directional field that is most orthogonal to the edge. This can be independent on both sides! forces agreement with matching
                double maxOrth = 2.0;
                int orthWhereLeft = -1, orthWhereRight = -1;
                RowVector3d normEdgeVector = (meshWhole.V.row(meshWhole.EV(intData.featureIndices(e),1)) - meshWhole.V.row(meshWhole.EV(intData.featureIndices(e),0))).normalized();
                //std::cout<<"normEdgeVector: "<<normEdgeVector<<std::endl;
                for (int i=0;i<intData.N;i++){
                    RowVector3d vecLeft = combedField.extField.block(faceLeft, i * 3, 1, 3).normalized();
                    //std::cout<<"vecLeft: "<<vecLeft<<std::endl;
                    double currOrth = abs(vecLeft.dot(normEdgeVector));
                    //std::cout<<"currOrth: "<<currOrth<<std::endl;
                    if (faceRight!=-1){
                        RowVector3d vecRight = combedField.extField.block(faceRight, ((i+combedField.matching(intData.featureIndices(e))+intData.N)%intData.N) * 3, 1, 3).normalized();
                        /*if (e==0){
                            std::cout<<"field.matching(intData.featureIndices(e)): "<<combedField.matching(intData.featureIndices(e))<<std::endl;
                            std::cout<<"normEdgeVector: "<<normEdgeVector<<std::endl;
                            std::cout<<"vecRight: "<<vecRight<<std::endl;
                            std::cout<<"vecLeft: "<<vecLeft<<std::endl;
                        }*/
                        //If this is sign symmetry, one would get chosen arbitrarily, and eventually this should not be important which
                        currOrth = (currOrth+abs(vecRight.dot(normEdgeVector)))/2.0;
                        //std::cout<<"currOrth: "<<currOrth<<std::endl;
                    }
                    if (currOrth<=maxOrth){
                        orthWhereLeft = i;
                        orthWhereRight = (i+combedField.matching(intData.featureIndices(e))+intData.N)%intData.N;
                        maxOrth = currOrth;
                    }
                }
                
                //std::cout<<"maxOrth: "<<maxOrth<<std::endl;
                
                /*if (e==0){
                    std::cout<<"orthWhereLeft: "<<orthWhereLeft<<std::endl;
                    std::cout<<"orthWhereRight: "<<orthWhereRight<<std::endl;
                }*/
                
                //adding a constraint for the proper differential on the edge to be zero. Is this wasteful?
                featureCutIntegerIndices.push_back(intData.N*meshCut.F(faceLeft, inFaceLeft) + orthWhereLeft);
                featureCutIntegerIndices.push_back(intData.N*meshCut.F(faceLeft, (inFaceLeft+1)%3) + orthWhereLeft);
                featureAlignMatTriplets.push_back(Triplet<double>(featureIndex, intData.N*meshCut.F(faceLeft, inFaceLeft) + orthWhereLeft, -1.0));
                featureAlignMatTriplets.push_back(Triplet<double>(featureIndex++, intData.N*meshCut.F(faceLeft, (inFaceLeft+1)%3) + orthWhereLeft, 1.0));
                featureFaceList.push_back(faceLeft);
                featureInFaceVectorList.push_back(orthWhereLeft);
                /*if (e==0){
                    std::cout<<"intData.N*meshCut.F(faceLeft, inFaceLeft) + orthWhereLeft: "<<intData.N*meshCut.F(faceLeft, inFaceLeft) + orthWhereLeft<<std::endl;
                    std::cout<<" intData.N*meshCut.F(faceLeft, (inFaceLeft+1)%3) + orthWhereLeft: "<< intData.N*meshCut.F(faceLeft, (inFaceLeft+1)%3) + orthWhereLeft<<std::endl;
                }*/
                if (faceRight!=-1){
                    featureCutIntegerIndices.push_back(intData.N*meshCut.F(faceRight, inFaceRight) + orthWhereRight);
                    featureCutIntegerIndices.push_back(intData.N*meshCut.F(faceRight, (inFaceRight+1)%3) + orthWhereRight);
                    featureAlignMatTriplets.push_back(Triplet<double>(featureIndex, intData.N*meshCut.F(faceRight, inFaceRight) + orthWhereRight, -1.0));
                    featureAlignMatTriplets.push_back(Triplet<double>(featureIndex++, intData.N*meshCut.F(faceRight, (inFaceRight+1)%3) + orthWhereRight, 1.0));
                    featureFaceList.push_back(faceRight);
                    featureInFaceVectorList.push_back(orthWhereRight);
                   /* if (e==0){
                        std::cout<<"intData.N*meshCut.F(faceRight, inFaceRight) + orthWhereRight: "<<intData.N*meshCut.F(faceRight, inFaceRight) + orthWhereRight<<std::endl;
                        std::cout<<"intData.N*meshCut.F(faceRight, (inFaceRight+1)%3) + orthWhereRight: "<< intData.N*meshCut.F(faceRight, (inFaceRight+1)%3) + orthWhereRight<<std::endl;
                        std::cout<<"<meshCut.V.row(meshCut.F(faceRight, inFaceRight)): "<<meshCut.V.row(meshCut.F(faceRight, inFaceRight))<<std::endl;
                        std::cout<<"<meshCut.V.row(meshCut.F(faceRight, (inFaceRight+1)%3)): "<<meshCut.V.row(meshCut.F(faceRight, (inFaceRight+1)%3))<<std::endl;
                        std::cout<<"meshWhole.V.row(meshWhole.EV(intData.featureIndices(e),0)): "<<meshWhole.V.row(meshWhole.EV(intData.featureIndices(e),0))<<std::endl;
                        std::cout<<"meshWhole.V.row(meshWhole.EV(intData.featureIndices(e),1)): "<<meshWhole.V.row(meshWhole.EV(intData.featureIndices(e),1))<<std::endl;
                        std::cout<<"combedField.extField.block(faceLeft, orthWhereLeft * 3, 1, 3): "<<combedField.extField.block(faceLeft, orthWhereLeft * 3, 1, 3)<<std::endl;
                        std::cout<<"combedField.extField.block(faceRight, orthWhereRight * 3, 1, 3): "<<combedField.extField.block(faceRight, orthWhereRight * 3, 1, 3)<<std::endl;
                        std::cout<<"meshCut.F(faceLeft, inFaceLeft): "<<meshCut.F(faceLeft, inFaceLeft)<<std::endl;
                        std::cout<<"meshCut.F(faceLeft, (inFaceLeft+1)%3): "<<meshCut.F(faceLeft, (inFaceLeft+1)%3)<<std::endl;
                        std::cout<<"meshCut.F(faceRight, inFaceRight): "<<meshCut.F(faceRight, inFaceRight)<<std::endl;
                        std::cout<<"meshCut.F(faceRight, (inFaceRight+1)%3): "<<meshCut.F(faceRight, (inFaceRight+1)%3)<<std::endl;
                        
                    }*/
                }
                    
            }
        }
        intData.featureAlignConstMat.resize(featureIndex, intData.N*meshCut.V.rows());
        intData.featureAlignConstMat.setFromTriplets(featureAlignMatTriplets.begin(), featureAlignMatTriplets.end());
        
        //Separating the DOFs of the features so they can be used as a layer matrix
        Eigen::VectorXi featureDofMap = dofsFromConstraintMatrix(intData.featureAlignConstMat);
        intData.featureAlignMat = featureProjectionMat(featureDofMap);
        //TODO: create integer version
        
        //Debug
        intData.featureFaces.resize(featureFaceList.size());
        for (int it=0;it<featureFaceList.size();it++)
            intData.featureFaces(it)=featureFaceList[it];
        
        intData.featureInFaceVectors.resize(featureInFaceVectorList.size());
        for (int it=0;it<featureInFaceVectorList.size();it++)
            intData.featureInFaceVectors(it)=featureInFaceVectorList[it];
        
        //trying to project out the constraint matrix
        //intData.constraintProjectMat = nullSpaceProjection(intData.constraintMat  *  intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat);
        
        //Test: doing it with constraints
        directional::sparse_identity(intData.N*meshCut.V.rows(), intData.N*meshCut.V.rows(), intData.featureAlignMat);
        Eigen::MatrixXi blockIndices(2,1); blockIndices<<0,1;
        vector<SparseMatrix<double>> bigC(2);
        bigC[0] = intData.constraintMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        bigC[1] = intData.featureAlignConstMat * intData.vertexTrans2CutMat * intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        SparseMatrix<double> fullC;
        directional::sparse_block(blockIndices, bigC, fullC);
        intData.constraintMat = fullC;
        
        //test: just making the relevant variables in the integer list
        //directional::sparse_identity(intData.N*meshCut.V.rows(), intData.N*meshCut.V.rows(), intData.featureAlignMat);
        //intData.constraintMat =  intData.constraintMat *  intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        //intData.constraintProjectMat = nullSpaceProjection(intData.constraintMat);
        //intData.featureAlignConstMat.resize(0,intData.N*meshCut.V.rows());
        
        //finding out all the variables that need to be integer
        /*Eigen::SparseMatrix<double, Eigen::RowMajor> FilterMat = intData.vertexTrans2CutMat *  intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        std::vector<int> featureIntegerList, mark(FilterMat.cols(), 0);
        for (int r : featureCutIntegerIndices){
            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(FilterMat, r); it; ++it){
                if (mark[it.col()] != 1){
                    mark[it.col()] = 1;
                    featureIntegerList.push_back(it.col());
                }
            }
        }
        
        
        std::unordered_set<int> S(featureIntegerList.begin(), featureIntegerList.end());
        //S.insert(intData.integerVars.data(), intData.integerVars.data() + intData.integerVars.size());
        std::vector<int> tmp(S.begin(), S.end());
        std::sort(tmp.begin(), tmp.end());
        intData.featureIntegerIndices.resize(tmp.size());
        for (int i = 0; i < tmp.size(); ++i)
            intData.featureIntegerIndices(i) = tmp[i];*/


    } else {
        directional::sparse_identity(intData.N*meshCut.V.rows(), intData.N*meshCut.V.rows(), intData.featureAlignMat);
        intData.constraintMat =  intData.constraintMat *  intData.linRedMat * intData.singIntSpanMat * intData.intSpanMat;
        //intData.constraintProjectMat = nullSpaceProjection(intData.constraintMat);
        intData.featureAlignConstMat.resize(0,intData.N*meshCut.V.rows());
       
    }
    
    
}
}

#endif


