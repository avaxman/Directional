//
//  branched_gradient.h
//  PatternsParam_bin
//
//  Created by Amir Vaxman on 01/06/2020.
//

#ifndef DIRECTIONAL_BRANCHED_GRADIENT_H
#define DIRECTIONAL_BRANCHED_GRADIENT_H


namespace directional{
  
  
  //assuming corner functions are arranged packets of N per vertex.
  inline void branched_gradient(const TriMesh& mesh,
                                const int N,
                                Eigen::SparseMatrix<double>& G)
  {
    
    using namespace Eigen;
    using namespace std;
    
    VectorXd dblA = mesh.faceAreas*2.0;
    Eigen::MatrixXd normals = mesh.faceNormals;
    vector<Triplet<double>> GTri;
    
    for (int i=0;i<mesh.F.rows();i++){
      RowVector3d currNormal=normals.row(i);
      for (int k=0;k<N;k++){
        RowVector3d localGradient(0.0,0.0,0.0);
        for (int j=0;j<3;j++){
          RowVector3d eVec = mesh.V.row(mesh.F(i,(j+2)%3))-mesh.V.row(mesh.F(i,(j+1)%3));
          RowVector3d gradComp = currNormal.cross(eVec)/dblA(i);
          
          //double oppCornerFunction = cornerFunctions(i, k+N*((j+2)%3));
          //localGradient=localGradient+(oppCornerFunction/dblA(i))*eVecRot;
          for (int l=0;l<3;l++)
            GTri.push_back(Triplet<double>(3*N*i+k*3+l, N*mesh.F(i,j)+k,gradComp(l)));
          
        }
      }
    }
    G.conservativeResize(3*N*mesh.F.rows(), N*mesh.V.rows());
    //cout<<"G.rows(): "<<3*N*F.rows()<<endl;;
    G.setFromTriplets(GTri.begin(), GTri.end());
  }
}




#endif /* branched_gradient_h */
