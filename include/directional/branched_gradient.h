//
//  branched_gradient.h
//  PatternsParam_bin
//
//  Created by Amir Vaxman on 01/06/2020.
//

#ifndef branched_gradient_h
#define branched_gradient_h


namespace directional{
  
  
  //assuming corner functions are arranged packets of N per vertex.
  IGL_INLINE void branched_gradient(const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const int N,
                                    Eigen::SparseMatrix<double>& G)
  {
    
    using namespace Eigen;
    using namespace std;
    
    VectorXd dblA;
    igl::doublearea(V,F,dblA);
    Eigen::MatrixXd normals;
    igl::per_face_normals(V, F, normals);
    vector<Triplet<double>> GTri;
    
    for (int i=0;i<F.rows();i++){
      RowVector3d currNormal=normals.row(i);
      for (int k=0;k<N;k++){
        RowVector3d localGradient(0.0,0.0,0.0);
        for (int j=0;j<3;j++){
          RowVector3d eVec = V.row(F(i,(j+2)%3))-V.row(F(i,(j+1)%3));
          RowVector3d gradComp = currNormal.cross(eVec)/dblA(i);
          
          //double oppCornerFunction = cornerFunctions(i, k+N*((j+2)%3));
          //localGradient=localGradient+(oppCornerFunction/dblA(i))*eVecRot;
          for (int l=0;l<3;l++)
            GTri.push_back(Triplet<double>(3*N*i+k*3+l, N*F(i,j)+k,gradComp(l)));
          
        }
      }
    }
    G.conservativeResize(3*N*F.rows(), N*V.rows());
    //cout<<"G.rows(): "<<3*N*F.rows()<<endl;;
    G.setFromTriplets(GTri.begin(), GTri.end());
  }
}




#endif /* branched_gradient_h */
