// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CIRCUMCIRCLE_H
#define DIRECTIONAL_CIRCUMCIRCLE_H


namespace directional {



IGL_INLINE void circumcircle(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             Eigen::MatrixXd& circumcenter,
                             Eigen::VectorXd& circumradius)
  {
    //Using the formulation from https://en.wikipedia.org/wiki/Circumscribed_circle#Cartesian_coordinates_from_cross-_and_dot-products
    using namespace Eigen;
    
    circumcenter.resize(F.rows(),3);
    circumradius.resize(F.rows());
    for (int i=0;i<F.rows();i++){
      RowVector3d P12 = V.row(F(i,1))-V.row(F(i,0));
      RowVector3d P23 = V.row(F(i,2))-V.row(F(i,1));
      RowVector3d P31 = V.row(F(i,0))-V.row(F(i,2));
      
      double denominator = 2.0*P12.cross(P23).squaredNorm();
      double alpha = -P12.dot(P31)*P23.squaredNorm();
      double beta = -P12.dot(P23)*P31.squaredNorm();
      double gamma = -P31.dot(P23)*P12.squaredNorm();
      circumcenter.row(i) = V.row(F(i,0)).array()*alpha+V.row(F(i,1)).array()*beta+V.row(F(i,2)).array()*gamma;
      
      //testing to see if face center is correct
      double radius1 = (circumcenter.row(i) - V.row(F(i,0))).squaredNorm();
      double radius2 = (circumcenter.row(i) - V.row(F(i,1))).squaredNorm();
      double radius3 = (circumcenter.row(i) - V.row(F(i,2))).squaredNorm();
      
      std::cout<<"radii1,2,3: "<<radius1<<" "<<radius2<<" "<<radius3<<std::endl;
      
      circumradius(i) =  (circumcenter.row(i) - V.row(F(i,0))).norm();
    }
  }


}


#endif /* circumcircle_h */
