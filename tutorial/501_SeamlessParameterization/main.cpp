#include <iostream>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/principal_combing.h>
#include <directional/curl_combing.h>
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/polyvector_field_cut_mesh_with_singularities.h>
#include <directional/integrable_polyvector_fields.h>


#include <directional/polyvector_field_matchings.h>
#include <directional/polyvector_field_singularities_from_matchings.h>
#include <directional/polyvector_field_cut_mesh_with_singularities.h>
#include <directional/polyvector_field_comb_from_matchings_and_cuts.h>
#include <directional/polyvector_field_poisson_reconstruction.h>



int N=4;
Eigen::MatrixXi wholeF, cutF;
Eigen::MatrixXd wholeV, cutV, rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::RowVector3d rawGlyphColor;
igl::viewer::Viewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singPositions;
Eigen::MatrixXd constPositions;
Eigen::MatrixXd glyphPrincipalColors(5,3);

Eigen::VectorXi cut2wholeIndices;  //map between cut vertices to whole vertices.
Eigen::VectorXi edge2TransitionIndices;  //map between all edges to transition variables (mostly -1; only relevant at cuts).
Eigen::MatrixXd cutUV;
Eigen::MatrixXi faceIsCut;

Eigen::SparseMatrix<double> vt2cMat;
Eigen::SparseMatrix<double> constraintMat, symmMat;
//between N vertex values and #cut transition variables to corner values (N*3 per face in same order)
Eigen::VectorXd edgeWeights;

Eigen::VectorXi constrainedVertices;

bool showCombedField=false, showSingularities=false;


void update_mesh()
{

  Eigen::MatrixXd fullV;
  Eigen::MatrixXi fullF;
  Eigen::MatrixXd fullC;
  
  if (!viewer.core.show_texture){
    fullC=Eigen::MatrixXd::Constant(cutF.rows(),3, 1.0);
    fullV=wholeV;
    fullF=wholeF;
    
    Eigen::MatrixXd fullGlyphColor(wholeF.rows(),3*N);
    for (int i=0;i<wholeF.rows();i++)
      for (int j=0;j<N;j++)
        fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
    
    directional::glyph_lines_raw(wholeV, wholeF, (showCombedField ? combedField : rawField), fullGlyphColor, false, true, fullV, fullF, fullC);
    
    if (showSingularities)
      directional::singularity_spheres(wholeV, wholeF, singPositions, singIndices, directional::defaultSingularityColors(N), false, true, fullV, fullF, fullC);
    
    double l = igl::avg_edge_length(wholeV, wholeF);
    Eigen::MatrixXd constColors(constPositions.rows(),3);
    constColors.setConstant(0.5);
    directional::point_spheres(constPositions, l/5.0, constColors, 8,  false, true, fullV, fullF, fullC);
    
    //drawing seam edges
    if (showCombedField){
      double l = igl::avg_edge_length(wholeV, wholeF);
      std::vector<int> seamEdges;
      /*for (int i=0;i<EV.rows();i++)
        if (combedMatching(i)!=0)
          seamEdges.push_back(i);*/
      
      for (int i=0;i<wholeF.rows();i++)
        for (int j=0;j<3;j++)
          if (faceIsCut(i,j))
            seamEdges.push_back(FE(i,j));
      
      Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
      for (int i=0;i<seamEdges.size();i++){
        P1.row(i)=wholeV.row(EV(seamEdges[i],0));
        P2.row(i)=wholeV.row(EV(seamEdges[i],1));
      }
      
      directional::line_cylinders(P1, P2, l/50.0, Eigen::MatrixXd::Constant(wholeF.rows(), 3, 0.0), 6, false, true, fullV, fullF, fullC);
    }
    viewer.data.clear();
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(fullV, fullF);
    viewer.data.set_colors(fullC);
  } else {
    fullV=cutV;
    fullF=cutF;
    
    
    
    viewer.data.clear();
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(fullV, fullF);
    viewer.data.set_uv(cutUV);
    
  }
  
  
}

// Handle keyboard input
bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': showCombedField=!showCombedField; update_mesh(); break;
    case '1': showSingularities=!showSingularities; update_mesh(); break;
    case '2': viewer.core.show_texture=!viewer.core.show_texture; update_mesh(); break;
  }
  return true;
}


int main()
{
  std::cout <<
  "  0  Toggle raw field/Combed field" << std::endl <<
  "  1  Show/hide singularities" << std::endl <<
  "  2  Show textured mesh/original mesh" << std::endl <<
  igl::readOFF(TUTORIAL_SHARED_PATH "/decimated-knight.off", wholeV, wholeF);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/decimated-knight.rawfield", N, rawField);
  igl::edge_topology(wholeV, wholeF, EV, FE, EF);
  
  Eigen::MatrixXd two_pv=rawField.block(0,0,rawField.rows(),6);
  Eigen::VectorXi b; b.resize(1); b<<0;
  Eigen::MatrixXd bc; bc.resize(1,6); bc<<two_pv.row(0);
  double constraint_percentage = 0.002;
  double rand_factor = 5;
  int iter = 0;
  Eigen::VectorXi blevel; blevel.resize(1); b<<1;
  igl::integrable_polyvector_fields_parameters params;
  igl::IntegrableFieldSolverData<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd> ipfdata;
  igl::integrable_polyvector_fields_precompute(wholeV, wholeF, b, bc, blevel, two_pv, ipfdata);
  for (int bi = 0; bi<50; ++bi)
  {
    printf("\n\n **** Batch %d ****\n", iter);
    igl::integrable_polyvector_fields_solve(ipfdata, params, two_pv, iter ==0);
    iter++;
    params.wSmooth *= params.redFactor_wsmooth;
  }
  
  printf("--Matchings and curl--\n");
  Eigen::MatrixXi match_ab, match_ba;  // matchings across interior edges
  Eigen::VectorXd curl; // curl per edge
  double avgCurl = igl::polyvector_field_matchings(two_pv, wholeV, wholeF, true, true, match_ab, match_ba, curl);
  double maxCurl = curl.maxCoeff();
  printf("curl -- max: %.5g, avg: %.5g\n", maxCurl,  avgCurl);
  // Compute singularities
  printf("--Singularities--\n");
  igl::polyvector_field_singularities_from_matchings(wholeV, wholeF, match_ab, match_ba, singPositions);
  printf("#singularities: %ld\n", singPositions.rows());
  // Get mesh cuts based on singularities
  printf("--Cuts--\n");
  igl::polyvector_field_cut_mesh_with_singularities(wholeV, wholeF, singPositions, faceIsCut);
  // Comb field
  printf("--Combing--\n");
  Eigen::MatrixXd combed;
  igl::polyvector_field_comb_from_matchings_and_cuts(wholeV, wholeF, two_pv, match_ab, match_ba, faceIsCut, combed);
  
  rawField.block(0,0,rawField.rows(),6)=two_pv;
  rawField.block(0,6,rawField.rows(),6)=-two_pv;
  
  singIndices.resize(singPositions.size());
  singIndices.setZero();
  
    
  //computing
  //directional::principal_matching(wholeV, wholeF,EV, EF, FE, rawField, matching, effort);
  directional::curl_matching(wholeV, wholeF,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(wholeV,wholeF,EV, EF, effort,N,prinIndices);
  
  //cutting and parameterizing
  
  std::vector<int> singPositionsList;
  std::vector<int> singIndicesList;
  for (int i=0;i<wholeV.rows();i++)
    if (prinIndices(i)!=0){
      singPositionsList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singPositions.resize(singPositionsList.size());
  singIndices.resize(singIndicesList.size());
  for (int i=0;i<singPositionsList.size();i++){
    singPositions(i)=singPositionsList[i];
    singIndices(i)=singIndicesList[i];
  }
  
  
  igl::polyvector_field_cut_mesh_with_singularities(wholeV, wholeF, singPositions, faceIsCut);
  directional::curl_combing(wholeV,wholeF, EV, EF, FE, faceIsCut, rawField, combedField, combedMatching, combedEffort);
  //std::cout<<"combedMatching: "<<combedMatching<<std::endl;
  
  igl::barycenter(wholeV, wholeF, barycenters);
  
  //Uniform weights for now
  edgeWeights=Eigen::VectorXd::Constant(EV.rows(), 1.0);
  cutF=wholeF;
  cutV=wholeV;
  Eigen::VectorXi integerVars;
  directional::setup_parameterization(N, wholeV, wholeF, combedMatching, singPositions, faceIsCut, vt2cMat, constraintMat, symmMat, constrainedVertices, cutV, cutF, integerVars);
  
  std::vector<int> constPositionsList;
  for (int i=0;i<wholeV.rows();i++)
    if (constrainedVertices(i)!=0){
      constPositionsList.push_back(i);
    }
  
  constPositions.resize(constPositionsList.size(),3);
  for (int i=0;i<constPositionsList.size();i++){
    constPositions.row(i)=wholeV.row(constPositionsList[i]);
  }
  
  double edgeLength=50.0;
  //integerVars.setZero();
  directional::parameterize(wholeV, wholeF, FE, combedField, edgeWeights, edgeLength, vt2cMat, constraintMat, symmMat, cutV, cutF, integerVars, cutUV);
  
  Eigen::MatrixXd emptyMat;
  igl::writeOBJ(TUTORIAL_SHARED_PATH "/decimated-knight-param.obj", cutV, cutF, emptyMat, emptyMat, cutUV, cutF);

  
  //testing vt2cMat
  
  
  //std::cout<<"cutUV: "<<cutUV<<std::endl;
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0;
  
  update_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


