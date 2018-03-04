#include <iostream>
#include <Eigen/Core>
#include <igl/viewer/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/get_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/principal_combing.h>


int currF=0, N;
Eigen::MatrixXi F;
Eigen::MatrixXd V, rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::RowVector3d rawGlyphColor;
igl::viewer::Viewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singPositions;
Eigen::MatrixXd positiveIndexColors, negativeIndexColors;

Eigen::MatrixXd glyphPrincipalColors(5,3);

bool showCombedField=false, showSingularities=false;


void update_mesh()
{
  
  Eigen::MatrixXd C(F.rows(),3);
  C.col(0)=Eigen::VectorXd::Constant(F.rows(),1.0);
  C.col(1)=Eigen::VectorXd::Constant(F.rows(),1.0);
  C.col(2)=Eigen::VectorXd::Constant(F.rows(),1.0);
  
  C.row(currF)<<0.5,0.1,0.1;
  
  Eigen::MatrixXd fullV=V;
  Eigen::MatrixXi fullF=F;
  Eigen::MatrixXd fullC=C;
  
  Eigen::MatrixXd fullGlyphColor(F.rows(),3*N);
  for (int i=0;i<F.rows();i++)
    for (int j=0;j<N;j++)
      fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
  
  directional::glyph_lines_raw(V, F, (showCombedField ? combedField : rawField), fullGlyphColor, false, true, fullV, fullF, fullC);
  
  if (showSingularities)
    directional::singularity_spheres(V, F, singPositions, singIndices, positiveIndexColors, negativeIndexColors, false, true, fullV, fullF, fullC);
  
  //drawing seam edges
  if (showCombedField){
    double l = igl::avg_edge_length(V, F);
    std::vector<int> seamEdges;
    for (int i=0;i<EV.rows();i++)
      if (combedMatching(i)!=0)
        seamEdges.push_back(i);
    
    Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
    for (int i=0;i<seamEdges.size();i++){
      P1.row(i)=V.row(EV(seamEdges[i],0));
      P2.row(i)=V.row(EV(seamEdges[i],1));
    }
  
    directional::line_cylinders(P1, P2, l/50.0, Eigen::MatrixXd::Constant(F.rows(), 3, 0.0), 6, false, true, fullV, fullF, fullC);
  }
  
  viewer.data.clear();
  viewer.data.set_face_based(true);
  viewer.data.set_mesh(fullV, fullF);
  viewer.data.set_colors(fullC);
}

// Handle keyboard input
bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': showCombedField=!showCombedField; update_mesh(); break;
    case '1': showSingularities=!showSingularities; update_mesh(); break;
  }
  return true;
}


int main()
{
  std::cout <<
  "  0  Toggle raw field/Combed field" << std::endl <<
  "  1  Show/hide singularities" << std::endl <<
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", N, rawField);
  igl::edge_topology(V, F, EV, FE, EF);
  
  //computing
  directional::principal_combing(V,F, rawField, combedField, combedMatching, combedEffort);
  directional::principal_matching(V, F,rawField, matching, effort);
  directional::get_indices(V,F,EV, EF, effort,N,prinIndices);
  
  std::vector<int> singPositionsList;
  std::vector<int> singIndicesList;
  for (int i=0;i<V.rows();i++)
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
  
  igl::barycenter(V, F, barycenters);
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0;
  
  // Set colors for Singularities
  positiveIndexColors.resize(4,3);
  positiveIndexColors << .25, 0, 0,
  .5,  0, 0,
  .75, 0, 0,
  1,   0, 0;
  
  negativeIndexColors.resize(4,3);
  negativeIndexColors << 0, .25, 0,
  0, .5,  0,
  0, .75, 0,
  0, 1,   0;
  
  update_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


