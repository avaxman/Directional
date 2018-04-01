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


int currF=0, N;
Eigen::MatrixXi F;
Eigen::MatrixXd V, rawField, barycenters;
Eigen::VectorXd effort;
Eigen::RowVector3d rawGlyphColor;
igl::viewer::Viewer viewer;
Eigen::VectorXi matching, indices;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singPositions;
Eigen::MatrixXd positiveIndexColors, negativeIndexColors;

Eigen::MatrixXd glyphPrincipalColors(5,3);

bool zeroPressed=false, showSingularities=false;

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
  
  Eigen::Vector3i otherFaces;
  Eigen::Vector3i zeroInFace;
  for (int i=0;i<3;i++){
    otherFaces(i)=(EF(FE(currF,i),0)==currF ? EF(FE(currF,i),1) : EF(FE(currF,i),0));
    zeroInFace(i)=(EF(FE(currF,i),0)==currF ? matching(FE(currF,i)) : -matching(FE(currF,i)));
  }
  
  Eigen::MatrixXd fullGlyphColor(F.rows(),3*N);
  for (int i=0;i<F.rows();i++)
    for (int j=0;j<N;j++)
      fullGlyphColor.block(i,3*j,1,3)<<rawGlyphColor;
  
  //for marked face and its neighbors, coloring adjacent matchings
  for (int i=0;i<N;i++){
    fullGlyphColor.block(currF,3*i,1,3)<<glyphPrincipalColors.row(i);
    for (int j=0;j<3;j++)
      if (otherFaces(j)!=-1)  //boundary edge
        fullGlyphColor.block(otherFaces(j),3*((i+zeroInFace(j)+N)%N),1,3)<<glyphPrincipalColors.row(i);
  }
  
  directional::glyph_lines_raw(V, F, rawField, fullGlyphColor, false, true, fullV, fullF, fullC);
  
  if (showSingularities)
    directional::singularity_spheres(V, F, singPositions, singIndices, positiveIndexColors, negativeIndexColors, false, true, fullV, fullF, fullC);
  
  viewer.data.clear();
  viewer.data.set_face_based(true);
  viewer.data.set_mesh(fullV, fullF);
  viewer.data.set_colors(fullC);
}

bool key_up(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=false; break;
  }
  return true;
}

// Handle keyboard input
bool key_down(igl::viewer::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
    case '1': showSingularities=!showSingularities; update_mesh(); break;
  }
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::viewer::Viewer& viewer, int button, int modifiers)
{
  if (!zeroPressed)
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view * viewer.core.model,
                               viewer.core.proj, viewer.core.viewport, V, F, fid, bc))
  {
    
    //choosing face
    if ((igl::viewer::Viewer::MouseButton)button==igl::viewer::Viewer::MouseButton::Left){
      currF=fid;
      update_mesh();
      return true;
    }
  }
  return false;
};

int main()
{
  std::cout <<
  "  0+Left button    Choose face" << std::endl <<
  "  1  Show/hide singularities" << std::endl <<
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", N, rawField);
  igl::edge_topology(V, F, EV, FE, EF);
  
  //computing
  directional::principal_matching(V, F,rawField, matching, effort);
  directional::effort_to_indices(V,F,EV, EF, effort,N,prinIndices);
  
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
  
  rawGlyphColor <<0.0, 0.2, 1.0;
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
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  update_mesh();
  viewer.launch();
}

