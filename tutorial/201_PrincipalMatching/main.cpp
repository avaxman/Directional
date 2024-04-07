#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOBJ.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/directional_viewer.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/pick.h>

int currF=0, N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field;
directional::DirectionalViewer viewer;


/*void update_raw_field_mesh()
{
  //configuring just part of the faces to have combed coloring
  Eigen::Vector3i otherFaces;
  Eigen::Vector3i zeroInFace;
  for (int i=0;i<3;i++){
    otherFaces(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? mesh.EF(mesh.FE(currF,i),1) : mesh.EF(mesh.FE(currF,i),0));
    zeroInFace(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? field.matching(mesh.FE(currF,i)) : -field.matching(mesh.FE(currF,i)));
  }
  
  Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(mesh.F.rows(),N);
  glyphColors.row(currF)=directional::DirectionalViewer::indexed_glyph_colors(field.extField.row(currF), false);
  for (int i=0;i<N;i++)
    for (int j=0;j<3;j++)
      glyphColors.block(otherFaces(j),3*((i+zeroInFace(j)+N)%N),1,3)<<glyphColors.row(currF).segment(3*i,3);
  
  viewer.set_field(field, glyphColors);
  
}*/

void callbackFunc()
{
    auto selection = polyscope::pick::getSelection();
    std::cout<<"selection: "<<selection.first<<" " <<selection.second<<std::endl;

}

int main()
{
  directional::readOBJ(TUTORIAL_DATA_PATH "/lilium.obj", mesh);
  ftb.init(mesh);
  viewer.init();
  directional::read_raw_field(TUTORIAL_DATA_PATH "/lilium.rawfield", ftb, N, field);
  directional::principal_matching(field);

  //triangle mesh setup
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  viewer.set_callback(callbackFunc);
  viewer.launch();
}

