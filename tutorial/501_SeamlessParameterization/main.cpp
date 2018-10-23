#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/principal_combing.h>
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/polyvector_field_cut_mesh_with_singularities.h>
#include <directional/integrable_polyvector_fields.h>


int N=6;
Eigen::MatrixXi FMeshWhole, FMeshCut;
Eigen::MatrixXd VMeshWhole, VCuteWhole, VField, VSings;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::RowVector3d rawGlyphColor;
igl::opengl::glfw::Viewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singVertices;
//Eigen::MatrixXd constPositions;
Eigen::MatrixXd glyphPrincipalColors(6,3);

Eigen::VectorXi cut2wholeIndices;  //map between cut vertices to whole vertices.
Eigen::VectorXi edge2TransitionIndices;  //map between all edges to transition variables (mostly -1; only relevant at cuts).
Eigen::MatrixXd cutUV; //, cutUVW;
Eigen::MatrixXi faceIsCut;

//parameterization variables
Eigen::SparseMatrix<double> vt2cMat;
Eigen::SparseMatrix<double> constraintMat, symmMat;

//between N vertex values and #cut transition variables to corner values (N*3 per face in same order)
Eigen::VectorXd edgeWeights;
Eigen::VectorXi constrainedVertices;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;


typedef enum {ROSY_FIELD, ROSY_PARAMETERIZATION, CF_FIELD, CF_PARAMETERIZATION} ViewingModes;
ViewingModes viewingMode=ROSY_FIELD;


void update_triangle_mesh()
{
  if ((viewingMode==ROSY_FIELD)||(viewingMode==CF_FIELD)){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
    CMesh=Eigen::MatrixXd::Constant(FMesh.rows(), 3, 1.0);
    viewer.data_list[0].set_colors(CMesh);
    viewer.data_list[0].show_texture=false;
  } else {
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshCut, FMeshCut);
    viewer.data_list[0].set_uv(cutUV);
    viewer.data_list[0].show_texture=true;
  }
}

void update_raw_field_mesh()
{
  if ((viewingMode==ROSY_PARAMETERIZATION)||(viewingMode==CF_PARAMETERIZATION)){
    for (int i=1;i<=3;i++){  //hide all other meshes
      viewer.data_list[i].show_faces=false;
      viewer.data_list[i].show_lines = false;
    }
  } else {
    Eigen::MatrixXd fullGlyphColor(FMeshWhole.rows(),3*N);
    for (int i=0;i<FMeshWhole.rows();i++)
      for (int j=0;j<N;j++)
        fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
    
    directional::glyph_lines_raw(VMeshWhole, FMeshWhole, (viewingMode==ROSY_FIELD ? combedFieldRosy : combedFieldCF), fullGlyphColor, VField, FField, CField);
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    //singularities mesh
    directional::singularity_spheres(VMeshWhole, FMeshWhole, singVertices, (viewingMode==ROSY_FIELD ? singIndicesRosy : singIndicesCF), directional::defaultSingularityColors(N) fullGlyphColor, VSings, FSings, CSings);
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings, FSings);
    viewer.data_list[2].set_colors(CSings);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    
    //seam mesh
    double avgEdgeLength = igl::avg_edge_length(VMesh, FMeshWhole);
    std::vector<int> seamEdges;
    for (int i=0;i<EV.rows();i++)
      if ((viewingMode==ROSY_FIELD ? combedMatchingRosy : combedMatchingCF)(i)!=0)
        seamEdges.push_back(i);
    
    Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
    for (int i=0;i<seamEdges.size();i++){
      P1.row(i)=VMeshWhole.row(EV(seamEdges[i],0));
      P2.row(i)=VMeshWhole.row(EV(seamEdges[i],1));
    }
    
    directional::line_cylinders(P1, P2, avgEdgeLength/25.0, Eigen::MatrixXd::Constant(FMesh.rows(), 3, 0.0), 6, VSeams, FSeams, CSeams);
    
    viewer.data_list[3].clear();
    viewer.data_list[3].set_mesh(VSeams, FSeams);
    viewer.data_list[3].set_colors(CSeams);
    viewer.data_list[3].show_faces = true;
    viewer.data_list[3].show_lines = false;
  }
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': showCombedField=!showCombedField; update_mesh(); break;
    case '1': showSingularities=!showSingularities; update_mesh(); break;
    case '2': viewer.data().show_texture=!viewer.data().show_texture; update_mesh(); break;
  }
  return true;
}


int main()
{
  std::cout <<
  "  0  Loaded 4-RoSy field" << std::endl <<
  "  1  Show/hide singularities" << std::endl <<
  "  2  Show textured mesh/original mesh" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka.rawfield", N, rawFieldRoSy);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-cf.rawfield", N, rawFieldCF);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  igl::barycenter(VMesh, FMesh, barycenters);
  
  //combing, cutting, and parameterizing according to RoSy field
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawFieldRosy, matchingRosy, effortRosy);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effortRosy,matchingRosy, N,prinIndices);
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<wholeV.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVerticesRosy.resize(singVerticesList.size());
  singIndicesRosy.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVerticesRosy(i)=singVerticesList[i];
    singIndicesRosy(i)=singIndicesList[i];
  }
  
  igl::polyvector_field_cut_mesh_with_singularities(VMesh, FMesh, singVerticesRosy, faceIsCutRosy);
  directional:principal_combing(VMesh,FMesh, EV, EF, FE, faceIsCutRosy, rawFieldRosy, combedFieldRosy, combedMatchingRosy, combedEffortRosy);
  
  //combing and cutting according to curl-free field
  directional::curl_matching(VMesh, FMesh,EV, EF, FE, rawFieldCF, matchingCF, effortCF);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effortCF,matchingCF, N,prinIndices);
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<wholeV.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVerticesCF.resize(singVerticesList.size());
  singIndicesCF.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVerticesCF(i)=singVerticesList[i];
    singIndicesCF(i)=singIndicesList[i];
  }
  
  igl::polyvector_field_cut_mesh_with_singularities(VMesh, FMesh, singVerticesCF, faceIsCutRosyCF);
  directional::curl_combing(VMesh,FMesh, EV, EF, FE, faceIsCutCF, rawFieldCF, combedFieldCF, combedMatchingCF, combedEffortCF);
  
  //Uniform weights for now
  edgeWeights=Eigen::VectorXd::Constant(EV.rows(), 1.0);
  cutF=wholeF;
  cutV=wholeV;
  Eigen::VectorXi integerVars;
  bool isBarycentric =true;  //to constrain u+v+w=1
  directional::setup_parameterization(N, wholeV, wholeF, combedMatching, singVertices, faceIsCut,  vt2cMat, constraintMat, symmMat, constrainedVertices, cutV, cutF, integerVars);
  
  std::vector<int> constPositionsList;
  for (int i=0;i<wholeV.rows();i++)
    if (constrainedVertices(i)!=0){
      constPositionsList.push_back(i);
    }
  
  constPositions.resize(constPositionsList.size(),3);
  for (int i=0;i<constPositionsList.size();i++){
    constPositions.row(i)=wholeV.row(constPositionsList[i]);
  }
  
  double edgeLength=200.0;
  bool isInteger = true;  //do not do translational seamless.
  integerVars.setZero();
  directional::parameterize(wholeV, wholeF, FE, combedField, edgeWeights, edgeLength, vt2cMat, constraintMat, symmMat, cutV, cutF, isInteger, integerVars, cutUVW, cutUV);
  
  //std::cout<<"cutUVW.col(0)-cutUVW.col(1)+cutUVW.col(2): "<<cutUVW.col(0)-cutUVW.col(1)+cutUVW.col(2)<<std::endl;
  
  Eigen::MatrixXd emptyMat;
  //igl::writeOBJ(TUTORIAL_SHARED_PATH "/torus-twosings-param.obj", cutV, cutF, emptyMat, emptyMat, cutUV, cutF);
  
  
  //testing vt2cMat
  
  
  //std::cout<<"cutUV: "<<cutUV<<std::endl;
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0,
  0.5,0.0,1.0;
  
  update_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


