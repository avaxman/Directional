#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/glyph_lines_raw.h>
#include <directional/line_cylinders.h>
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


int N;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams;
Eigen::MatrixXd CField, CSeams, CSings;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::RowVector3d rawGlyphColor;
igl::opengl::glfw::Viewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi prinIndices;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd glyphPrincipalColors(6,3);
Eigen::MatrixXd cutUVFull, cutUVRot;


typedef enum {FIELD, ROT_PARAMETERIZATION, FULL_PARAMETERIZATION} ViewingModes;
ViewingModes viewingMode=FIELD;


void update_triangle_mesh()
{
  if (viewingMode==FIELD){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
    viewer.data_list[0].set_colors(Eigen::RowVector3d(1.0,1.0,1.0));
    viewer.data_list[0].show_texture=false;
  } else {
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshCut, FMeshCut);
    viewer.data_list[0].set_uv(viewingMode==ROT_PARAMETERIZATION ? cutUVRot : cutUVFull);
    viewer.data_list[0].show_texture=true;
  }
}

void update_raw_field_mesh()
{
  if ((viewingMode=ROT_PARAMETERIZATION)||(viewingMode=FULL_PARAMETERIZATION)){
    for (int i=1;i<=3;i++){  //hide all other meshes
      viewer.data_list[i].show_faces=false;
      viewer.data_list[i].show_lines = false;
    }
  } else {
    Eigen::MatrixXd fullGlyphColor(FMeshWhole.rows(),3*N);
    for (int i=0;i<FMeshWhole.rows();i++)
      for (int j=0;j<N;j++)
        fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
    
    directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, fullGlyphColor, VField, FField, CField);
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    //singularities mesh
    directional::singularity_spheres(VMeshWhole, FMeshWhole, singVertices, singIndices, directional::defaultSingularityColors(N), VSings, FSings, CSings);
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings, FSings);
    viewer.data_list[2].set_colors(CSings);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    
    //seam mesh
    double avgEdgeLength = igl::avg_edge_length(VMeshWhole, FMeshWhole);
    std::vector<int> seamEdges;
    for (int i=0;i<EV.rows();i++)
      if (combedMatching(i)!=0)
        seamEdges.push_back(i);
    
    Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
    for (int i=0;i<seamEdges.size();i++){
      P1.row(i)=VMeshWhole.row(EV(seamEdges[i],0));
      P2.row(i)=VMeshWhole.row(EV(seamEdges[i],1));
    }
    
    directional::line_cylinders(P1, P2, avgEdgeLength/25.0, Eigen::MatrixXd::Constant(FMeshWhole.rows(), 3, 0.0), 6, VSeams, FSeams, CSeams);
    
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
    case '0': viewingMode = FIELD; break;
    case '1': viewingMode = ROT_PARAMETERIZATION; break;
    case '2': viewingMode = FULL_PARAMETERIZATION; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/dragon-param-rot-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVRot, FMeshCut);
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/dragon-param-full-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVFull, FMeshCut);
      break;
    //case '2': viewer.data().show_texture=!viewer.data().show_texture; update_mesh(); break;
  }
  update_triangle_mesh();
  update_raw_field_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  0  Loaded field" << std::endl <<
  "  1  Show textured rotationally-seamless parameterization mesh" << std::endl <<
  "  2  Show textured fully-seamless parameterization mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  

  igl::readOFF(TUTORIAL_SHARED_PATH "/dragon4.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/dragon4.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  //combing and cutting
  directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,prinIndices);
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<VMeshWhole.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVertices.resize(singVerticesList.size());
  singIndices.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVertices(i)=singVerticesList[i];
    singIndices(i)=singIndicesList[i];
  }
  
  directional::ParameterizationData pd;
  igl::polyvector_field_cut_mesh_with_singularities(VMeshWhole, FMeshWhole, singVertices, pd.face2cut);
  directional::principal_combing(VMeshWhole,FMeshWhole, EV, EF, FE, pd.face2cut, rawField, combedField, combedMatching, combedEffort);
  
  std::cout<<"Setting up parameterization"<<std::endl;
  
  directional::setup_parameterization(N, VMeshWhole, FMeshWhole, combedMatching, singVertices, pd, VMeshCut, FMeshCut);
  
  double lengthRatio=0.01;
  bool isInteger = false;  //do not do translational seamless.
  std::cout<<"Solving rotationally-seamless parameterization"<<std::endl;
  directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, cutUVRot);
  std::cout<<"Done!"<<std::endl;
  
  isInteger = true;  //do not do translational seamless.
  std::cout<<"Solving fully-seamless parameterization"<<std::endl;
  directional::parameterize(VMeshWhole, FMeshWhole, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, cutUVFull);
  std::cout<<"Done!"<<std::endl;
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0,
  0.5,0.0,1.0;
  
  //apending and updating raw field mesh
  viewer.append_mesh();
  
  //singularity mesh
  viewer.append_mesh();
  
  //seams mesh
  viewer.append_mesh();
  
  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.data_list[0].show_lines=false;
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


