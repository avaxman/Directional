#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/line_cylinders.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/isolines.h>
#include <directional/is_order_preserving.h>
#include <directional/gradient.h>


int N;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams;
Eigen::MatrixXd CField, CSeams, CSings;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd cutUVReg, cutUVLI, cornerWholeUV, cutReducedUV;
Eigen::MatrixXd funcColors(6,3), VIsoLinesReg, CIsoLinesReg,VIsoLinesLI, CIsoLinesLI;
Eigen::MatrixXi FIsoLinesReg, FIsoLinesLI;
igl::opengl::glfw::Viewer viewer;

typedef enum {FIELD, INTEGRATION_REGULAR, INTEGRATION_LI} ViewingModes;
ViewingModes viewingMode=FIELD;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;

// Create a texture that hides the integer translation in the parametrization



void update_view()
{
  for (int i=1;i<=3;i++)  //hide all other meshes
    viewer.data_list[i].show_faces=(viewingMode==FIELD);
  
  viewer.data_list[4].show_faces=(viewingMode==INTEGRATION_REGULAR);
  viewer.data_list[5].show_faces=(viewingMode==INTEGRATION_LI);
}

void append_meshes(const Eigen::MatrixXd& VAdd, const Eigen::MatrixXi& FAdd, const Eigen::MatrixXd& CAdd, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& C){
  int oldVSize = V.rows();
  int oldFSize = F.rows();
  int oldCSize = C.rows();
  
  V.conservativeResize(V.rows()+VAdd.rows(),3);
  F.conservativeResize(F.rows()+FAdd.rows(),3);
  C.conservativeResize(C.rows()+CAdd.rows(),3);
  
  V.block(oldVSize, 0, VAdd.rows(),3)=VAdd;
  F.block(oldFSize, 0, FAdd.rows(),3)=FAdd.array()+oldVSize;
  C.block(oldCSize, 0, CAdd.rows(),3)=CAdd;
}

void create_isoline_mesh(const Eigen::MatrixXd& paramFuncsN,
                         Eigen::MatrixXd& VIsoLines,
                         Eigen::MatrixXi& FIsoLines,
                         Eigen::MatrixXd& CIsoLines)
{
  double l = 1.25*igl::avg_edge_length(VMeshWhole, FMeshWhole);
  double isolineRadius=0.02;
  int jumps = (N%2 == 0 ? 2 : 1);
  std::cout<<"paramFuncsN.rows(), paramFuncsN.cols(): "<<paramFuncsN.rows()<<","<<paramFuncsN.cols()<<std::endl;
  for (int i=0;i<paramFuncsN.cols()/jumps;i++){
    Eigen::VectorXd d = paramFuncsN.col(i);
    Eigen::MatrixXd isoV;
    Eigen::MatrixXi isoE;
    igl::isolines(VMeshCut,FMeshCut, d, 100, isoV, isoE);
    
    Eigen::MatrixXd P1(isoE.rows(),3), P2(isoE.rows(),3);
    for (int i=0;i<isoE.rows();i++){
      P1.row(i)=isoV.row(isoE(i,0));
      P2.row(i)=isoV.row(isoE(i,1));
    }
    
    Eigen::MatrixXd VIsoLinesTemp, CIsoLinesTemp;
    Eigen::MatrixXi FIsoLinesTemp;
    directional::line_cylinders(P1, P2, l*isolineRadius,funcColors.row(i).replicate(P1.rows(),1),4, VIsoLinesTemp, FIsoLinesTemp, CIsoLinesTemp);
    
    append_meshes(VIsoLinesTemp, FIsoLinesTemp, CIsoLinesTemp, VIsoLines, FIsoLines, CIsoLines);
    
  }
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = INTEGRATION_REGULAR; break;
    case '3': viewingMode = INTEGRATION_LI; break;
    /*case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-rot-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVRot, FMeshCut);
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-full-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVFull, FMeshCut);
      break;*/
  }
  update_view();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show regular integration" << std::endl <<
  "  3  Show locally-injective integration" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/botanic-garden-bubble.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/botanic-garden-bubble.rawfield", N, rawField);
  
  //normalize the field
  for (int i=0;i<rawField.rows();i++)
    for (int j=0;j<rawField.cols();j+=3)
      rawField.block(i,j,1,3).normalize();
  
  std::cout<<"N: "<<N<<std::endl;
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  //combing and cutting
  Eigen::VectorXd curlNorm;
  directional::curl_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort, curlNorm);
  std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;
  directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,singVertices, singIndices);
  

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  
  intData.verbose=true;
  intData.localInjectivity=false;
  intData.integralSeamless = true;
  intData.roundSeams=false;
  intData.lengthRatio=0.04;
  std::cout<<"Solving regular integration"<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut, cutReducedUV,  cutUVReg,cornerWholeUV);
  std::cout<<"Done!"<<std::endl;
  
  
  Eigen::MatrixXd gradientField;
  directional::gradient(VMeshWhole, FMeshWhole, N, cornerWholeUV, gradientField);
  
  Eigen::VectorXi isOrderPreserving;
  directional::is_order_preserving(VMeshWhole, FMeshWhole, gradientField, isOrderPreserving);
  
  intData.verbose=true;
  intData.localInjectivity=true;
  intData.integralSeamless = true;
  std::cout<<"Solving locally-injective integration"<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut, cutReducedUV,  cutUVLI,cornerWholeUV);
  std::cout<<"Done!"<<std::endl;
  
  
  //Mesh visualization
  Eigen::MatrixXd meshColors(FMeshWhole.rows(),3);
  meshColors=directional::default_mesh_color().replicate(FMeshWhole.rows(),1);
  for (int i=0;i<isOrderPreserving.size();i++)
    if (!isOrderPreserving(i))
      meshColors.row(i)<<1.0,0.0,0.0;
      
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
  viewer.data_list[0].set_colors(meshColors);
  viewer.data_list[0].set_face_based(true);
  viewer.data_list[0].show_texture=false;
  viewer.data_list[0].show_lines=false;
  
  //raw field mesh
  viewer.append_mesh();
  directional::glyph_lines_raw(VMeshWhole, FMeshWhole, combedField, directional::indexed_glyph_colors(combedField), VField, FField, CField,1.0);
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_faces = true;
  viewer.data_list[1].show_lines = false;
  
  //singularity mesh
  viewer.append_mesh();
  directional::singularity_spheres(VMeshWhole, FMeshWhole, N, singVertices, singIndices, VSings, FSings, CSings,2.5);
  
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].set_colors(CSings);
  viewer.data_list[2].show_faces = true;
  viewer.data_list[2].show_lines = false;
  
  //seams mesh
  viewer.append_mesh();
  
  Eigen::VectorXi isSeam=Eigen::VectorXi::Zero(EV.rows());
  for (int i=0;i<FE.rows();i++)
    for (int j=0;j<3;j++)
      if (intData.face2cut(i,j))
        isSeam(FE(i,j))=1;
  directional::seam_lines(VMeshWhole, FMeshWhole, EV, combedMatching, VSeams, FSeams, CSeams,2.5);
  
  viewer.data_list[3].clear();
  viewer.data_list[3].set_mesh(VSeams, FSeams);
  viewer.data_list[3].set_colors(CSeams);
  viewer.data_list[3].show_faces = true;
  viewer.data_list[3].show_lines = false;
  
  funcColors<<1.0,0.0,0.0,
  0.0,1.0,0.0,
  0.0,0.0,1.0,
  0.5,0.5,0.0,
  0.0,0.5,0.5,
  0.5,0.0,0.5;
  
  //creating isoline mesh
  create_isoline_mesh(cutUVReg, VIsoLinesReg, FIsoLinesReg, CIsoLinesReg);
  create_isoline_mesh(cutUVLI, VIsoLinesLI, FIsoLinesLI, CIsoLinesLI);
  

  //Isolines of regular parameterization
  viewer.append_mesh();
  viewer.data_list[4].clear();
  viewer.data_list[4].set_mesh(VIsoLinesReg, FIsoLinesReg);
  viewer.data_list[4].set_colors(CIsoLinesReg);
  viewer.data_list[4].show_faces = false;
  viewer.data_list[4].show_lines = false;
  
  //Isolines of regular parameterization
  viewer.append_mesh();
  viewer.data_list[5].clear();
  viewer.data_list[5].set_mesh(VIsoLinesLI, FIsoLinesLI);
  viewer.data_list[5].set_colors(CIsoLinesLI);
  viewer.data_list[5].show_faces = false;
  viewer.data_list[5].show_lines = false;
  
  
  update_view();
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


