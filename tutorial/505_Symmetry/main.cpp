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


int N;
Eigen::MatrixXi FMeshWhole, FMeshCut, FField, FSings, FSeams;
Eigen::MatrixXd VMeshWhole, VMeshCut, VField, VSings, VSeams;
Eigen::MatrixXd CField, CSeams, CSings;
Eigen::MatrixXd rawField, combedField, barycenters;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd cutUVTri, cutUVSign, cornerWholeUV, cutReducedUV;
igl::opengl::glfw::Viewer viewer;

typedef enum {FIELD, SIGN_SYMMETRY, TRI_SYMMETRY} ViewingModes;
ViewingModes viewingMode=FIELD;


void trace_isolines(const Eigen::MatrixXd& paramFuncs,
                    Eigen::MatrixXd& P1,
                    Eigen::MatrixXd& P2,
                    Eigen::MatrixXd& C)
{
  
  Eigen::MatrixXd funcColors(8,3);
  funcColors<<1.0,0.0,0.0,
  0.0,1.0,0.0,
  0.0,0.0,1.0,
  1.0,0.0,0.5,
  0.5,1.0,0.0,
  0.0,0.5,1.0,
  1.0,0.5,0.0,
  0.0,1.0,0.5;
  funcColors.array()/=2.0;
  int jumps = (N%2 == 0 ? 2 : 1);
  P1.resize(0,3), P2.resize(0,3), C.resize(0,3);
  for (int i=0;i<paramFuncs.cols()/jumps;i++){
    Eigen::VectorXd d = paramFuncs.col(i);
    
    /*std::cout<<"d.min(): "<<d.minCoeff()<<std::endl;
     std::cout<<"d.max(): "<<d.maxCoeff()<<std::endl;*/
    Eigen::MatrixXd isoV;
    Eigen::MatrixXi isoE;
    igl::isolines(VMeshCut,FMeshCut, d, 100, isoV, isoE);
    
    int prevIndex = P1.rows();
    P1.conservativeResize(P1.rows()+isoE.rows(),3);
    P2.conservativeResize(P2.rows()+isoE.rows(),3);
    C.conservativeResize(C.rows()+isoE.rows(),3);
    for (int i=0;i<isoE.rows();i++){
      P1.row(prevIndex+i)=isoV.row(isoE(i,0));
      P2.row(prevIndex+i)=isoV.row(isoE(i,1));
    }
    
    C.block(prevIndex, 0, isoE.rows(),3)=funcColors.row(i).replicate(isoE.rows(),1);
    
  }
}



void update_triangle_mesh()
{
  if (viewingMode==FIELD){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshWhole, FMeshWhole);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_face_based(false);
    viewer.data_list[0].show_lines=false;
  } else if ((viewingMode==SIGN_SYMMETRY) || (viewingMode==TRI_SYMMETRY)){
    viewer.data_list[0].clear();
    viewer.data_list[0].set_mesh(VMeshCut, FMeshCut);
    viewer.data_list[0].set_colors(directional::default_mesh_color());
    viewer.data_list[0].set_face_based(false);
    viewer.data_list[0].show_lines=false;
  }
  
  viewer.data_list[4].show_faces=(viewingMode==SIGN_SYMMETRY);
  viewer.data_list[5].show_faces=(viewingMode==TRI_SYMMETRY);
}

void update_raw_field_mesh()
{
  for (int i=1;i<=4;i++)  //hide all other meshes
    viewer.data_list[i].show_faces=(viewingMode==FIELD);
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = SIGN_SYMMETRY; break;
    case '3': viewingMode = TRI_SYMMETRY; break;
  }
  update_triangle_mesh();
  update_raw_field_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show only sign-symmetric integrated functions" << std::endl <<
  "  3  Show triangular-symmetric integrated functions" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/teddy.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/teddy.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);
  igl::barycenter(VMeshWhole, FMeshWhole, barycenters);
  
  //combing and cutting
  directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMeshWhole,FMeshWhole,EV, EF, effort,matching, N,singVertices, singIndices);
  
  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  
  intData.verbose=true;
  intData.integralSeamless=false;
  intData.localInjectivity=false;
    
  std::cout<<"Free (sign-symmetric) Integrating..."<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut, cutReducedUV,  cutUVSign,cornerWholeUV);
  

  std::cout<<"Done!"<<std::endl;
  
  //Triangular symmetric  + sign symmetric for N=6
  intData.symmFunc.resize(6,2);
  intData.symmFunc.row(0)<<1,0;
  intData.symmFunc.row(1)<<0,1;
  intData.symmFunc.row(2)<<-1,1;
  intData.symmFunc.block(3,0,3,2)=-intData.symmFunc.block(0,0,3,2);
  std::cout<<"Solving triangular-constrained integration..."<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut, cutReducedUV,  cutUVTri,cornerWholeUV);
  std::cout<<"Done!"<<std::endl;
  
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
  
  //sign-symmetric isolines mesh

  viewer.append_mesh();
  Eigen::MatrixXd P1,P2,C;
  trace_isolines(cutUVTri, P1, P2, C);
  viewer.data_list[4].clear();
  Eigen::MatrixXd VSign, CSign;
  Eigen::MatrixXi TSign;
  directional::line_cylinders(P1, P2, 2.5, C, 6, VSign, TSign, CSign);
  viewer.data_list[4].set_mesh(VSign, TSign);
  viewer.data_list[4].set_colors(CSign);
  viewer.data_list[4].show_faces = true;
  viewer.data_list[4].show_lines = false;
  
  //tri-symmetric isolines mesh
  viewer.append_mesh();
  trace_isolines(cutUVSign, P1, P2, C);
  viewer.data_list[5].clear();
  Eigen::MatrixXd VTri, CTri;
  Eigen::MatrixXi TTri;
  directional::line_cylinders(P1, P2, 2.5, C, 6, VTri, TTri, CTri);
  viewer.data_list[5].set_mesh(VTri, TTri);
  viewer.data_list[5].set_colors(CTri);
  viewer.data_list[5].show_faces = true;
  viewer.data_list[5].show_lines = false;
  
  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.data_list[0].show_lines=false;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


