#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/principal_matching.h>
#include <directional/conjugate_frame_fields.h>
#include <directional/ConjugateFFSolverData.h>
#include <vector>
#include <cstdlib>

#include "tutorial_shared_path.h"


Eigen::VectorXi matching, indices;
Eigen::VectorXd effort;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField,representative, cValues, barycenters;
Eigen::MatrixXcd pvField;
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd glyphPrincipalColors(5,3);

Eigen::MatrixXd rawFieldOrig, rawFieldConjugate;
Eigen::VectorXd conjugacyOrig, conjugacyConjugate;
Eigen::VectorXi singVerticesOrig, singVerticesConj;
Eigen::VectorXi singIndicesOrig, singIndicesConj;

int N=4;
double conjMax;


// Face barycenters
Eigen::MatrixXd B;

// Scale for visualizing the fields
double global_scale;

// Input constraints
Eigen::VectorXi b;
Eigen::MatrixXd bc;



typedef enum {ORIGINAL_FIELD, ORIGINAL_CONJUGACY, OPTIMIZED_FIELD, OPTIMIZED_CONJUGACY} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;


void update_triangle_mesh()
{
  if ((viewingMode ==ORIGINAL_FIELD)||(viewingMode ==OPTIMIZED_FIELD)){
    CMesh = Eigen::MatrixXd::Constant(FMesh.rows(), 3, 1.0);
    for (int i = 0; i < b.rows(); i++)
      CMesh.row(b(i)) = Eigen::RowVector3d(0.5,0.1,0.1);
  }else{  
    Eigen::VectorXd currConjugacy = (viewingMode==ORIGINAL_CONJUGACY ? conjugacyOrig: conjugacyConjugate);
    igl::jet(currConjugacy, 0.0,conjMax, CMesh);
    
  }
  viewer.data_list[0].set_colors(CMesh);
}


void update_raw_field_mesh()
{
  using namespace std;
  using namespace Eigen;
  
  if ((viewingMode==ORIGINAL_CONJUGACY) || (viewingMode==OPTIMIZED_CONJUGACY)){
    for (int i=1;i<=2;i++){  //hide all other meshes
      viewer.data_list[i].show_faces=false;
      viewer.data_list[i].show_lines = false;
    }
  } else {
    
    Eigen::MatrixXd fullGlyphColor(FMesh.rows(),3*N);
    for (int i=0;i<FMesh.rows();i++)
      for (int j=0;j<N;j++)
        fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
    
    
    directional::glyph_lines_raw(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? rawFieldOrig : rawFieldConjugate), fullGlyphColor,VField, FField, CField);
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    //singularity mesh
    directional::singularity_spheres(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? singVerticesOrig : singVerticesConj), (viewingMode==ORIGINAL_FIELD ? singIndicesOrig : singIndicesConj), directional::defaultSingularityColors(N), VSings, FSings, CSings);
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings, FSings);
    viewer.data_list[2].set_colors(CSings);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
  }
  
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  
  if ((key >= '1') && (key <='4'))
    viewingMode = (ViewingModes)(key - '1');
  
  update_triangle_mesh();
  update_raw_field_mesh();
  return false;
}


int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  std::cout <<
  "  1      Original field" << std::endl <<
  "  2      Conjugacy of original field" << std::endl <<
  "  3      Conjugate field" << std::endl <<
  "  4      Conjugacy of optimized field." << std::endl;
  
  // Load a mesh in OBJ format
  igl::readOBJ(TUTORIAL_SHARED_PATH "/inspired_mesh.obj", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  // Compute face barycenters
  igl::barycenter(VMesh, FMesh, barycenters);
  
  // Load constraints
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_b.dmat",b);
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_bc.dmat",bc);
  
  bc.conservativeResize(bc.rows(), 2*bc.cols());
  bc.block(0,6,bc.rows(),6) = -bc.block(0,0,bc.rows(),6);
  
  //initial solution
  Eigen::MatrixXcd pvField;
  directional::polyvector_field(VMesh, FMesh, b, bc, N, pvField);
  directional::polyvector_to_raw(VMesh, FMesh, pvField, N, rawFieldOrig);
  
  //cout<<"rawFieldOrig: "<<bc<<endl;
  
  Eigen::VectorXi prinIndices;
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawFieldOrig, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N,prinIndices);
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<VMesh.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVerticesOrig.resize(singVerticesList.size());
  singIndicesOrig.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVerticesOrig(i)=singVerticesList[i];
    singIndicesOrig(i)=singIndicesList[i];
  }
  
  // Initialize conjugate field with smooth field
  directional::ConjugateFFSolverData csdata(VMesh,FMesh);
  
  double lambdaOut = directional::conjugate_frame_fields(csdata, b, rawFieldOrig, rawFieldConjugate);
  
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawFieldConjugate, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N,prinIndices);
  
  singVerticesList.clear();
  singIndicesList.clear();
  for (int i=0;i<VMesh.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVerticesConj.resize(singVerticesList.size());
  singIndicesConj.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVerticesConj(i)=singVerticesList[i];
    singIndicesConj(i)=singIndicesList[i];
  }
  
  csdata.evaluateConjugacy(rawFieldOrig, conjugacyOrig);
  conjMax = conjugacyOrig.lpNorm<Infinity>();
  
  csdata.evaluateConjugacy(rawFieldConjugate, conjugacyConjugate);
  
  //triangle mesh setup
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(Eigen::RowVector3d::Constant(3,1.0));
  viewer.data().show_lines = false;
  
  //apending and updating raw field mesh
  viewer.append_mesh();
  
  //singularity mesh
  viewer.append_mesh();
  
  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.selected_data_index=0;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  
}
