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
Eigen::MatrixXd rawField,representative, barycenters;
Eigen::MatrixXcd pvField;
igl::opengl::glfw::Viewer viewer;

Eigen::MatrixXd rawFieldOrig, rawFieldConjugate;
Eigen::VectorXd conjugacyOrig, conjugacyConjugate;
Eigen::VectorXi singVerticesOrig, singVerticesConj;
Eigen::VectorXi singIndicesOrig, singIndicesConj;

int N=4;
double conjMaxOrig;

// Input constraints
Eigen::VectorXi b;
Eigen::MatrixXd bc;

typedef enum {ORIGINAL_FIELD, ORIGINAL_CONJUGACY, OPTIMIZED_FIELD, OPTIMIZED_CONJUGACY} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;


void update_triangle_mesh()
{
  if ((viewingMode==ORIGINAL_FIELD)||(viewingMode==OPTIMIZED_FIELD)){
    CMesh=directional::default_mesh_color().replicate(FMesh.rows(),1);
    for (int i = 0; i < b.rows(); i++)
      CMesh.row(b(i)) = directional::selected_face_color();
  }else{
    Eigen::VectorXd currConjugacy = (viewingMode==ORIGINAL_CONJUGACY ? conjugacyOrig: conjugacyConjugate);
    igl::jet(currConjugacy, 0.0,conjMaxOrig, CMesh);
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
    
    directional::glyph_lines_raw(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? rawFieldOrig : rawFieldConjugate),
                                 directional::indexed_glyph_colors((viewingMode==ORIGINAL_FIELD ? rawFieldOrig : rawFieldConjugate)),VField, FField, CField, 1.5);
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    //singularity mesh
    directional::singularity_spheres(VMesh, FMesh,N,  (viewingMode==ORIGINAL_FIELD ? singVerticesOrig : singVerticesConj), (viewingMode==ORIGINAL_FIELD ? singIndicesOrig : singIndicesConj), VSings, FSings, CSings, 2.0);
    
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
  
  Eigen::VectorXi prinIndices;
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawFieldOrig, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N,singVerticesOrig, singIndicesOrig);
  
  // Initialize conjugate field with smooth field
  directional::ConjugateFFSolverData csdata(VMesh,FMesh);
  
  directional::conjugate_frame_fields(csdata, b, rawFieldOrig, rawFieldConjugate);
  
  directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawFieldConjugate, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N,singVerticesConj, singIndicesConj);
  
  csdata.evaluateConjugacy(rawFieldOrig, conjugacyOrig);
  conjMaxOrig = conjugacyOrig.lpNorm<Infinity>();
  
  csdata.evaluateConjugacy(rawFieldConjugate, conjugacyConjugate);
  
  //triangle mesh setup
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(directional::default_mesh_color());
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
