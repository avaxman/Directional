#include <vector>
#include <cstdlib>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/readDMAT.h>
#include <igl/readOBJ.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/conjugate_frame_fields.h>
#include <directional/ConjugateFFSolverData.h>
#include <directional/directional_viewer.h>
#include "tutorial_shared_path.h"


Eigen::VectorXi matching, indices;
Eigen::VectorXd effort;
Eigen::MatrixXi F, EV, EF, FE;
Eigen::MatrixXd V, CMesh;
Eigen::MatrixXd rawField,representative, barycenters;
Eigen::MatrixXcd pvField;
directional::DirectionalViewer viewer;

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
    CMesh=directional::default_mesh_color().replicate(F.rows(),1);
    for (int i = 0; i < b.rows(); i++)
      CMesh.row(b(i)) = directional::selected_face_color();
  }else{
    Eigen::VectorXd currConjugacy = (viewingMode==ORIGINAL_CONJUGACY ? conjugacyOrig: conjugacyConjugate);
    igl::jet(currConjugacy, 0.0,conjMaxOrig, CMesh);
  }
  viewer.set_mesh_colors(CMesh);
}


void update_raw_field_mesh()
{
  using namespace std;
  using namespace Eigen;
  
  if ((viewingMode==ORIGINAL_CONJUGACY) || (viewingMode==OPTIMIZED_CONJUGACY)){
    viewer.toggle_field();
    viewer.toggle_singularities();
    viewer.toggle_seams();
  } else {
    viewer.set_field(viewingMode==ORIGINAL_FIELD ? rawFieldOrig : rawFieldConjugate, directional::indexed_glyph_colors(viewingMode==ORIGINAL_FIELD ? rawFieldOrig : rawFieldConjugate));
    viewer.set_singularities(N, viewingMode==ORIGINAL_FIELD ? singVerticesOrig : singVerticesConj, viewingMode==ORIGINAL_FIELD ? singIndicesOrig : singIndicesConj);
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
  igl::readOBJ(TUTORIAL_SHARED_PATH "/inspired_mesh.obj", V, F);
  igl::edge_topology(V, F, EV, FE, EF);
  
  // Compute face barycenters
  igl::barycenter(V, F, barycenters);
  
  // Load constraints
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_b.dmat",b);
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_bc.dmat",bc);
  
  bc.conservativeResize(bc.rows(), 2*bc.cols());
  bc.block(0,6,bc.rows(),6) = -bc.block(0,0,bc.rows(),6);
  
  //initial solution
  Eigen::MatrixXcd pvField;
  directional::polyvector_field(V, F, b, bc, N, pvField);
  directional::polyvector_to_raw(V, F, pvField, N, rawFieldOrig);
  
  Eigen::VectorXi prinIndices;
  directional::principal_matching(V, F,EV, EF, FE, rawFieldOrig, matching, effort,singVerticesOrig, singIndicesOrig);
  
  // Initialize conjugate field with smooth field
  directional::ConjugateFFSolverData csdata(V,F);
  
  directional::conjugate_frame_fields(csdata, b, rawFieldOrig, rawFieldConjugate);
  
  directional::principal_matching(V, F,EV, EF, FE, rawFieldConjugate, matching, effort,singVerticesConj, singIndicesConj);
  
  csdata.evaluateConjugacy(rawFieldOrig, conjugacyOrig);
  conjMaxOrig = conjugacyOrig.lpNorm<Infinity>();
  
  csdata.evaluateConjugacy(rawFieldConjugate, conjugacyConjugate);
  
  //triangle mesh setup
  viewer.set_mesh(V, F);
  viewer.data().show_lines=false;
  update_triangle_mesh();
  update_raw_field_mesh();

  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  
}
