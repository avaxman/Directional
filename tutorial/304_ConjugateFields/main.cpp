#include <vector>
#include <cstdlib>
#include <igl/readDMAT.h>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/polyvector_field.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/conjugate_frame_fields.h>
#include <directional/ConjugateFFSolverData.h>
#include <directional/directional_viewer.h>
#include "tutorial_shared_path.h"

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawFieldOrig, rawFieldConjugate;
directional::DirectionalViewer viewer;
Eigen::VectorXd conjugacyOrig, conjugacyConjugate;

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
    viewer.set_selected_faces(b);
  }else{
    Eigen::VectorXd currConjugacy = (viewingMode==ORIGINAL_CONJUGACY ? conjugacyOrig: conjugacyConjugate);
    viewer.set_face_data(currConjugacy,0.0,conjMaxOrig);
  }

}

void update_raw_field_mesh()
{
  using namespace std;
  using namespace Eigen;
  
  if ((viewingMode==ORIGINAL_CONJUGACY) || (viewingMode==OPTIMIZED_CONJUGACY)){
    viewer.toggle_field(false);
    viewer.toggle_singularities(false);
    viewer.toggle_seams(false);
  } else {
    viewer.set_field(viewingMode==ORIGINAL_FIELD ? rawFieldOrig : rawFieldConjugate, directional::DirectionalViewer::indexed_glyph_colors(viewingMode==ORIGINAL_FIELD ? rawFieldOrig.extField : rawFieldConjugate.extField));
    viewer.toggle_field(true);
    viewer.toggle_singularities(true);
    viewer.toggle_seams(true);
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
  directional::readOBJ(TUTORIAL_SHARED_PATH "/inspired_mesh.obj", mesh);
  ftb.init(mesh);
  
  // Load constraints
  Eigen::VectorXi bFull;
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_b.dmat",bFull);
  Eigen::MatrixXd bcFull;
  igl::readDMAT(TUTORIAL_SHARED_PATH "/inspired_mesh_bc.dmat",bcFull);
  
  bcFull.conservativeResize(bcFull.rows(), 3*N);
  bcFull.block(0,6,bcFull.rows(),6) = -bcFull.block(0,0,bcFull.rows(),6);
  
  //putting the constraints in single-vector matrices
  bc.resize(N*bcFull.rows(),3);
  b.resize(N*bFull.rows());
  for (int i=0;i<bFull.rows();i++){
    for (int j=0;j<N;j++){
      b(i*N+j)=bFull(i);
      bc.row(i*N+j)=bcFull.block(i,3*j,1,3);
    }
  }
  
  //initial solution
  directional::CartesianField pvField;
  directional::polyvector_field(ftb, b, bc, N, pvField);
  directional::polyvector_to_raw(pvField, rawFieldOrig);
  
  directional::principal_matching(rawFieldOrig);
  
  // Initialize conjugate field with smooth field
  directional::ConjugateFFSolverData csdata(mesh);
  
  directional::conjugate_frame_fields(csdata, b, rawFieldOrig, rawFieldConjugate);
  
  directional::principal_matching(rawFieldConjugate);
  
  csdata.evaluateConjugacy(rawFieldOrig, conjugacyOrig);
  conjMaxOrig = conjugacyOrig.lpNorm<Infinity>();
  
  csdata.evaluateConjugacy(rawFieldConjugate, conjugacyConjugate);
  
  //triangle mesh setup
  viewer.set_mesh(mesh);
  viewer.set_field(rawFieldOrig);
  update_triangle_mesh();
  update_raw_field_mesh();

  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  
}
