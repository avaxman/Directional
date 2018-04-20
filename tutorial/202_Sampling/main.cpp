#include <igl/viewer/Viewer.h>
#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <directional/dual_cycles.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/index_prescription.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/euler_characteristic.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/power_to_representative.h>
#include <directional/power_field.h>
#include <directional/singularity_spheres.h>
#include <directional/glyph_lines_raw.h>

Eigen::VectorXi singVertices;
Eigen::VectorXi singIndices;

Eigen::SparseMatrix<double> basisCycles;
Eigen::VectorXi prinSingIndices;
Eigen::RowVector3d rawGlyphColor;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, BC, FN;
Eigen::VectorXi matching;
Eigen::MatrixXd rawField;
Eigen::VectorXd effort;
Eigen::VectorXi constFaces;
Eigen::MatrixXd constVecMat;

Eigen::VectorXd cycleCurvature;
Eigen::VectorXi vertex2cycle;
Eigen::VectorXi innerEdges;

int N=2;
double vfScale=0.01;

bool showSmoothness=false;
double globalRotation=0.0;

typedef enum {TRIVIAL_ONE_SING, TRIVIAL_PRINCIPAL_MATCHING, IMPLICIT_FIELD} ViewingModes;
ViewingModes viewingMode=TRIVIAL_ONE_SING;

igl::viewer::Viewer viewer;


void update_directional_field()
{
  
  using namespace Eigen;
  using namespace std;
  VectorXd rotationAngles;
  prinSingIndices=VectorXi::Zero(basisCycles.rows());
  for (int i=0;i<singVertices.size();i++)
    prinSingIndices(singVertices[i])=singIndices[i];
  
  double TCError;
  directional::index_prescription(V,F,innerEdges, basisCycles,cycleCurvature,prinSingIndices,N,rotationAngles, TCError);
  
  Eigen::MatrixXd representative;
  directional::rotation_to_representative(V, F,EV,EF,rotationAngles,N,globalRotation, representative);
  directional::representative_to_raw(V,F,representative,N, rawField);
  
  if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
    Eigen::VectorXd effort;
    directional::principal_matching(V, F,EV, EF, FE, rawField, matching, effort);
    directional::effort_to_indices(V,F,EV, EF, effort,N,prinSingIndices);
  }
  
  if (viewingMode==IMPLICIT_FIELD){
    constVecMat.conservativeResize(constFaces.rows(),3);
    for (int i=0;i<constFaces.size();i++)
      constVecMat.row(i)<<rawField.block(constFaces(i),0,1,3).normalized();
    
    Eigen::VectorXd effort;
    Eigen::MatrixXcd powerField;
    directional::power_field(V, F, constFaces, constVecMat, N, powerField);
    directional::power_to_representative(V,F, powerField,N,representative);
    representative.rowwise().normalize();
    directional::representative_to_raw(V,F,representative,N, rawField);
    directional::principal_matching(V, F,EV,EF,FE,rawField, matching, effort);
    directional::effort_to_indices(V,F,EV,EF,effort,N,prinSingIndices);
  }
}


void update_mesh()
{
  
  Eigen::MatrixXd C(F.rows(),3);
  C.col(0)=Eigen::VectorXd::Constant(F.rows(),1.0);
  C.col(1)=Eigen::VectorXd::Constant(F.rows(),1.0);
  C.col(2)=Eigen::VectorXd::Constant(F.rows(),1.0);
  
  for (int i = 0; i < constFaces.rows(); i++)
    C.row(constFaces(i)) = Eigen::RowVector3d(1, 0, 0);
  
  Eigen::MatrixXd fullV=V;
  Eigen::MatrixXi fullF=F;
  Eigen::MatrixXd fullC=C;
  
  Eigen::VectorXi currIndices;
  
  if (viewingMode==TRIVIAL_ONE_SING)
      directional::singularity_spheres(V, F, singVertices, singIndices, directional::defaultSingularityColors(N), false, true, fullV, fullF, fullC);
  
  if ((viewingMode==TRIVIAL_PRINCIPAL_MATCHING)||(viewingMode==IMPLICIT_FIELD))
    directional::singularity_spheres(V, F, prinSingIndices, directional::defaultSingularityColors(N), false, true, fullV, fullF, fullC);
  
  directional::glyph_lines_raw(V, F, rawField, rawGlyphColor, false, true, fullV, fullF, fullC);
  
  viewer.data.clear();
  viewer.data.set_face_based(true);
  viewer.data.set_mesh(fullV, fullF);
  viewer.data.set_colors(fullC);
}


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifiers)
{
  using namespace std;
  switch(key)
  {
    case '1': viewingMode=TRIVIAL_ONE_SING;
      break;
    case '2': viewingMode=TRIVIAL_PRINCIPAL_MATCHING;
      break;
    case '3': viewingMode=IMPLICIT_FIELD;
      break;
      
    case 'A':{
      singIndices[0]++;
      singIndices[1]--;
      cout<<"singularity index: "<<singIndices[0]<<std::endl;
      break;
    }
    case 'S':{
      singIndices[0]--;
      singIndices[1]++;
      cout<<"singularity index: "<<singIndices[0]<<std::endl;
      break;
    }
      
    case 'D':{
      globalRotation+=igl::PI/32;
      std::cout<<"globalRotation" <<globalRotation<<std::endl;
      break;
    }
      
    default: break;  //dunno why this is needed but it helps...
      
  }
  update_directional_field();
  update_mesh();
  return true;
}

double sign(double x){
  if (x>0) return 1.0;
  if (x<0) return -1.0;
  return 0.0;
}



int main()
{
  using namespace Eigen;
  using namespace std;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/spherers.obj", V, F);
  igl::edge_topology(V, F, EV,FE,EF);
  igl::barycenter(V,F,BC);
  igl::per_face_normals(V,F,FN);
  
  directional::dual_cycles(V, F,EV, EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
  
  //taking midway faces as constraints for the implicit field interpolation
  vector<int> constFacesList;
  for (int i=0;i<F.rows();i++){
    for (int j=0;j<3;j++)
      if (sign(V.row(F(i,j))(2))!=sign(V.row(F(i,(j+1)%3))(2))){
        constFacesList.push_back(i);
        break;
      }
  }
  constFaces.resize(constFacesList.size());
  for (int i=0;i<constFacesList.size();i++)
    constFaces(i)=constFacesList[i];
  
  rawGlyphColor <<0.0, 0.2, 1.0;
  
  singVertices.resize(2);
  singIndices.resize(2);
  singVertices(0)=35;
  singVertices(1)=36;
  singIndices(0)=N;
  singIndices(1)=N;
  update_directional_field();
  update_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}
