#include <igl/viewer/Viewer.h>
#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <directional/dual_cycles.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/get_indices.h>
#include <directional/trivial_connection.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/euler_characteristic.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/complex_to_representative.h>
#include <directional/complex_field.h>

std::vector<int> singVertices;
std::vector<int> singIndices;

Eigen::VectorXi prinSingIndices;
Eigen::MatrixXi meshF, EV, FE, EF;
Eigen::MatrixXd meshV, BC, FN;
Eigen::SparseMatrix<double, Eigen::RowMajor> basisCycleMat;
Eigen::VectorXi indices, matching;
Eigen::MatrixXd directionalField;
Eigen::VectorXd effort;
Eigen::VectorXi constFaces;
Eigen::MatrixXd constVecMat;
Eigen::MatrixXd positiveColors(4,3), negativeColors(4,3);

int N=2;
double vfScale=0.01;



bool showSmoothness=false;
double globalRotation=0.0;

typedef enum {TRIVIAL_ONE_SING, TRIVIAL_PRINCIPAL_MATCHING, IMPLICIT_FIELD} ViewingModes;
ViewingModes viewingMode=TRIVIAL_ONE_SING;

igl::viewer::Viewer viewer;

void ConcatMeshes(const Eigen::MatrixXd &VA, const Eigen::MatrixXi &FA, const Eigen::MatrixXd &VB, const Eigen::MatrixXi &FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F)
{
    V.resize(VA.rows() + VB.rows(), VA.cols());
    V << VA, VB;
    F.resize(FA.rows() + FB.rows(), FA.cols());
    F << FA, (FB.array() + VA.rows());
}



void UpdateDirectionalField()
{
    
    using namespace Eigen;
    using namespace std;
    typedef complex<double> Complex;
    VectorXd rotationAngles;
    indices=VectorXi::Zero(basisCycleMat.rows());
    for (int i=0;i<singVertices.size();i++)
        indices(singVertices[i])=singIndices[i];
    
    if (indices.sum()!=N*igl::euler_characteristic(meshV, meshF)){
        std::cout<<"Warning! the prescribed singularities are not compatible with topology."<<std::endl;
        std::cout<<"chi = "<<igl::euler_characteristic(meshV, meshF)<<", indices.sum()="<<indices.sum()<<std::endl;
        return;
    }
    
    double TCError;
    directional::trivial_connection(meshV,meshF,basisCycleMat,indices,N,rotationAngles, TCError);
    cout<<"Trivial connection error: "<<TCError<<std::endl;
    
    Eigen::MatrixXd representative;
    
    directional::rotation_to_representative(meshV, meshF,EV,EF,rotationAngles,N,globalRotation, representative);
    directional::representative_to_raw(meshV,meshF,representative,N, directionalField);
    
    
    if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
        Eigen::VectorXd effort;
        directional::principal_matching(meshV, meshF,directionalField,N, effort);
        directional::get_indices(meshV,meshF,basisCycleMat,effort,N,prinSingIndices);
        std::cout<<"prinSingIndices sum: "<<prinSingIndices.sum()<<std::endl;
    }

    //overriding current field
    if (viewingMode==IMPLICIT_FIELD){
        constVecMat.conservativeResize(constFaces.rows(),3);
        for (int i=0;i<constFaces.size();i++)
            constVecMat.row(i)<<directionalField.block(constFaces(i),0,1,3).normalized();
        
        Eigen::VectorXd effort;
        Eigen::MatrixXcd complexField;
        directional::complex_field(meshV, meshF, constFaces, constVecMat, N, complexField);
        directional::complex_to_representative(meshV,meshF, complexField,N,representative);
        representative.rowwise().normalize();
        directional::representative_to_raw(meshV,meshF,representative,N, directionalField);
        directional::principal_matching(meshV, meshF,directionalField,N, effort);
        directional::get_indices(meshV,meshF,basisCycleMat,effort,N,prinSingIndices);
    }
}




void UpdateCurrentView()
{
    using namespace Eigen;
    using namespace std;
    MatrixXd singPoints;
    MatrixXd singColors;
    VectorXi currIndices;
    if (viewingMode==TRIVIAL_ONE_SING)
        currIndices=indices;
    
    if ((viewingMode==TRIVIAL_PRINCIPAL_MATCHING)||(viewingMode==IMPLICIT_FIELD))
        currIndices=prinSingIndices;
    
    
    std::cout<<"starting to draw"<<std::endl;
    
    Eigen::MatrixXd fieldV, fieldC;
    Eigen::MatrixXi fieldF;
    directional::drawable_field(meshV, meshF, directionalField, Eigen::RowVector3d(0,0,1), N, directional::field_draw_flags::NONE, fieldV, fieldF, fieldC);
    Eigen::MatrixXd meshC = Eigen::RowVector3d(1, 1, 1).replicate(meshF.rows(), 1);
    
    for (int i = 0; i < constFaces.rows(); i++)
        meshC.row(constFaces(i)) = Eigen::RowVector3d(1, 0, 0);
    
    std::cout<<"created drawable field"<<std::endl;
    
    
    Eigen::MatrixXd singV, singC;
    Eigen::MatrixXi singF;
    
    directional::draw_singularities(meshV, currIndices, positiveColors, negativeColors, .01, singV, singF, singC);
    
    std::cout<<"created drawable singularities"<<std::endl;
    
    Eigen::MatrixXd a, V;
    Eigen::MatrixXi b, F;
    ConcatMeshes(meshV, meshF, fieldV, fieldF, a, b);
    ConcatMeshes(a, b, singV, singF, V, F);
    Eigen::MatrixXd C(F.rows(), 3);
    C<< meshC, fieldC, singC;
    
    // Update the viewer
    viewer.data.clear();
    viewer.data.set_face_based(true);
    viewer.data.set_mesh(V, F);
    viewer.data.set_colors(C);

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
    UpdateDirectionalField();
    UpdateCurrentView();
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
    igl::readOBJ(TUTORIAL_SHARED_PATH "/spherers.obj", meshV, meshF);
    igl::edge_topology(meshV, meshF, EV,FE,EF);
    igl::barycenter(meshV,meshF,BC);
    igl::per_face_normals(meshV,meshF,FN);
    
    VectorXi primalTreeEdges, dualTreeEdges;
    directional::dual_cycles(meshF,EV, EF, basisCycleMat);
  
    //taking midway faces as constraints for the implicit field interpolation
    vector<int> constFacesList;
    for (int i=0;i<meshF.rows();i++){
        for (int j=0;j<3;j++)
            if (sign(meshV.row(meshF(i,j))(2))!=sign(meshV.row(meshF(i,(j+1)%3))(2))){
                constFacesList.push_back(i);
                break;
            }
    }
    constFaces.resize(constFacesList.size());
    for (int i=0;i<constFacesList.size();i++)
        constFaces(i)=constFacesList[i];
    
    //cout<<"constFaces: "<<constFaces<<endl;
    
    positiveColors << 0, 0,0.25,
                      0, 0, 0.5,
                      0, 0,0.75,
                      0, 0,1.0;
    
    negativeColors << 0.25, 0, 0,
                      0.5, 0, 0,
                      0.75, 0, 0,
                      1.0, 0, 0;
    
    singVertices.resize(2);
    singIndices.resize(2);
    singVertices[0]=35;
    singVertices[1]=36;
    singIndices[0]=N;
    singIndices[1]=N;
    UpdateDirectionalField();
    UpdateCurrentView();
    viewer.callback_key_down = &key_down;
    viewer.launch();
}
