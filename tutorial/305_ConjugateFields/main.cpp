#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/write_raw_field.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/polyvector_iteration_functions.h>
#include <directional/PolyVectorData.h>
#include <directional/polyvector_field.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/curl_matching.h>
#include <directional/extrinsic_intrinsic_matrices.h>

Eigen::VectorXi constVertices;
directional::TriMesh mesh;
directional::IntrinsicVertexTangentBundle vtb;
directional::CartesianField pvFieldConjugate, rawFieldConjugate,constraintsField, rawFieldOrig, pvFieldOrig;
Eigen::MatrixXd constVectors;
double smoothWeight, roSyWeight, alignWeight;
directional::DirectionalViewer viewer;
int N = 4;



int main()
{
    // Load mesh
    directional::readOFF(TUTORIAL_DATA_PATH "/botanic-garden-bubble.off", mesh);
    vtb.init(mesh);
    pvFieldConjugate.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    pvFieldOrig.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    
    //discovering and constraining sharp edges
    std::vector<int> constVerticeslist;
    std::vector<Eigen::Vector3d> constVectorslist;
    for (int i=0;i<mesh.boundEdges.rows();i++){
        int vertexIndex = (mesh.EF(mesh.boundEdges(i),0)==-1 ? mesh.EV(mesh.boundEdges(i),1) : mesh.EV(mesh.boundEdges(i),0));
        constVerticeslist.push_back(vertexIndex);
        constVectorslist.push_back((mesh.V.row(mesh.EV(mesh.boundEdges(i),0))-mesh.V.row(mesh.EV(mesh.boundEdges(i),1))).normalized());
    }
    
    constVertices.resize(constVerticeslist.size());
    constVectors.resize(constVectorslist.size(),3);
    Eigen::MatrixXd constSources(constVectorslist.size(),3);  //For visualization
    for (int i=0;i<constVertices.size();i++){
        constVertices(i)=constVerticeslist[i];
        constVectors.row(i)=constVectorslist[i];
        constSources.row(i) = mesh.V.row(constVertices(i));
    }
    
    smoothWeight = 1.0;
    roSyWeight = 10.0;
    alignWeight = 0.1;
    
    directional::PolyVectorData pvData;
    pvData.N = N;
    pvData.tb = &vtb;
    pvData.verbose = true;
    pvData.constSpaces = constVertices;
    pvData.constVectors = constVectors;
    pvData.wAlignment = alignWeight*Eigen::VectorXd::Constant(constVertices.size(),1.0);
    pvData.wSmooth = smoothWeight;
    pvData.wRoSy = roSyWeight;
    
    //Computing regular PolyuVector field without iterations
    directional::polyvector_field(pvData, pvFieldOrig);
    directional::polyvector_to_raw(pvFieldOrig, rawFieldOrig, N%2==0);
    directional::principal_matching(rawFieldOrig);
    
    //Iterating for a conjugate field
    pvData.numIterations = 25;
    pvData.initImplicitFactor = 1.0;
    std::vector<directional::PvIterationFunction> iterationFunctions;
    iterationFunctions.push_back(directional::conjugate);
    directional::polyvector_field(pvData, pvFieldConjugate, iterationFunctions);
    directional::polyvector_to_raw(pvFieldConjugate, rawFieldConjugate, N%2==0);
    directional::principal_matching(rawFieldConjugate);
    
    //Visualization
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_raw_field(constSources, constVectors, mesh.avgEdgeLength, "Constraints",  0);
    viewer.set_field_color(directional::DirectionalViewer::default_vector_constraint_color());
    //viewer.highlight_vertices(constVertices,"Const Vertices", 0);
    //viewer.toggle_ver(true,0);
    viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
   
    //directional::CartesianField prinField;
   // prinField.init(ftb, directional::fieldTypeEnum::RAW_FIELD,N);
    Eigen::MatrixXd extField(mesh.V.rows(), 3*N);
    extField<<mesh.minVertexPrincipalDirections, mesh.maxVertexPrincipalDirections, -mesh.minVertexPrincipalDirections, -mesh.maxVertexPrincipalDirections;
   // prinField.set_extrinsic_field(extField);
    //testing conjugate function
    //directional::CartesianField rawField, conjugatePVField, prinPVField;
    //prinPVField.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD, N);
    //rawField.init(vtb, directional::fieldTypeEnum::RAW_FIELD, N);
    //rawField.set_extrinsic_field(extField);
    //for (int i=0;i<mesh.V.rows();i++)
    //    std::cout<<"conjugacy of raw field before conversion: "<<rawField.extField.row(i).head(3)*rawField.extField.row(i).segment(3,3).transpose()<<std::endl;
    //directional::raw_to_polyvector(rawField, prinPVField);
    //std::cout<<"prinPVField: "<<prinPVField.intField<<std::endl;
    //directional::polyvector_to_raw(prinPVField, rawField, N%2==0);
    //for (int i=0;i<mesh.V.rows();i++)
    //    std::cout<<"conjugacy of raw field after conversion: "<<rawField.extField.row(i).head(3)*mesh.Sv[i]*rawField.extField.row(i).segment(3,3).transpose()<<std::endl;
    //conjugatePVField = directional::conjugate(prinPVField, pvData);
    //directional::polyvector_to_raw(conjugatePVField, rawFieldConjugate, N%2==0);
    viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 2);
    viewer.set_raw_field(mesh.V, extField,1.0, "Principal directions Field", 2);
    viewer.set_vertex_data(mesh.vertexPrincipalCurvatures.col(0), mesh.vertexPrincipalCurvatures.col(0).minCoeff(), mesh.vertexPrincipalCurvatures.col(0).maxCoeff(), "Min curvature", 0);
    viewer.set_vertex_data(mesh.vertexPrincipalCurvatures.col(1), mesh.vertexPrincipalCurvatures.col(1).minCoeff(), mesh.vertexPrincipalCurvatures.col(1).maxCoeff(), "Max curvature", 0);
    //viewer.set_face_data((mesh.facePrincipalCurvatures.col(1)-mesh.facePrincipalCurvatures.col(0)).cwiseAbs(), mesh.facePrincipalCurvatures.col(0).minCoeff(), mesh.facePrincipalCurvatures.col(1).maxCoeff(), "umbilic", 0);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 2);
    viewer.launch();
}
