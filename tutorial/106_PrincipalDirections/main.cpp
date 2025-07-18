#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/directional_viewer.h>
#include <directional/IntrinsicVertexTangentBundle.h>


directional::TriMesh mesh;
directional::IntrinsicVertexTangentBundle vtb;
directional::DirectionalViewer viewer;

int main()
{
    // Load mesh
    directional::readOFF(TUTORIAL_DATA_PATH "/botanic-garden-bubble.off", mesh);
    vtb.init(mesh);
    
    Eigen::MatrixXd minField(mesh.V.rows(), 6);
    minField<<mesh.minVertexPrincipalDirections, -mesh.minVertexPrincipalDirections;
    Eigen::MatrixXd maxField(mesh.V.rows(), 6);
    maxField<<mesh.maxVertexPrincipalDirections, -mesh.maxVertexPrincipalDirections;
     
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_raw_field(mesh.V, minField, "Min Principal Direction",  0, 0.3*mesh.avgEdgeLength);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0/255.0}, 0);
    viewer.set_raw_field(mesh.V, maxField, "Max Principal Direction",  1, 0.3*mesh.avgEdgeLength);
    viewer.set_field_color({125.0/255.0, 107.0/255.0, 8.0/255.0}, 1);
    viewer.set_raw_field(mesh.V, mesh.vertexNormals, "Vertex normals",  2, 0.3*mesh.avgEdgeLength);
    viewer.set_field_color({8.0/255.0, 125.0/255.0, 107.0/255.0}, 2);
    
    viewer.set_surface_vertex_data(mesh.vertexPrincipalCurvatures.col(0), "Min curvature");
    viewer.set_surface_vertex_data(mesh.vertexPrincipalCurvatures.col(1), "Max curvature");
    viewer.launch();
}
