#include <directional/directional_viewer.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>

int N;
directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField field;
directional::DirectionalViewer viewer;

int main()
{
  directional::readOFF(TUTORIAL_DATA_PATH "/bumpy.off",mesh);
  ftb.init(mesh);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/bumpy.rawfield", ftb, N, field);
  directional::read_singularities(TUTORIAL_DATA_PATH "/bumpy.sings", field);
  
  viewer.init();
  viewer.set_surface_mesh(mesh);
  viewer.set_cartesian_field(field);
  viewer.launch();
}

