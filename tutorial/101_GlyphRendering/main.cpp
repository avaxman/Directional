#include <directional/directional_viewer.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>

int N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field;
directional::DirectionalViewer viewer;

int main()
{
  directional::readOFF(TUTORIAL_DATA_PATH "/bumpy.off",mesh);
  ftb.init(mesh);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/bumpy.rawfield", ftb, N, field);
  directional::read_singularities(TUTORIAL_DATA_PATH "/bumpy.sings", field);
  directional::DirectionalViewer viewer;
  viewer.init();
  
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  //viewer.toggle_mesh_edges(false);

  viewer.launch();
}

