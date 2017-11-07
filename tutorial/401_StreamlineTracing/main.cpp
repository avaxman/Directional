#include <igl/barycenter.h>
#include <igl/edge_topology.h>
#include <igl/local_basis.h>
#include <igl/parula.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/polyvector_field_matchings.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/sort_vectors_ccw.h>
#include <igl/streamlines.h>
//#include <igl/copyleft/comiso/nrosy.h>
#include <igl/viewer/Viewer.h>
#include <directional/complex_field.h>
#include <directional/complex_to_raw.h>

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>


// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;

Eigen::VectorXi cIDs;
Eigen::MatrixXd cValues;
Eigen::MatrixXcd complexField;
Eigen::MatrixXd raw;

igl::StreamlineData sl_data;
igl::StreamlineState sl_state;

int degree;         // degree of the vector field
int half_degree;    // degree/2 if treat_as_symmetric
bool treat_as_symmetric = true;

int anim_t = 0;
int anim_t_dir = 1;


bool pre_draw(igl::viewer::Viewer &viewer)
{
    using namespace Eigen;
    using namespace std;

    if (!viewer.core.is_animating)
        return false;

    igl::streamlines_next(V, F, sl_data, sl_state);
    Eigen::RowVector3d color = Eigen::RowVector3d::Zero();
    double value = ((anim_t) % 100) / 100.;

    if (value > 0.5)
        value = 1 - value;
    value = value / 0.5;
    igl::parula(value, color[0], color[1], color[2]);

    viewer.data.add_edges(sl_state.start_point, sl_state.end_point, color);

    anim_t += anim_t_dir;

    return false;
}

bool key_down(igl::viewer::Viewer &viewer, unsigned char key, int modifier)
{
    if (key == ' ')
    {
        viewer.core.is_animating = !viewer.core.is_animating;
        return true;
    }
    return false;
}

int main(int argc, char *argv[])
{
    using namespace Eigen;
    using namespace std;


    // Load a mesh in OFF format
    igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", V, F);
    // Create a Vector Field
    Eigen::VectorXi b;
    Eigen::MatrixXd bc;
    Eigen::VectorXd S; // unused

    b.resize(1);
    b << 0;
    bc.resize(1, 3);
    bc << 1, 1, 1;

    half_degree = 3;
    treat_as_symmetric = true;

    Eigen::MatrixXd temp_field, temp_field2;
    //igl::copyleft::comiso::nrosy(V, F, b, bc, VectorXi(), VectorXd(), MatrixXd(), 1, 0.5, temp_field, S);
  
    directional::complex_field(V, F, b,  bc , 4, complexField);
  
    // Convert it to raw field
    directional::complex_to_raw(V,F,complexField,4,raw, true);

    igl::streamlines_init(V, F, raw, false, sl_data, sl_state);

    // Viewer Settings
    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V, F);
    viewer.callback_pre_draw = &pre_draw;
    viewer.callback_key_down = &key_down;

    viewer.core.show_lines = false;

    viewer.core.is_animating = false;
    viewer.core.animation_max_fps = 30.;

    // Paint mesh grayish
    Eigen::MatrixXd C;
    C.setConstant(viewer.data.V.rows(), 3, 0.1);
    viewer.data.set_colors(C);


    // Draw vector field on sample points
    igl::StreamlineState sl_state0;
    sl_state0 = sl_state;
    igl::streamlines_next(V, F, sl_data, sl_state0);
    Eigen::MatrixXd v = sl_state0.end_point - sl_state0.start_point;
    v.rowwise().normalize();

    viewer.data.add_edges(sl_state0.start_point,
                          sl_state0.start_point + 0.059 * v,
                          Eigen::RowVector3d::Constant(1.0f));

    cout <<
    "Press [space] to toggle animation" << endl;
    viewer.launch();
}
