#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/combing.h>
#include <directional/directional_viewer.h>
#include <directional/readOBJ.h>

int N;
directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField origField, combedField;
directional::DirectionalViewer viewer;

bool showCombed=false;


void callbackFunc() {
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
    if (showCombed) {
        if (ImGui::Button("Show Original Field")) {
            showCombed=false;
            viewer.set_cartesian_field(origField);
        }
    } else {
        if (ImGui::Button("Show Combed Field")) {
            showCombed=true;
            viewer.set_cartesian_field(combedField);
        }
        
    }
    viewer.toggle_combed_colors(true, false);
    ImGui::PopItemWidth();
}



int main()
{
    
    directional::readOBJ(TUTORIAL_DATA_PATH "/lilium.obj", mesh);
    ftb.init(mesh);
    directional::read_raw_field(TUTORIAL_DATA_PATH "/lilium.rawfield", ftb, N, origField);
    
    //computing the principal matching and then combing the field to trivialize the matching anywhere but a set of seams
    directional::principal_matching(origField);
    directional::combing(origField, combedField);
    Eigen::VectorXd edgeData=origField.effort;
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.set_cartesian_field(origField);
    showCombed = false;
    std::vector<int> HighlightedEdgesList;
    std::for_each(combedField.matching.data(), combedField.matching.data() + combedField.matching.size(), [&, i = 0](int val) mutable {
        if (val != 0) HighlightedEdgesList.push_back(i);
        ++i;
    });
    Eigen::VectorXi HighlightedEdges = Eigen::Map<Eigen::VectorXi>(HighlightedEdgesList.data(), HighlightedEdgesList.size());
    viewer.highlight_edges(HighlightedEdges, "Combing seams");
    viewer.set_surface_edge_data(edgeData, "Matching effort");
    viewer.set_callback(callbackFunc);
    viewer.launch();
}


