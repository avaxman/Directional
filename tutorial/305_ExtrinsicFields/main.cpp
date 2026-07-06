#include <iostream>
#include <fstream>
#include <directional/readOBJ.h>
#include <directional/readOFF.h>
#include <directional/CartesianField.h>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/read_raw_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/principal_matching.h>
#include <polyscope/polyscope.h>
#include <directional/ExtrinsicVertexTangentBundle.h>

using namespace std;

int N = 1;
double interpResolution = 20;
directional::TriMesh mesh;
directional::ExtrinsicVertexTangentBundle vtb;
directional::PCFaceTangentBundle ftb;
directional::CartesianField rawFaceField, powerFaceField;
directional::CartesianField rawVertexField, powerVertexField;;
directional::DirectionalViewer viewer;


typedef enum {FACE_FIELD, VERTEX_FIELD, INTERPOLATED_FIELD} ViewingModes;
ViewingModes viewingMode=FACE_FIELD;



// Compute interpolated normal at barycentric coordinates
Eigen::Vector3d compute_normal(const Eigen::Vector3d& baryCoords,
                               const Eigen::Vector3d& ni,
                               const Eigen::Vector3d& nj,
                               const Eigen::Vector3d& nk) {
    Eigen::Vector3d m = baryCoords(0)*ni + baryCoords(1)*nj + baryCoords(2)*nk;
    return m.normalized();
}

// Compute parallel-transported vector at point. The vectors are considered to be representatives
Eigen::Vector3d interpolate_vector(const Eigen::Vector3d& baryCoords,
                                      const Eigen::Vector3d& ni,
                                      const Eigen::Vector3d& nj,
                                      const Eigen::Vector3d& nk,
                                      const Eigen::Vector3d& vi,
                                      const Eigen::Vector3d& vj,
                                      const Eigen::Vector3d& vk,
                                      const int N) {
    Eigen::Vector3d np = compute_normal(baryCoords, ni,nj,nk);
    Eigen::Matrix3d Ri = directional::shortest_rotation_matrix(ni, np);
    Eigen::Matrix3d Rj = directional::shortest_rotation_matrix(nj, np);
    Eigen::Matrix3d Rk = directional::shortest_rotation_matrix(nk, np);
    
    Eigen::Vector3d Vxp = Ri*vi;
    Vxp.normalize();
    Eigen::Vector3d Vyp = np.cross(Vxp);
    
    Eigen::Vector3d Vxi = vi;
    Vxi.normalize();
    Eigen::Vector3d Vyi = ni.cross(Vxi);
    
    Eigen::Vector3d Vxj = vj;
    Vxj.normalize();
    Eigen::Vector3d Vyj = nj.cross(Vxj);
    
    Eigen::Vector3d Vxk = vk;
    Vxk.normalize();
    Eigen::Vector3d Vyk = nk.cross(Vxk);
    
    Eigen::Vector3d transVxi = Ri*Vxi;
    std::complex<double> connip = std::complex(transVxi.dot(Vxp), transVxi.dot(Vyp));
    Eigen::Vector3d transVxj = Rj*Vxj;
    std::complex<double> connjp = std::complex(transVxj.dot(Vxp), transVxj.dot(Vyp));
    Eigen::Vector3d transVxk = Rk*Vxk;
    std::complex<double> connkp = std::complex(transVxk.dot(Vxp), transVxk.dot(Vyp));
    
    std::complex<double> cvi(vi.dot(Vxi), vi.dot(Vyi));
    std::complex<double> cvj(vj.dot(Vxj), vj.dot(Vyj));
    std::complex<double> cvk(vk.dot(Vxk), vk.dot(Vyk));
    
    std::complex<double> pvi = pow(cvi, N);
    std::complex<double> pvj = pow(cvj, N);
    std::complex<double> pvk = pow(cvk, N);
    
    std::complex<double> pvp = baryCoords(0)*pvi*pow(connip, N) + baryCoords(1)*pvj*pow(connjp, N)+ baryCoords(2)*pvk*pow(connkp, N);
    
    std::complex<double> cvp = pow(pvp, 1.0/(double)N);
    Eigen::Vector3d vp = cvp.real()*Vxp + cvp.imag()*Vyp;
    
    return vp;
    
    //return baryCoords(0)*(Ri*vi) + baryCoords(1)*(Rj*vj)+ baryCoords(2)*(Rk*vk);
}


Eigen::MatrixXd barycentric_samples(int interpResolution)
{
    int n = interpResolution + 1;
    int numSamples = (n + 1) * (n + 2) / 2;
    
    Eigen::MatrixXd B(numSamples, 3);
    
    int row = 0;
    for (int i = 0; i <= n; ++i)
    {
        for (int j = 0; j <= n - i; ++j)
        {
            int k = n - i - j;
            
            B(row, 0) = double(i) / n;
            B(row, 1) = double(j) / n;
            B(row, 2) = double(k) / n;
            
            row++;
        }
    }
    
    return B;
}

//so as to not lose precision
Eigen::Matrix3d normalize_rotation_matrix(const Eigen::Matrix3d& R)
{
    // Compute SVD
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(R, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();
    
    // Compute nearest rotation
    Eigen::Matrix3d Rn = U * V.transpose();
    
    // Ensure right-handedness
    if (Rn.determinant() < 0)
    {
        U.col(2) *= -1;  // Flip last column
        Rn = U * V.transpose();
    }
    
    return Rn;
}


Eigen::Matrix3d inner_triangle_transport(const Eigen::Vector3d& coordSource,
                                         const Eigen::Vector3d& coordTarget,
                                         const Eigen::Vector3d& ni,
                                         const Eigen::Vector3d& nj,
                                         const Eigen::Vector3d& nk,
                                         int numSamples = 100){
    Eigen::Vector3d prevNormal = compute_normal(coordSource, ni, nj, nk);
    Eigen::Vector3d currNormal;
    Eigen::Matrix3d accumRot = Eigen::Matrix3d::Identity();
    for (int i=1;i<numSamples;i++){
        double t = (double)(i)/(numSamples-1);
        Eigen::Vector3d currCorrds = coordSource*(1-t) + coordTarget*t;
        currNormal = compute_normal(currCorrds, ni, nj, nk);
        Eigen::Matrix3d currR = directional::shortest_rotation_matrix(prevNormal,currNormal);
        prevNormal = currNormal;
        accumRot = currR*accumRot;
    }
    return accumRot;
}


void locate_singularities(const directional::CartesianField& vertexField, Eigen::MatrixXd& singLocations, Eigen::VectorXi& singIndices){
    using namespace Eigen;
    using namespace std;
    using namespace directional;
    singLocations.resize(vertexField.singIndices.size(), 3);
    singIndices = vertexField.singIndices;
    Eigen::VectorXi singFaces = vertexField.singLocalCycles;
    const directional::TriMesh& mesh = *((directional::ExtrinsicVertexTangentBundle*)vertexField.tb)->mesh;
    for (int i=0;i<singIndices.size();i++){
        std::cout<<"singFaces(i): "<<singFaces(i)<<std::endl;
        std::cout<<"singIndices(i): "<<singIndices(i)<<std::endl;
        //finding the index by line bisection and then on the line itself
        Vector3d ni = mesh.vertexNormals.row(mesh.F(singFaces(i),0)).transpose();
        Vector3d nj = mesh.vertexNormals.row(mesh.F(singFaces(i),1)).transpose();
        Vector3d nk = mesh.vertexNormals.row(mesh.F(singFaces(i),2)).transpose();
        Vector3d vi = vertexField.extField.row(mesh.F(singFaces(i),0)).head(3).transpose();
        Vector3d vj = vertexField.extField.row(mesh.F(singFaces(i),1)).head(3).transpose();
        Vector3d vk = vertexField.extField.row(mesh.F(singFaces(i),2)).head(3).transpose();
        vi.normalize(); vj.normalize(); vk.normalize();
        double baryBegin = 0.0;
        double baryEnd = 1.0;
        double fullHolonomy = vertexField.tb->cycleCurvatures(singFaces(i));
        //std::cout<<"fullHolonomy: "<<fullHolonomy<<std::endl;
        double currBaryLine;
        std::cout<<"Searching Bi"<<std::endl;
        while (baryEnd-baryBegin>1e-7){
            std::cout<<"baryBegin, baryEnd: "<<baryBegin<<", "<<baryEnd<<std::endl;
            currBaryLine = (baryBegin+baryEnd)/2.0;
            Eigen::Vector3d coordSource{currBaryLine, 1.0 - currBaryLine, 0.0};
            Eigen::Vector3d coordTarget{currBaryLine, 0.0, 1.0 - currBaryLine};
            Eigen::Vector3d vij = interpolate_vector(coordSource, ni, nj, nk, vi, vj, vk, N);
            Eigen::Vector3d vki = interpolate_vector(coordTarget, ni, nj, nk, vi, vj, vk, N);
            Eigen::Vector3d nij = compute_normal(coordSource, ni,nj,nk);
            Eigen::Vector3d nki = compute_normal(coordTarget, ni,nj,nk);
            Eigen::Matrix3d Rijki = shortest_rotation_matrix(nij, nki);
            Eigen::Matrix3d Riij = shortest_rotation_matrix(ni, nij);
            Eigen::Matrix3d Rkii = shortest_rotation_matrix(nki, ni);
            Eigen::Matrix3d loopR = Rkii*Rijki*Riij;
            Eigen::Vector3d w;
            w << loopR(2,1) - loopR(1,2),
                 loopR(0,2) - loopR(2,0),
                 loopR(1,0) - loopR(0,1);

            double sin_theta_signed = 0.5 * w.dot(ni);  // signed sin(theta)
            double cos_theta = 0.5 * (loopR.trace() - 1.0);
            double holonomy = std::atan2(sin_theta_signed, cos_theta);
            auto wrap = [](double x) {
                x = std::fmod(x + std::numbers::pi, 2.0 * std::numbers::pi);
                return x < 0 ? x + 2.0 * std::numbers::pi - std::numbers::pi : x - std::numbers::pi;
            };
            double effortiij = wrap((double)N*std::atan2(nij.dot((Riij * vi).cross(vij)), (Riij * vi).dot(vij)));
            double effortijki = wrap((double)N*std::atan2(nki.dot((Rijki * vij).cross(vki)), (Rijki * vij).dot(vki)));
            double effortkii = wrap((double)N*std::atan2(ni.dot((Rkii * vki).cross(vi)), (Rkii * vki).dot(vi)));
            double index = (effortiij+effortijki+effortkii + (double)N*holonomy)/(2*std::numbers::pi);
            if (std::abs(index)<10e-5)
                baryEnd = currBaryLine;
            else
                baryBegin = currBaryLine;
            
        }
        double Bi = currBaryLine;
        //Doing the same for Bj
        baryBegin = 0.0; baryEnd = 1.0;
        while (baryEnd-baryBegin>1e-7){
            currBaryLine = (baryBegin+baryEnd)/2.0;
            Eigen::Vector3d coordSource{0.0, currBaryLine, 1.0 - currBaryLine};
            Eigen::Vector3d coordTarget{ 1.0 - currBaryLine, currBaryLine, 0.0};
            Eigen::Vector3d vjk = interpolate_vector(coordSource, ni, nj, nk, vi, vj, vk, N);
            Eigen::Vector3d vij = interpolate_vector(coordTarget, ni, nj, nk, vi, vj, vk, N);
            Eigen::Vector3d njk = compute_normal(coordSource, ni,nj,nk);
            Eigen::Vector3d nij = compute_normal(coordTarget, ni,nj,nk);
            Eigen::Matrix3d Rjkij = shortest_rotation_matrix(njk, nij);
            Eigen::Matrix3d Rjjk = shortest_rotation_matrix(nj, njk);
            Eigen::Matrix3d Rijj = shortest_rotation_matrix(nij, nj);
            Eigen::Matrix3d loopR = Rijj*Rjkij*Rjjk;
          
            
            Eigen::Vector3d w;
            w << loopR(2,1) - loopR(1,2),
                 loopR(0,2) - loopR(2,0),
                 loopR(1,0) - loopR(0,1);

            double sin_theta_signed = 0.5 * w.dot(nj);  // signed sin(theta)
            double cos_theta = 0.5 * (loopR.trace() - 1.0);
            double holonomy = std::atan2(sin_theta_signed, cos_theta);
            auto wrap = [](double x) {
                x = std::fmod(x + std::numbers::pi, 2.0 * std::numbers::pi);
                return x < 0 ? x + 2.0 * std::numbers::pi - std::numbers::pi : x - std::numbers::pi;
            };
            double effortjjk = wrap((double)N*std::atan2(njk.dot((Rjjk * vj).cross(vjk)), (Rjjk * vj).dot(vjk)));
            double effortjkij = wrap((double)N*std::atan2(nij.dot((Rjkij * vjk).cross(vij)), (Rjkij * vjk).dot(vij)));
            double effortijj = wrap((double)N*std::atan2(nj.dot((Rijj * vij).cross(vj)), (Rijj * vij).dot(vj)));
            double index = (effortjjk+effortjkij+effortijj + N*holonomy)/(2*std::numbers::pi);
            if (std::abs(index)<10e-5)
                baryEnd = currBaryLine;
            else
                baryBegin = currBaryLine;
        
        }
        double Bj = currBaryLine;
        double Bk = 1.0-Bi-Bj;
        singLocations.row(i) = Bi*mesh.V.row(mesh.F(singFaces(i),0)) + Bj*mesh.V.row(mesh.F(singFaces(i),1)) + Bk*mesh.V.row(mesh.F(singFaces(i),2));
    }
    
}

Eigen::Vector3d rotateAroundAxis(
    const Eigen::Vector3d& v,
    const Eigen::Vector3d& n,
    double t)
{
    const double norm = n.norm();
    if (norm == 0.0) return v; // or assert

    Eigen::Vector3d axis = n / norm;

    // Use quaternion (very stable for small angles)
    Eigen::Quaterniond q(Eigen::AngleAxisd(t, axis));
    return q * v;
}


void extrinsic_interpolation(const directional::CartesianField& vertexField, const int interpResolution, const bool normalize, Eigen::MatrixXd& newSources, Eigen::MatrixXd& newRawField){
    using namespace Eigen;
    using namespace std;
    using namespace directional;
    MatrixXd barySamples = barycentric_samples(interpResolution);
    newSources.resize(barySamples.rows()*mesh.F.rows(), 3);
    newRawField.resize(barySamples.rows()*mesh.F.rows(), 3*N);
    const directional::TriMesh& mesh = *((directional::ExtrinsicVertexTangentBundle*)vertexField.tb)->mesh;
    for (int i=0;i<mesh.F.rows();i++){
        Vector3d ni = mesh.vertexNormals.row(mesh.F(i,0)).transpose();
        Vector3d nj = mesh.vertexNormals.row(mesh.F(i,1)).transpose();
        Vector3d nk = mesh.vertexNormals.row(mesh.F(i,2)).transpose();
        Vector3d repvi = vertexField.extField.row(mesh.F(i,0)).head(3).transpose();
        Vector3d repvj = vertexField.extField.row(mesh.F(i,1)).head(3).transpose();
        Vector3d repvk = vertexField.extField.row(mesh.F(i,2)).head(3).transpose();
        
        //parallel transport and then merge
        for (int j=0;j<barySamples.rows();j++){
            Vector3d np = barySamples(j,0)*ni+barySamples(j,1)*nj+barySamples(j,2)*nk;
            RowVector3d p = barySamples(j,0)*mesh.V.row(mesh.F(i,0))+barySamples(j,1)*mesh.V.row(mesh.F(i,1))+barySamples(j,2)*mesh.V.row(mesh.F(i,2));
            np.normalize();
            Eigen::Vector3d vp = interpolate_vector(barySamples.row(j),ni, nj, nk, repvi, repvj, repvk,N);
            if (normalize) vp.normalize();
            newRawField.row(barySamples.rows()*i+j).head(3) = vp.transpose();
            //reproducing the rest of the field by rotation around the normal
            for (int n=1;n<N;n++)
                newRawField.row(barySamples.rows()*i+j).segment(3*n,3) = rotateAroundAxis(vp, np, (double)n*2.0*(std::numbers::pi)/(double)N).transpose();

            newSources.row(barySamples.rows()*i+j) = p;
        }
    }
}


void update_visualization()
{
    viewer.toggle_singularities(viewingMode==FACE_FIELD, 0);
    viewer.toggle_singularities(viewingMode==VERTEX_FIELD, 1);
    viewer.toggle_cartesian_field(viewingMode==FACE_FIELD,0);
    viewer.toggle_cartesian_field(viewingMode==VERTEX_FIELD,1);
    viewer.toggle_cartesian_field(viewingMode==INTERPOLATED_FIELD,2);
}

void callbackFunc()
{
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
    
    const char* items[] = { "Face-based field", "Extrinsic Vertex-based field", "Interpolated Field"};
    static const char* current_item = NULL;
    
    float combo_width = 0.0f;
    if (combo_width == 0.0f) {
        ImGuiStyle& style = ImGui::GetStyle();
        for (auto& item : items)
            combo_width = std::max(combo_width, ImGui::CalcTextSize(item).x);
        combo_width += style.FramePadding.x * 5.0 + ImGui::GetFontSize() + style.ItemInnerSpacing.x;
    }
    
    ImGui::PushItemWidth(combo_width);
    if (ImGui::BeginCombo("Viewing mode", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)){
                current_item = items[n];
                switch (n){
                    case 0:
                        viewingMode = FACE_FIELD;
                        break;
                    case 1:
                        viewingMode = VERTEX_FIELD;
                        break;
                    case 2:
                        viewingMode = INTERPOLATED_FIELD;
                        break;
                }
                update_visualization();
            }
            if (is_selected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }
    
    //if (ImGui::Button("Save Raw Field"))
    //    directional::write_raw_field(TUTORIAL_OUTPUT_PATH "/polyvector.rawfield", (viewingMode==HARD_PRESCRIPTION ? rawFieldHard : rawFieldSoft));
    
    ImGui::PopItemWidth();
}


int main(int argc, char *argv[])
{
    directional::readOBJ(TUTORIAL_DATA_PATH "/bunny1k.obj", mesh);
    viewer.init();
    viewer.set_callback(&callbackFunc);
    
    ftb.init(mesh);
    vtb.init(mesh);
    
    Eigen::VectorXi constFaces, constVertices;
    Eigen::MatrixXd constVectors;
    constFaces.resize(1);
    constFaces<<0;
    constVectors.resize(1,3);
    constVectors<<mesh.V.row(mesh.F(0,2))-mesh.V.row(mesh.F(0,1));
    constVertices.resize(1);
    constVertices<<mesh.F(0,1);
    
    bool normalizeField = true;
    directional::power_field(vtb, constVertices, constVectors, Eigen::VectorXd::Constant(constVertices.size(),-1.0), N, powerVertexField, normalizeField);
    directional::power_field(ftb, constFaces, constVectors, Eigen::VectorXd::Constant(constFaces.size(),-1.0), N, powerFaceField, normalizeField);
    
    //computing power fields
    directional::power_to_raw(powerFaceField, N, rawFaceField,normalizeField);
    directional::power_to_raw(powerVertexField, N, rawVertexField,normalizeField);
    
    directional::principal_matching(rawFaceField);
    directional::principal_matching(rawVertexField);
    
    //Interpolated field
    bool normalizeInterpolation = true;
    Eigen::MatrixXd newSources, newRawField;
    extrinsic_interpolation(rawVertexField, interpResolution, normalizeInterpolation, newSources, newRawField);
    Eigen::MatrixXd singLocations;
    Eigen::VectorXi singIndices;
    locate_singularities(rawVertexField, singLocations, singIndices);
    viewer.set_surface_mesh(mesh);
    viewer.set_cartesian_field(rawFaceField, "Face-Based Field",  0);
    viewer.set_cartesian_field(rawVertexField, "Vertex-Based Field", 1);
    viewer.set_raw_field(newSources, newRawField, "Interpolated Field", 2, 0.3*mesh.avgEdgeLength/(double)(interpResolution));
    viewer.update_singularities_locations(rawVertexField, singLocations, 1, "Exact Singularities");
    update_visualization();
    viewer.launch();
    
    return 0;
}
