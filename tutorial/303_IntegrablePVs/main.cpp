#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/local_basis.h>
#include <igl/avg_edge_length.h>
#include <igl/is_border_vertex.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
//#include <directional/polyvector_field_matchings.h>
//#include <directional/polyvector_field_singularities_from_matchings.h>
//#include <directional/polyvector_field_cut_mesh_with_singularities.h>
//#include <directional/polyvector_field_comb_from_matchings_and_cuts.h>
#//include <directional/polyvector_field_poisson_reconstruction.h>
//#include <directional/integrable_polyvector_fields.h>
#include <directional/glyph_lines_raw.h>
#include <directional/line_cylinders.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/curl_combing.h>
#include <directional/integrable_polyvector_fields.h>
#include <igl/cut_mesh.h>
#include <igl/slice.h>
#include <igl/false_barycentric_subdivision.h>

#include <iostream>
#include <fstream>
#include <igl/matlab_format.h>

#include "tutorial_shared_path.h"

using namespace std;

// Input mesh
//std::vector<bool> V_border;
//std::vector<std::vector<int> > VF, VFi;
//std::vector<std::vector<int> > VV;
//Eigen::MatrixXi TT, TTi;
//Eigen::MatrixXi E, E2F, F2E;

Eigen::VectorXi matchingOrig, matchingCF, combedMatchingOrig, combedMatchingCF;
Eigen::VectorXd effortOrig, effortCF, combedEffortOrig, combedEffortCF;
Eigen::MatrixXi FMesh, FField, FSings, FSeams;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings, VSeams, barycenters;
Eigen::MatrixXd CMesh, CField, CSings, CSeams;
Eigen::MatrixXd rawFieldOrig, rawFieldCF;
Eigen::MatrixXd combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge
igl::opengl::glfw::Viewer viewer;
Eigen::MatrixXd glyphPrincipalColors(5,3);

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;

// "Subdivided" mesh obtained by splitting each triangle into 3 (only needed for display)
Eigen::MatrixXd Vbs;
Eigen::MatrixXi Fbs;


//for averagin curl to faces, for visualization
Eigen::SparseMatrix<double> AE2F;

// Scale for visualizing textures
double uv_scale;
double curlMax;
int N;


//Eigen::MatrixXd scalars_ori;

// Data for curl-free PolyVector field
//Eigen::MatrixXd scalars;


// percentage of constrained faces
double constraint_percentage = 0.002;

// Random length factor
double rand_factor = 5;

// The set of parameters for calculating the curl-free fields
igl::integrable_polyvector_fields_parameters params;

// Solver data (needed for precomputation)
igl::IntegrableFieldSolverData<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::MatrixXd, Eigen::MatrixXd> ipfdata;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;


typedef enum {ORIGINAL_FIELD, ORIGINAL_CURL, OPTIMIZED_FIELD, OPTIMIZED_CURL} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;

int display_mode = 1;

int iter = 0;

// Create a texture that hides the integer translation in the parametrization
void line_texture(Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_R,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_G,
                  Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> &texture_B)
{
  unsigned size = 128;
  unsigned size2 = size/2;
  unsigned lineWidth = 3;
  texture_R.setConstant(size, size, 255);
  for (unsigned i=0; i<size; ++i)
    for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
      texture_R(i,j) = 0;
  for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
    for (unsigned j=0; j<size; ++j)
      texture_R(i,j) = 0;

  texture_G = texture_R;
  texture_B = texture_R;
}

// Create a random set of tangent vectors
/*void generate_constraints()
{
  b.resize(42,1);
  b<<
  663,
  513,
  3872,
  2601,
  3549,
  2721,
  3796,
  594,
  868,
  1730,
  1581,
  3081,
  1471,
  1650,
  454,
  2740,
  2945,
  3808,
  3679,
  3589,
  450,
  2656,
  1791,
  1792,
  2917,
  3744,
  1536,
  2809,
  3866,
  3658,
  1665,
  2670,
  1601,
  1793,
  3614,
  524,
  2877,
  449,
  455,
  3867,
  3871,
  2592;

blevel.setOnes(b.rows(),1);
bc.resize(b.rows(),6);
bc<<
-0.88046298335147721,0.27309862654264377,0.38755912468723769,-0.350632259447135,-0.92528970817792766,-0.14455440005410564,
0.91471470003012889,0.392936119054695,-0.094330397492144599,0.32234487030777614,-0.85027369799342767,-0.41608703787410195,
0.94335566040683105,0.073867667925654024,-0.32345581709658111,0.19950360079371404,-0.90525435056476755,0.37511714710727789,
-0.92054671613540229,0.15077598183983737,0.36036141124232496,-0.27998315313687211,-0.89796618385425386,-0.33950871360506074,
0.88944399663239937,0.23035525634795684,-0.39474780902172396,0.27297422303141039,-0.96047177712172194,0.054580572670497061,
-0.83112706922096102,-0.55599928943162547,0.0096221078617792517,0.52546831822766438,-0.79091522174894457,-0.31358596675362826,
0.90724658517664569,-0.41046292080872998,-0.091781394228251156,-0.34252813327252363,-0.84767620917196618,0.40511667741613094,
-0.8932101465465786,0.23975524191487588,0.38038540729184012,-0.33645713296414853,-0.91759364410558497,-0.21170380718016926,
-0.87839308390284521,0.27039404931731387,0.39409725734320344,-0.29712518405497651,-0.95481177255558192,-0.0071487054467167244,
-0.91448048788760505,-0.17055891298176448,0.36692655188106316,0.29811257890714044,-0.89715315396744022,0.32595261714489804,
0.82798126471567035,-0.56074230404745851,0.003885065171440813,-0.53510484459763941,-0.78801608401899037,0.30445600111594384,
-0.87831929581593793,0.25312706437601257,0.40556368658667746,-0.26531767440854004,-0.9637845762158106,0.026941089342378936,
-0.87482003689209031,-0.27011021313654948,0.4021571531272935,0.32303198334357713,-0.94388894366288889,0.0687313594225408,
0.87408456883093666,-0.48487387939766946,-0.029554823793924323,-0.43846604347950752,-0.81368808449189478,0.38165328489525413,
0.94988212941972827,-0.041936960956176939,-0.30978255521381903,-0.16637246118261648,-0.90677959514398765,-0.3873899456497869,
0.87516493768857484,-0.27181042881473483,-0.40025669591913515,-0.36755520380602424,-0.91147911093961553,-0.18468622708756641,
-0.87064244687577641,0.27922257819020818,0.40498948323008854,-0.32176729617260508,-0.94599403842079244,-0.039510585747255161,
-0.91274615133859638,-0.1832657459904497,0.36511385835536858,0.29782083933521708,-0.91026141603074595,0.28762284704690655,
0.875611546674125,0.28258715176515403,-0.39172556846369444,0.36000975242683186,-0.92250014843287764,0.13923524804764778,
0.76763693171870195,-0.64088483679642994,0.00040868803559811201,-0.63058113310112196,-0.75518119878562417,0.17907761327874747,
0.963064265211517,0.17044712473620016,-0.20845862597111031,0.061832174999749308,-0.89345471128813481,-0.44487690546019126,
-0.88228998376691692,-0.46837234310148523,0.046815945597605227,0.41604986062280985,-0.82249303168905052,-0.38782434980116298,
-0.96608602970701829,0.11121907649833783,0.23304098400879364,0.010641270548624404,-0.88457418950525291,0.46627810008860171,
-0.96329451047686887,0.055809409647239343,0.26258140810033831,0.07182051046944142,-0.88891411988025926,0.45240855623364673,
-0.71244584326772997,-0.70122065397026967,-0.026655484539588895,0.70046172163981768,-0.70836773631021255,-0.086997279682342638,
0.88646445996853696,0.2549240118236365,-0.38625705094979518,0.35132981358631576,-0.91395520354543514,0.20310895597591658,
-0.86109327343809683,-0.30822574449366841,0.40437020769461601,0.37896596246993836,-0.91928725525816557,0.10628142645421024,
0.86443027504389158,-0.29669958642983363,-0.40586901212079213,-0.37200560813855077,-0.92052106924988175,-0.11938504337027039,
0.95370728000967508,-0.24689991217686594,-0.17170572915195079,-0.14736898596800915,-0.88138264597997584,0.4488284898935197,
-0.81439393313167019,0.57995723960933832,0.020300785774083896,-0.55494919604589421,-0.78855235001798585,0.26498411478639117,
0.89527216270596455,0.22395367264061938,-0.38513959442592183,0.33680943342191538,-0.90609008785063272,0.25604717974787594,
-0.9003647006267198,0.20802946062196581,0.38218732236782926,-0.32431023000528064,-0.90640636884236436,-0.27064805418831556,
-0.87050937437709508,-0.28614105672408718,0.40042068475344922,0.37746788793940733,-0.91025870352880611,0.17013843253251126,
-0.95715079751532439,0.0030851788865496879,0.28957353554324744,0.12908381923211401,-0.89056292562302997,0.43615942397041058,
-0.87324677619075319,-0.28591869514051466,0.39457644080913162,0.3438918663433696,-0.93530416305293196,0.083333707698197687,
0.91999856277124803,-0.1621255206103257,-0.35681642348085474,-0.27672206872177485,-0.91342693749618353,-0.2984562389005877,
-0.8669467282645521,0.29036174243712859,0.40508447128995645,-0.34873789125620602,-0.93406639205959408,-0.07682355385522964,
-0.9570365266718136,-0.22821899053183647,0.17887755302603078,0.12590409644663494,-0.88275887883510706,-0.45264215483728532,
-0.94033215083998489,0.087395510869996196,0.32884262311388451,-0.2131320783418921,-0.90465024471116184,-0.36902933748646671,
-0.96131014054749453,0.18866284908038999,0.20072155603578434,-0.08260532909072589,-0.89255302833360861,0.44331191188407898,
-0.95240414686152941,-0.02752900142620229,0.30359264668538755,0.15128346580527452,-0.9073021943457209,0.39232134929083828,
-0.94070423353276911,-0.31552769387286655,0.12457053990729766,0.22959741970407915,-0.86253407908715607,-0.45091017650802745;

}*/

/*void drawCuts(igl::opengl::glfw::Viewer& viewer,
              const Eigen::MatrixXi &cuts)
{
  int maxCutNum = cuts.sum();
  Eigen::MatrixXd start(maxCutNum,3);
  Eigen::MatrixXd end(maxCutNum,3);
  int ind = 0;
  for (unsigned int i=0;i<F.rows();i++)
    for (int j=0;j<3;j++)
      if (cuts(i,j))
      {
        start.row(ind) = V.row(F(i,j));
        end.row(ind) = V.row(F(i,(j+1)%3));
        ind++;
      }
  viewer.data().add_edges(start, end , Eigen::RowVector3d(1.,0,1.));
}*/

/*void drawField(igl::opengl::glfw::Viewer &viewer,
               const Eigen::MatrixXd &field,
               const Eigen::RowVector3d &color)
{
  for (int n=0; n<2; ++n)
  {
    Eigen::MatrixXd VF = field.block(0,n*3,F.rows(),3);
    Eigen::VectorXd c = VF.rowwise().norm();
    viewer.data().add_edges(B - global_scale*VF, B + global_scale*VF , color);
  }
}*/

/*void drawConstraints(igl::opengl::glfw::Viewer &viewer)
{
  for (int n=0; n<2; ++n)
  {
    Eigen::MatrixXd Bc = igl::slice(B, b, 1);
    Eigen::MatrixXd color;
    color.setZero(b.rows(),3);
    color.col(2).setOnes();
    for (int i =0; i<b.rows(); ++i)
      if (blevel[i] ==1 && n>0)
        color.row(i)<<0.7,0.7,0.7;
    // Eigen::RowVector3d color; color<<0.5,0.5,0.5;
    viewer.data().add_edges(Bc - global_scale*bc.block(0,n*3,bc.rows(),3), Bc + global_scale*bc.block(0,n*3,bc.rows(),3) , color);
  }

}*/


/*void colorEdgeMeshFaces(const Eigen::VectorXd &values,
                        const double &minimum,
                        const double &maximum,
                        Eigen::MatrixXd &C)
{
  C.setConstant(Fbs.rows(),3,1);

  Eigen::MatrixXd colors;
  igl::jet(values, minimum, maximum, colors);

  for (int ei = 0; ei<E.rows(); ++ei)
  {
    const Eigen::RowVector3d &this_color = colors.row(ei);
    int f0 = E2F(ei,0);
    int f1 = E2F(ei,1);
    if(f0 != -1)
    {
      int i0 = -1;
      for (int k = 0; k<3; k++)
        if (F2E(f0,k)== ei)
        {
          i0 = k;
          break;
        }
      C.row(3*f0+i0) = this_color;
    }
    if(f1 != -1)
    {
      int i1 = -1;
      for (int k = 0; k<3; k++)
        if (F2E(f1,k)== ei)
        {
          i1 = k;
          break;
        }
      C.row(3*f1+i1) = this_color;
    }
  }

}*/

void update_triangle_mesh()
{

  if ((viewingMode ==ORIGINAL_FIELD)||(viewingMode ==OPTIMIZED_FIELD))
    viewer.data_list[0].set_colors(Eigen::RowVector3d::Constant(3,1.0));
  else{  //curl viewing - currently averaged to the face
    igl::jet(AE2F*(viewingMode==ORIGINAL_CURL ? rawFieldOrig: rawFieldCF), 0,curlMax, CMesh);
    viewer.data_list[0].set_colors(CMesh);
  }
}


void update_raw_field_mesh()
{
  using namespace std;
  using namespace Eigen;
  
  if ((viewingMode==ORIGINAL_CURL) || (viewingMode==OPTIMIZED_CURL)){
    for (int i=1;i<=3;i++){  //hide all other meshes
      viewer.data_list[i].show_faces=false;
      viewer.data_list[i].show_lines = false;
    }
  } else {
    Eigen::MatrixXd fullGlyphColor(FMesh.rows(),3*N);
    for (int i=0;i<FMesh.rows();i++)
      for (int j=0;j<N;j++)
        fullGlyphColor.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
    
    directional::glyph_lines_raw(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF), fullGlyphColor,VField, FField, CField);
    
    viewer.data_list[1].clear();
    viewer.data_list[1].set_mesh(VField, FField);
    viewer.data_list[1].set_colors(CField);
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    //singularity mesh
    directional::singularity_spheres(VMesh, FMesh, (viewingMode==ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode==ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), directional::defaultSingularityColors(N), VSings, FSings, CSings);
    
    viewer.data_list[2].clear();
    viewer.data_list[2].set_mesh(VSings, FSings);
    viewer.data_list[2].set_colors(CSings);
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    //seam mesh
    double avgEdgeLength = igl::avg_edge_length(VMesh, FMesh);
    std::vector<int> seamEdges;
    for (int i=0;i<EV.rows();i++)
      if ((viewingMode==ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF)(i)!=0)
        seamEdges.push_back(i);
    
    Eigen::MatrixXd P1(seamEdges.size(),3), P2(seamEdges.size(),3);
    for (int i=0;i<seamEdges.size();i++){
      P1.row(i)=VMesh.row(EV(seamEdges[i],0));
      P2.row(i)=VMesh.row(EV(seamEdges[i],1));
    }
    
    directional::line_cylinders(P1, P2, avgEdgeLength/25.0, Eigen::MatrixXd::Constant(FMesh.rows(), 3, 0.0), 6, VSeams, FSeams, CSeams);
    
    viewer.data_list[3].clear();
    viewer.data_list[3].set_mesh(VSeams, FSeams);
    viewer.data_list[3].set_colors(CSeams);
    viewer.data_list[3].show_faces = true;
    viewer.data_list[3].show_lines = false;
  }

}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{

  if ((key >= '1') && (key <='4'))
    display_mode = key - '1';
  
  if (key == 'A')
  {
    //do a batch of iterations
    printf("--Improving Curl--\n");
    for (int bi = 0; bi<5; ++bi)
    {
      Eigen::MatrixXd two_pv=rawFieldCF.block(0,0,rawFieldCF.rows(),6);
      printf("\n\n **** Batch %d ****\n", iter);
      igl::integrable_polyvector_fields_solve(ipfdata, params, two_pv, iter ==0);
      iter++;
      params.wSmooth *= params.redFactor_wsmooth;
      rawFieldCF.block(0,0,rawFieldCF.rows(),6)=two_pv;
      rawFieldCF.block(0,6,rawFieldCF.rows(),6)=-two_pv;
    }
    
    /*// Post process current field
    // Compute curl_minimizing matchings and curl
    printf("--Matchings and curl--\n");
    Eigen::MatrixXi match_ab, match_ba;  // matchings across interior edges
    double avgCurl = igl::polyvector_field_matchings(two_pv, V, F, true, true, match_ab, match_ba, curl);
    double maxCurl = curl.maxCoeff();
    printf("curl -- max: %.5g, avg: %.5g\n", maxCurl,  avgCurl);
    // Compute singularities
    printf("--Singularities--\n");
    igl::polyvector_field_singularities_from_matchings(V, F, match_ab, match_ba, singularities);
    printf("#singularities: %ld\n", singularities.rows());
    // Get mesh cuts based on singularities
    printf("--Cuts--\n");
    igl::polyvector_field_cut_mesh_with_singularities(V, F, singularities, cuts);
    // Comb field
    printf("--Combing--\n");
    Eigen::MatrixXd combed;
    igl::polyvector_field_comb_from_matchings_and_cuts(V, F, two_pv, match_ab, match_ba, cuts, combed);
    // Reconstruct integrable vector fields from combed field*/
    
    Eigen::VectorXi prinIndices;
    directional::curl_matching(VMesh, FMesh,EV, EF, FE, rawFieldCF, matchingCF, effortCF, curlCF);
    directional::effort_to_indices(VMesh,FMesh,EV, EF, effortCF,matchingCF, N,prinIndices);
    
    std::vector<int> singVerticesList;
    std::vector<int> singIndicesList;
    for (int i=0;i<VMesh.rows();i++)
      if (prinIndices(i)!=0){
        singVerticesList.push_back(i);
        singIndicesList.push_back(prinIndices(i));
      }
    
    singVerticesCF.resize(singVerticesList.size());
    singIndicesCF.resize(singIndicesList.size());
    for (int i=0;i<singVerticesList.size();i++){
      singVerticesCF(i)=singVerticesList[i];
      singIndicesCF(i)=singIndicesList[i];
    }
    
    directional::curl_combing(VMesh,FMesh, EV, EF, FE,rawFieldCF, combedFieldCF, combedMatchingCF, combedEffortCF);
  }

  update_raw_field_mesh();
  return false;
}

int main(int argc, char *argv[])
{

   cerr<<"Press keys 1-0 for various visualizations, 'A' to optimize for less curl." <<endl;
  
  // Load a mesh
  igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka.off", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka.rawfield", N, rawFieldOrig);

  /*V_border = igl::is_border_vertex(VMesh,FMesh);
  igl::adjacency_list(F, VV);
  igl::vertex_triangle_adjacency(VMesh,FMesh,VF,VFi);
  igl::triangle_triangle_adjacency(FMesh,TT,TTi);*/

  igl::edge_topology(VMesh, FMesh, EV, FE, EF);

  igl::false_barycentric_subdivision(VMesh, FMesh, Vbs, Fbs);

  //global_scale =  .2*igl::avg_edge_length(VMesh, FMesh);
  //uv_scale = 0.6/igl::avg_edge_length(V, F);
  igl::barycenter(VMesh, FMesh, barycenters);

  //igl::local_basis(VMesh,FMesh,B1,B2,B3);

  Eigen::VectorXi prinIndices;
  directional::curl_matching(VMesh, FMesh,EV, EF, FE, rawFieldOrig, matchingOrig, effortOrig, curlOrig);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effortOrig,matchingOrig, N,prinIndices);
  
  curlMax= curlOrig.maxCoeff();
  
  std::vector<int> singVerticesList;
  std::vector<int> singIndicesList;
  for (int i=0;i<VMesh.rows();i++)
    if (prinIndices(i)!=0){
      singVerticesList.push_back(i);
      singIndicesList.push_back(prinIndices(i));
    }
  
  singVerticesOrig.resize(singVerticesList.size());
  singIndicesOrig.resize(singIndicesList.size());
  for (int i=0;i<singVerticesList.size();i++){
    singVerticesOrig(i)=singVerticesList[i];
    singIndicesOrig(i)=singIndicesList[i];
  }
  
  directional::curl_combing(VMesh,FMesh, EV, EF, FE,rawFieldOrig, combedFieldOrig, combedMatchingOrig, combedEffortOrig);

  //trivial constraints
  Eigen::VectorXi b; b.resize(1); b<<0;
  Eigen::MatrixXd bc; bc.resize(1,6); bc<<rawFieldOrig.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b<<1;
  Eigen::MatrixXd twoVectorMat=rawFieldOrig.block(0,0,rawFieldOrig.rows(),6);
  igl::integrable_polyvector_fields_precompute(VMesh, FMesh, b, bc, blevel, twoVectorMat , ipfdata);
  
  rawFieldCF = rawFieldOrig;
  matchingCF = matchingOrig;
  effortCF = effortOrig;
  combedFieldCF =combedFieldOrig;
  combedMatchingCF =combedMatchingOrig;
  combedEffortCF =combedEffortOrig;
  curlCF = curlOrig;
  singVerticesCF = singVerticesOrig;
  singIndicesCF = singIndicesOrig;

  // Replace the standard texture with an integer shift invariant texture
  line_texture(texture_R, texture_G, texture_B);

  //triangle mesh setup
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(Eigen::RowVector3d::Constant(3,1.0));
  
  //apending and updating raw field mesh
  viewer.append_mesh();
 
  //singularity mesh
  viewer.append_mesh();
 
  //seam mesh
  viewer.append_mesh();
  

  update_triangle_mesh();
  update_raw_field_mesh();
  viewer.selected_data_index=0;
  
  glyphPrincipalColors<<1.0,0.0,0.5,
  0.0,1.0,0.5,
  1.0,0.5,0.0,
  0.0,0.5,1.0,
  0.5,1.0,0.0;
  
  viewer.callback_key_down = &key_down;
  viewer.data().show_lines = false;
  viewer.launch();

  return 0;
}
