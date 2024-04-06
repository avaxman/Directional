// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_VIEWER_H
#define DIRECTIONAL_VIEWER_H

#include <Eigen/Core>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <directional/streamlines.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>

/***
 This class implements the Directional viewer, as an extension of the libigl viewer (as a wrapper). This
 viewer providers specialized functionality for outputting directional fields and their combinatorial and geometic properties.
 The numerous tutorial examples highlight its functionality.
 ***/


namespace directional
{

    class DirectionalViewer{
    private:
        std::vector<const TriMesh*> meshList;  //meshes that are being viewed
        std::vector<const CartesianField*> fieldList;

        std::vector<polyscope::SurfaceMesh*> psSurfaceMeshList;
        std::vector<polyscope::PointCloud*> psFieldSourceList;
        std::vector<std::vector<polyscope::PointCloudVectorQuantity*>> psFieldList;
        std::vector<polyscope::PointCloud*> psSingList;
        std::vector<polyscope::CurveNetwork*> psStreamlineList;
        //std::vector<Eigen::MatrixXd> fieldColors;
        std::vector<directional::StreamlineData> slData;
        std::vector<directional::StreamlineState> slState;

        //std::vector<Eigen::MatrixXd> edgeVList;  //edge-diamond vertices list
        //std::vector<Eigen::MatrixXi> edgeFList;  //edge-diamond faces list
        //std::vector<Eigen::VectorXi> edgeFEList;  //edge-diamond faces->original mesh edges list

    public:
        DirectionalViewer(){}
        ~DirectionalViewer(){}

        void init(){
            polyscope::init();
        }

        void launch(){
            polyscope::show();
        }

        void set_callback(void (*callbackFunc)()){
            polyscope::state::userCallback = callbackFunc;
        }

        void inline set_mesh(const TriMesh& mesh,
                             const int meshNum=0,
                             const std::string meshName="Mesh")
        {
            if (meshList.size()<meshNum+1) {
                meshList.resize(meshNum + 1);
                psSurfaceMeshList.resize(meshNum + 1);
            }

            meshList[meshNum]=&mesh;
            psSurfaceMeshList[meshNum]=polyscope::registerSurfaceMesh(meshName + std::to_string(meshNum), mesh.V, mesh.F);
        }

        /*void inline set_mesh_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                                        const int meshNum=0)
        {
            Eigen::MatrixXd meshColors;
            if (C.rows()==0)
                meshColors=default_mesh_color();
            else
                meshColors=C;

            if ((meshList[meshNum]->V.rows()==meshList[meshNum]->F.rows())&&(meshColors.rows()!=meshList[meshNum]->F.rows()))  //assume face_based was meant
                data_list[NUMBER_OF_SUBMESHES*meshNum].set_face_based(true);
            data_list[NUMBER_OF_SUBMESHES*meshNum].set_colors(meshColors);
            data_list[NUMBER_OF_SUBMESHES*meshNum].show_faces=true;
            if (edgeVList[meshNum].size()!=0){
                data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_faces=false;
                data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_lines=false;
            }
            selected_data_index=NUMBER_OF_SUBMESHES*meshNum;
            //CList[meshNum]=C;
        }*/

        void inline set_vertex_data(const Eigen::VectorXd& vertexData,
                                    const double minRange,
                                    const double maxRange,
                                    const std::string dataName="",
                                    const int meshNum=0)
        {
            psSurfaceMeshList[meshNum]->addVertexScalarQuantity("vertex data" + std::to_string(meshNum), vertexData);
        }

        void inline set_face_data(const Eigen::VectorXd& faceData,
                                  const double minRange,
                                  const double maxRange,
                                  const std::string dataName="",
                                  const int meshNum=0)
        {
            psSurfaceMeshList[meshNum]->addFaceScalarQuantity("face data" + std::to_string(meshNum), faceData);
        }

        void inline set_edge_data(const Eigen::VectorXd& edgeData,
                                  const double minRange,
                                  const double maxRange,
                                  const std::string dataName="",
                                  const int meshNum=0)
        {
            std::vector<int> permArr(edgeData.size());
            for (int i=0;i<permArr.size();i++)
                permArr[i]=i;
            psSurfaceMeshList[meshNum]->setEdgePermutation(permArr);
            psSurfaceMeshList[meshNum]->addEdgeScalarQuantity("edge data" + std::to_string(meshNum), edgeData);
        }



        //this function assumes face-based fields

        /*void inline set_selected_faces(const Eigen::VectorXi& selectedFaces, const int meshNum=0){
            assert(fieldList[meshNum]->tb->discTangType()==discTangTypeEnum::FACE_SPACES);
            Eigen::MatrixXd CMesh=directional::DirectionalViewer::default_mesh_color().replicate(meshList[meshNum]->F.rows(),1);
            for (int i=0;i<selectedFaces.size();i++)
                CMesh.row(selectedFaces(i))=selected_face_color();
            set_mesh_colors(CMesh,meshNum);

            //coloring field
            Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->F.rows(),fieldList[meshNum]->N);
            for (int i=0;i<selectedFaces.rows();i++)
                glyphColors.row(selectedFaces(i))=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);

            set_field_colors(glyphColors,meshNum);
        }*/

        //This function assumes vertex-based fields
        /*void inline set_selected_vertices(const Eigen::VectorXi& selectedVertices, const int meshNum=0){
            assert(fieldList[meshNum]->tb->discTangType()==discTangTypeEnum::VERTEX_SPACES);
            std::vector<int> selectedFacesList;
            for (int i=0;i<selectedVertices.size();i++)
                for (int j=0;j<meshList[meshNum]->vertexValence(selectedVertices(i))-(meshList[meshNum]->isBoundaryVertex(selectedVertices(i)) ? 1 : 0);j++)
                    selectedFacesList.push_back(meshList[meshNum]->VF(selectedVertices(i),j));

            Eigen::VectorXi selectedFaces(selectedFacesList.size());
            selectedFaces=Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(selectedFacesList.data(), selectedFacesList.size());

            Eigen::MatrixXd CMesh=directional::DirectionalViewer::default_mesh_color().replicate(meshList[meshNum]->F.rows(),1);
            for (int i=0;i<selectedFaces.size();i++)
                CMesh.row(selectedFaces(i))=selected_face_color();
            set_mesh_colors(CMesh,meshNum);

            //coloring field
            Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->V.rows(),fieldList[meshNum]->N);
            for (int i=0;i<selectedVertices.rows();i++)
                glyphColors.row(selectedVertices(i))=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);

            set_field_colors(glyphColors,meshNum);
        }*/

        /*void inline set_selected_vector(const int selectedFace, const int selectedVector, const int meshNum=0)
        {
            Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->F.rows(),fieldList[meshNum]->N);
            glyphColors.row(selectedFace)=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);
            glyphColors.block(selectedFace,3*selectedVector,1,3)=directional::DirectionalViewer::selected_vector_glyph_color();
            set_field_colors(glyphColors, meshNum);
        }*/

        void inline set_field(const CartesianField& _field,
                              const std::string fieldName = "field",
                              const int meshNum=0,
                              const int fieldNum=0,
                              const double sizeRatio = 0.3,
                              const int sparsity=0,
                              const double offsetRatio = 0.2)

        {
            if (fieldList.size()<fieldNum+1) {
                fieldList.resize(fieldNum + 1);
                psFieldSourceList.resize(fieldNum + 1);
                psFieldList.resize(fieldNum+1);
                psSingList.resize(fieldNum+1);
            }

            fieldList[fieldNum]=&_field;

            psFieldSourceList[fieldNum] = polyscope::registerPointCloud("sources" + std::to_string(fieldNum), _field.tb->sources);
            psFieldList[fieldNum].resize(_field.N);
            psFieldSourceList[fieldNum]->setPointRadius(10e-6);
            psFieldSourceList[fieldNum]->setPointRenderMode(polyscope::PointRenderMode::Quad);
            for (int i=0;i<_field.N;i++) {
                psFieldList[fieldNum][i] = psFieldSourceList[fieldNum]->addVectorQuantity(std::string("field ") + std::to_string(fieldNum) + std::string("-") + std::to_string(i),
                                                                                          _field.extField.block(0, 3 * i, _field.extField.rows(),
                                                                                                                3));
                psFieldList[fieldNum][i]->setVectorLengthScale(sizeRatio*meshList[meshNum]->avgEdgeLength, false);
            }

            set_singularities(fieldList[fieldNum]->singLocalCycles,
                              fieldList[fieldNum]->singIndices,
                              fieldNum);
        }

        /*void inline set_field_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                                         const int meshNum=0,
                                         const double sizeRatio = 0.9,
                                         const int sparsity=0)
        {
            if (fieldColors.size()<meshNum+1)
                fieldColors.resize(meshNum+1);
            if (fieldList.size()<meshNum+1)
                fieldList.resize(meshNum+1);
            fieldColors[meshNum]=C;
            if (C.rows()==0)
                fieldColors[meshNum]=default_glyph_color();

            //directional::glyph_lines_mesh(meshList[meshNum]->F, N[meshNum], fieldColors, CField);

            //TODO: something more efficient than feeding the entire field again
            Eigen::MatrixXd VField, CField;
            Eigen::MatrixXi FField;
            directional::glyph_lines_mesh(fieldList[meshNum]->tb->sources, fieldList[meshNum]->tb->normals, fieldList[meshNum]->tb->adjSpaces, fieldList[meshNum]->extField, fieldColors[meshNum], sizeRatio,meshList[meshNum]->avgEdgeLength, VField, FField, CField, sparsity);

            data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_mesh(VField,FField);
            data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_colors(CField);
            data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].show_lines=false;
        }*/


        void inline set_singularities(const Eigen::VectorXi& singElements,
                                      const Eigen::VectorXi& singIndices,
                                      const int fieldNum=0,
                                      const double radiusRatio=1.25)
        {

            Eigen::MatrixXd singSources(singElements.rows(),3);
            for (int i=0;i<singElements.size();i++)
                singSources.row(i) = fieldList[fieldNum]->tb->cycleSources.row(singElements(i));

            psSingList[fieldNum] = polyscope::registerPointCloud("sings" + std::to_string(fieldNum), singSources);
            psSingList[fieldNum]->addScalarQuantity("indices", singIndices.cast<double>());
        }

        void inline set_seams(const Eigen::VectorXi& combedMatching,
                              const int meshNum=0,
                              const double widthRatio = 0.05)
        {


            Eigen::MatrixXd VSeams1, CSeams1,VSeams2, CSeams2, VSeams, CSeams;
            Eigen::MatrixXi FSeams1, FSeams2, FSeams;

            Eigen::MatrixXi hlHalfedges = Eigen::MatrixXi::Constant(meshList[meshNum]->F.rows(),3,-1);

            std::set<int> seamVertexSet;

            //figuring out the highlighted halfedges
            for (int i=0;i<combedMatching.size();i++){
                if (combedMatching(i)==0)
                    continue;

                int inFaceIndex=-1;

                if (meshList[meshNum]->EF(i,0)!=-1){
                    for (int j=0;j<3;j++)
                        if (meshList[meshNum]->FE(meshList[meshNum]->EF(i,0),j)==i)
                            inFaceIndex=j;
                    hlHalfedges(meshList[meshNum]->EF(i,0), inFaceIndex)=0;
                    seamVertexSet.insert(meshList[meshNum]->F(meshList[meshNum]->EF(i,0), inFaceIndex));
                }
                inFaceIndex=-1;
                if (meshList[meshNum]->EF(i,1)!=-1){
                    for (int j=0;j<3;j++)
                        if (meshList[meshNum]->FE(meshList[meshNum]->EF(i,1),j)==i)
                            inFaceIndex=j;
                    hlHalfedges(meshList[meshNum]->EF(i,1), inFaceIndex)=0;
                    seamVertexSet.insert(meshList[meshNum]->F(meshList[meshNum]->EF(i,1), inFaceIndex));
                }
            }

            Eigen::VectorXi seamVertices(seamVertexSet.size());
            int currPos=0;
            for (std::set<int>::iterator si=seamVertexSet.begin();si!=seamVertexSet.end();si++)
                seamVertices(currPos++)=*si;
            //directional::halfedge_highlights(meshList[meshNum]->V, meshList[meshNum]->F, hlHalfedges, default_seam_color(),VSeams1,FSeams1,CSeams1, widthRatio, 1e-4);
            //directional::vertex_highlights(meshList[meshNum]->V, meshList[meshNum]->F, seamVertices, default_seam_color().replicate(seamVertexSet.size(),1),VSeams2,FSeams2,CSeams2, widthRatio, 1e-4);

            //uniting both meshes
            VSeams.resize(VSeams1.rows()+VSeams2.rows(),3);
            VSeams.topRows(VSeams1.rows())=VSeams1;
            VSeams.bottomRows(VSeams2.rows())=VSeams2;
            FSeams.resize(FSeams1.rows()+FSeams2.rows(),3);
            FSeams.topRows(FSeams1.rows())=FSeams1;
            FSeams.bottomRows(FSeams2.rows())=FSeams2.array()+VSeams1.rows();
            CSeams.resize(CSeams1.rows()+CSeams2.rows(),3);
            CSeams.topRows(CSeams1.rows())=CSeams1;
            CSeams.bottomRows(CSeams2.rows())=CSeams2;

            //data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].clear();
            //data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].set_mesh(VSeams, FSeams);
            //data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].set_colors(CSeams);
            //data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].show_lines = false;
        }



        void inline init_streamlines(const int meshNum=0,
                                     const Eigen::VectorXi& seedLocations=Eigen::VectorXi(),
                                     const double distRatio=3.0)
        {
            if (slData.size()<meshNum+1){
                slData.resize(meshNum+1);
                slState.resize(meshNum+1);
                psStreamlineList.resize(meshNum+1);
            }
            //assert(fieldList[meshNum]->tb->discTangType()==discTangTypeEnum::FACE_SPACES);
            directional::streamlines_init(*fieldList[meshNum], seedLocations,distRatio,slData[meshNum], slState[meshNum]);
        }

        void inline advance_streamlines(const double dTimeRatio,
                                            const int meshNum=0,
                                            const double widthRatio=0.05,
                                            const double colorAttenuationRate = 0.9){

            //double avgEdgeLength = igl::avg_edge_length(meshList[meshNum]->V, meshList[meshNum]->F);  //inefficient!
            double dTime = dTimeRatio*meshList[meshNum]->avgEdgeLength;
            directional::streamlines_next(slData[meshNum], slState[meshNum],dTime);
            double width = widthRatio*meshList[meshNum]->avgEdgeLength;

            //generating colors according to original elements and their time signature
            /*Eigen::MatrixXd slColors(slState[meshNum].segStart.size(),3);

            //problem: if the field is vertex-faced, "orig face" is invalid!
            for (int i=0;i<slState[meshNum].segOrigFace.size();i++){
                if (fieldColors[meshNum].rows()==1)
                    slColors.row(i)=fieldColors[meshNum];
                else{
                    double blendFactor = pow(colorAttenuationRate,(double)slState[meshNum].segTimeSignatures[i]/meshList[meshNum]->avgEdgeLength);
                    //HACK: currently not supporting different colors for vertex-based fields
                    if(fieldList[meshNum]->tb->discTangType()==discTangTypeEnum::FACE_SPACES)
                        slColors.row(i)=fieldColors[meshNum].block(slState[meshNum].segOrigFace[i], 3*slState[meshNum].segOrigVector[i], 1,3);
                    else
                        slColors.row(i)=fieldColors[meshNum].block(slState[meshNum].segOrigFace[0], 3*slState[meshNum].segOrigVector[i], 1,3);
                    slColors.row(i).array()=slColors.row(i).array()*blendFactor+default_mesh_color().array()*(1.0-blendFactor);
                }
            }*/

            Eigen::MatrixXd VStream, CStream;
            Eigen::MatrixXi FStream;
            Eigen::MatrixXd P1(slState[meshNum].segStart.size(),3), P2(slState[meshNum].segEnd.size(),3);
            for (int i=0;i<slState[meshNum].segStart.size();i++){
                P1.row(i)=slState[meshNum].segStart[i];
                P2.row(i)=slState[meshNum].segEnd[i];
            }
            Eigen::MatrixXd nodes(P1.rows()+P2.rows(),3);
            nodes<<P1, P2;
            Eigen::MatrixXi edges(P1.rows(),2);
            edges.col(0) = Eigen::VectorXi::LinSpaced(P1.rows(), 0, P1.rows()-1);
            edges.col(1) = Eigen::VectorXi::LinSpaced(P2.rows(), P1.rows(), P1.rows()+P2.rows()-1);
            psStreamlineList[meshNum] = polyscope::registerCurveNetwork("streamlines" + std::to_string(meshNum), nodes, edges);
            /*directional::line_cylinders(P1,P2, width, slColors, 4, VStream, FStream, CStream);
            data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].clear();
            data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].set_mesh(VStream, FStream);
            data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].set_colors(CStream);
            data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].show_lines = false;*/



        }

        /*void inline set_isolines(const directional::TriMesh& cutMesh,
                                     const Eigen::MatrixXd& vertexFunction,
                                     const int meshNum=0,
                                     const double sizeRatio=0.1)
        {


            Eigen::MatrixXd isoV, isoN;
            Eigen::MatrixXi isoE, isoOrigE;
            Eigen::VectorXi funcNum;

            directional::branched_isolines(cutMesh.V, cutMesh.F, vertexFunction, isoV, isoE, isoOrigE, isoN, funcNum);

            double l = sizeRatio*cutMesh.avgEdgeLength;

            Eigen::MatrixXd VIso, CIso;
            Eigen::MatrixXi FIso;

            Eigen::MatrixXd funcColors = isoline_colors();
            Eigen::MatrixXd CFunc;
            CFunc.resize(funcNum.size(),3);
            for (int i=0;i<funcNum.size();i++)
                CFunc.row(i)=funcColors.row(funcNum(i));

            directional::bar_chains(cutMesh.V, cutMesh.F, isoV,isoE,isoOrigE, isoN,l,(funcNum.template cast<double>().array()+1.0)*l/1000.0,CFunc, VIso, FIso, CIso);

            data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].clear();
            data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].set_mesh(VIso, FIso);
            data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].set_colors(CIso);
            data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_lines = false;
        }

        void inline set_uv(const Eigen::MatrixXd UV,
                               const int meshNum=0)
        {
            data_list[NUMBER_OF_SUBMESHES*meshNum].set_uv(UV);
        }

        void inline set_texture(const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
                                    const int meshNum=0)
        {
            data_list[NUMBER_OF_SUBMESHES*meshNum].set_texture(R,G,B);
        }

        void inline set_active(const bool active, const int meshNum=0){
            for (int i=NUMBER_OF_SUBMESHES*meshNum;i<NUMBER_OF_SUBMESHES*meshNum+NUMBER_OF_SUBMESHES;i++)
                data_list[i].show_faces=active;
        }*/

        void inline toggle_mesh(const bool active, const int meshNum=0){
            psSurfaceMeshList[meshNum]->setEnabled(active);
        }

        void inline toggle_mesh_edges(const bool active, const int meshNum=0){
            //data_list[NUMBER_OF_SUBMESHES*meshNum].show_lines=active;
        }

        void inline toggle_field(const bool active, const int fieldNum=0){
            //data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].show_faces=active;
            for (int i=0;i<psFieldList[fieldNum].size();i++)
                psFieldList[fieldNum][i]->setEnabled(active);
        }

        void inline toggle_singularities(const bool active, const int fieldNum=0){
            //data_list[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].show_faces=active;
            psSingList[fieldNum]->setEnabled(active);
        }

        /*void inline toggle_seams(const bool active, const int meshNum=0){
            data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].show_faces=active;
        }

        void inline toggle_streamlines(const bool active, const int meshNum=0){
            data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].show_faces=active;
        }

        void inline toggle_isolines(const bool active, const int meshNum=0){
            data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_faces=active;
        }

        void inline toggle_texture(const bool active, const int meshNum=0){
            data_list[NUMBER_OF_SUBMESHES*meshNum].show_texture=active;
        }

        //disabling the original mesh
        void inline toggle_edge_data(const bool active, const int meshNum=0){
            data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_faces=active;
            //data_list[NUMBER_OF_SUBMESHES*meshNum].show_faces=!active;
        }

        //static functions for default values
        //Mesh colors
        static Eigen::RowVector3d inline default_mesh_color(){
            return Eigen::RowVector3d::Constant(1.0);
        }

        //Color for faces that are selected for editing and constraints
        static Eigen::RowVector3d inline selected_face_color(){
            return Eigen::RowVector3d(0.7,0.2,0.2);
        }

        //Glyph colors
        static Eigen::RowVector3d inline default_glyph_color(){
            return Eigen::RowVector3d(0.0,0.2,1.0);
        }

        //Glyphs in selected faces
        static Eigen::RowVector3d inline selected_face_glyph_color(){
            return Eigen::RowVector3d(223.0/255.0, 210.0/255.0, 16.0/255.0);
        }

        //The selected glyph currently edited from a selected face
        static Eigen::RowVector3d inline selected_vector_glyph_color(){
            return Eigen::RowVector3d(0.0,1.0,0.5);
        }

        //Colors by indices in each directional object.
        static Eigen::MatrixXd inline isoline_colors(){

            Eigen::Matrix<double, 15,3> glyphPrincipalColors;
            glyphPrincipalColors<< 0.0,0.5,1.0,
                    1.0,0.5,0.0,
                    0.0,1.0,0.5,
                    1.0,0.0,0.5,
                    0.5,0.0,1.0,
                    0.5,1.0,0.0,
                    1.0,0.5,0.5,
                    0.5,1.0,0.5,
                    0.5,0.5,1.0,
                    0.5,1.0,1.0,
                    1.0,0.5,1.0,
                    1.0,1.0,0.5,
                    0.0,1.0,1.0,
                    1.0,0.0,1.0,
                    1.0,1.0,0.0;

            return glyphPrincipalColors;
        }


        //Colors by indices in each directional object. If the field is combed they will appear coherent across faces.
        static Eigen::MatrixXd inline indexed_glyph_colors(const Eigen::MatrixXd& field, bool signSymmetry=true){

            Eigen::Matrix<double, 15,3> glyphPrincipalColors;
            glyphPrincipalColors<< 0.0,0.5,1.0,
                    1.0,0.5,0.0,
                    0.0,1.0,0.5,
                    1.0,0.0,0.5,
                    0.5,0.0,1.0,
                    0.5,1.0,0.0,
                    1.0,0.5,0.5,
                    0.5,1.0,0.5,
                    0.5,0.5,1.0,
                    0.5,1.0,1.0,
                    1.0,0.5,1.0,
                    1.0,1.0,0.5,
                    0.0,1.0,1.0,
                    1.0,0.0,1.0,
                    1.0,1.0,0.0;

            Eigen::MatrixXd fullGlyphColors(field.rows(),field.cols());
            int N = field.cols()/3;
            for (int i=0;i<field.rows();i++)
                for (int j=0;j<N;j++)
                    fullGlyphColors.block(i,3*j,1,3)<< (signSymmetry && (N%2==0) ? glyphPrincipalColors.row(j%(N/2)) : glyphPrincipalColors.row(j));

            return fullGlyphColors;
        }

        //Jet-based singularity colors
        static Eigen::MatrixXd inline default_singularity_colors(const int N){
            Eigen::MatrixXd fullColors;
            Eigen::VectorXd NList(2*N);
            for (int i=0;i<N;i++){
                NList(i)=-N+i;
                NList(N+i)=i+1;
            }
            igl::jet(-NList,true,fullColors);
            return fullColors;
        }

        //Colors for emphasized edges, mostly seams and cuts
        static Eigen::RowVector3d inline default_seam_color(){
            return Eigen::RowVector3d(0.0,0.0,0.0);
        }

        static Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> inline default_texture(){
            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R,texture_G,texture_B;
            unsigned size = 128;
            unsigned size2 = size/2;
            unsigned lineWidth = 5;
            texture_B.setConstant(size, size, 0);
            texture_G.setConstant(size, size, 0);
            texture_R.setConstant(size, size, 0);
            for (unsigned i=0; i<size; ++i)
                for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
                    texture_B(i,j) = texture_G(i,j) = texture_R(i,j) = 255;
            for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
                for (unsigned j=0; j<size; ++j)
                    texture_B(i,j) = texture_G(i,j) = texture_R(i,j) = 255;

            Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> fullTexMat(size*3,size);
            fullTexMat<<texture_R, texture_G, texture_B;
            return fullTexMat;
        }*/

    };  //of DirectionalViewer class

}


#endif