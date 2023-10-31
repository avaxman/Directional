// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_VIEWER_H
#define DIRECTIONAL_VIEWER_H

#include <Eigen/Core>
#include <directional/OpenGLViewer.h>
#include <directional/glyph_lines_mesh.h>
#include <directional/singularity_spheres.h>
#include <directional/seam_lines.h>
#include <directional/edge_diamond_mesh.h>
#include <directional/branched_isolines.h>
#include <directional/bar_chain.h>
#include <directional/halfedge_highlights.h>
#include <directional/vertex_highlights.h>
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

#define NUMBER_OF_SUBMESHES 7
#define FIELD_MESH 1
#define SING_MESH 2
#define SEAMS_MESH 3
#define STREAMLINE_MESH 4
#define EDGE_DIAMOND_MESH 5
#define ISOLINES_MESH 6

    static const char* defaultVertexShaderText =
            "#version 110\n"
            "uniform mat4 MVP;\n"
            "attribute vec3 vCol;\n"
            "attribute vec2 vPos;\n"
            "varying vec3 color;\n"
            "void main()\n"
            "{\n"
            "    gl_Position = MVP * vec4(vPos, 0.0, 1.0);\n"
            "    color = vCol;\n"
            "}\n";

    static const char* defaultFragmentShaderText =
            "#version 110\n"
            "varying vec3 color;\n"
            "void main()\n"
            "{\n"
            "    gl_FragColor = vec4(color, 1.0);\n"
            "}\n";


    struct viewerMeshStatus{
    public:
        typedef enum colorInterpolationType{FACE_BASED, CORNER_BASED, VERTEX_BASED};

        colorInterpolationType ciType;
        bool showEdges;
        bool showFaces;

        Eigen::MatrixXd V;  //vertices
        Eigen::MatrixXi F;  //#Fx3 faces
        Eigen::MatrixXd C;  //colors one row per X={face,corner,vertex} based according to order (corner order is within mesh.F row by row)

        const char* vertexShaderText;
        const char* fragmentShaderText;

        void clear(){
            V=Eigen::MatrixXd(0,3);
            F=Eigen::MatrixXi(0,3);
            C=Eigen::MatrixXd(0,3);
            vertexShaderText=defaultVertexShaderText;
            fragmentShaderText=defaultFragmentShaderText;
        }

        void set_mesh(const Eigen::MatrixXd& _V, const Eigen::MatrixXi& _F){
            V=_V;
            F=_F;
            C=Eigen::MatrixXd(0,3);
            vertexShaderText=defaultVertexShaderText;
            fragmentShaderText=defaultFragmentShaderText;
            ciType=FACE_BASED;
            showEdges=false;
            showFaces=true;
        }

        void set_mesh(const TriMesh& mesh){
            set_mesh(mesh.V, mesh.F);
        }

        viewerMeshStatus(){clear();}
        ~viewerMeshStatus();
    };

    class DirectionalViewer{
    private:
        OpenGLViewer glViewer;
        //Geometry
        std::vector<const TriMesh*> meshList;  //meshes that are being viewed

        std::vector<const CartesianField*> fieldList;
        std::vector<Eigen::MatrixXd> fieldColors;
        std::vector<directional::StreamlineData> slData;
        std::vector<directional::StreamlineState> slState;

        std::vector<Eigen::MatrixXd> edgeVList;  //edge-diamond vertices list
        std::vector<Eigen::MatrixXi> edgeFList;  //edge-diamond faces list
        std::vector<Eigen::VectorXi> edgeFEList;  //edge-diamond faces->original mesh edges list

        std::vector<Eigen::MatrixXd> fieldVList;
        std::vector<Eigen::MatrixXi> fieldFList;

        //Handling the different meshes
        int currViewMesh;
        std::list<viewerMeshStatus> vmList;


    public:
        DirectionalViewer(){
            currViewMesh=0;
            glViewer.init(defaultVertexShaderText, defaultFragmentShaderText);
        }
        ~DirectionalViewer(){
            glViewer.terminate();
        }

        bool launch(){
            glViewer.launch();
        }

        void inline set_mesh(const TriMesh& mesh,
                                 const int meshNum=0)
        {
            Eigen::MatrixXd meshColors;
            meshColors=default_mesh_color();

            if (NUMBER_OF_SUBMESHES*(meshNum+1)>vmList.size()){  //allocating until there
                int currDLSize=vmList.size();
                for (int i=currDLSize;i<NUMBER_OF_SUBMESHES*(meshNum+1);i++)
                    vmList.push_back(viewerMeshStatus());
            }

            currViewMesh=NUMBER_OF_SUBMESHES*meshNum;  //the last triangle mesh
            vmList[NUMBER_OF_SUBMESHES*meshNum].clear();
            vmList[NUMBER_OF_SUBMESHES*meshNum].set_mesh(mesh);
            vmList[NUMBER_OF_SUBMESHES*meshNum].C=meshColors;

            if (meshList.size()<meshNum+1){
                meshList.resize(meshNum+1);
                fieldList.resize(meshNum+1);
                edgeVList.resize(meshNum+1);
                edgeFList.resize(meshNum+1);
                edgeFEList.resize(meshNum+1);
            }
            meshList[meshNum]=&mesh;
        }

        void inline set_mesh_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                                        const int meshNum=0)
        {
            Eigen::MatrixXd meshColors;
            if (C.rows()==0)
                meshColors=default_mesh_color();
            else
                meshColors=C;

            vmList[NUMBER_OF_SUBMESHES*meshNum].ciType=viewerMeshStatus::colorInterpolationType::FACE_BASED;
            vmList[NUMBER_OF_SUBMESHES*meshNum].C=meshColors;
            if (edgeVList[meshNum].size()!=0){
                vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].ciType =viewerMeshStatus::colorInterpolationType::VERTEX_BASED;
                vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_edges=false;
            }
            currViewMesh=NUMBER_OF_SUBMESHES*meshNum;
            //CList[meshNum]=C;
        }

        void inline set_vertex_data(const Eigen::VectorXd& vertexData,
                                    const double minRange,
                                    const double maxRange,
                                    const int meshNum=0)
        {
            Eigen::MatrixXd C;
            parula(vertexData, minRange,maxRange, C);
            set_mesh_colors(C, meshNum);
        }

        //STOPPED HERE

        void inline set_face_data(const Eigen::VectorXd& faceData,
                                  const double minRange,
                                  const double maxRange,
                                  const int meshNum=0)
        {
            Eigen::MatrixXd C;
            parula(faceData, minRange,maxRange, C);
            set_mesh_colors(C, meshNum);
        }

        void inline set_edge_data(const Eigen::VectorXd& edgeData,
                                      const double minRange,
                                      const double maxRange,
                                      const int meshNum=0)
        {

            if (edgeVList[meshNum].size()==0){  //allocate
                edge_diamond_mesh(meshList[meshNum]->V,meshList[meshNum]->F,meshList[meshNum]->EV,meshList[meshNum]->EF,edgeVList[meshNum],edgeFList[meshNum],edgeFEList[meshNum]);
                vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].clear();
                vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].set_mesh(edgeVList[meshNum],edgeFList[meshNum]);
            }

            Eigen::VectorXd edgeFData(edgeFList[meshNum].rows());
            for (int i=0;i<edgeFList[meshNum].rows();i++)
                edgeFData(i)=edgeData(edgeFEList[meshNum](i));

            Eigen::MatrixXd C;
            parula(edgeFData, minRange,maxRange, C);
            vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].set_colors(C);
            vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_faces=false;
            vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_lines=false;

           currViewMesh=NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH;
        }

        //this function assumes face-based fields
        void inline set_selected_faces(const Eigen::VectorXi& selectedFaces, const int meshNum=0){
            assert(fieldList[meshNum]->tb->discTangType()==discTangTypeEnum::FACE_SPACES);
            Eigen::MatrixXd CMesh=directional::DirectionalViewer::default_mesh_color().replicate(meshList[meshNum]->F.rows(),1);
            for (int i=0;i<selectedFaces.size();i++)
                CMesh.row(selectedFaces(i))=selected_face_color();
            vmList[NUMBER_OF_SUBMESHES*meshNum].set_mesh_colors(CMesh);

            //coloring field
            Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->F.rows(),fieldList[meshNum]->N);
            for (int i=0;i<selectedFaces.rows();i++)
                glyphColors.row(selectedFaces(i))=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);

            set_field_colors(glyphColors,meshNum);
        }

        //This function assumes vertex-based fields
        void inline set_selected_vertices(const Eigen::VectorXi& selectedVertices, const int meshNum=0){
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
        }

        void inline set_selected_vector(const int selectedFace, const int selectedVector, const int meshNum=0)
        {
            Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->F.rows(),fieldList[meshNum]->N);
            glyphColors.row(selectedFace)=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);
            glyphColors.block(selectedFace,3*selectedVector,1,3)=directional::DirectionalViewer::selected_vector_glyph_color();
            set_field_colors(glyphColors, meshNum);
        }

        void inline set_field(const CartesianField& _field,
                                  const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                                  const int meshNum=0,
                                  const double sizeRatio = 0.9,
                                  const int sparsity=0,
                                  const double offsetRatio = 0.2)

        {
            if (fieldList.size()<meshNum+1)
                fieldList.resize(meshNum+1);
            if (fieldColors.size()<meshNum+1)
                fieldColors.resize(meshNum+1);
            fieldList[meshNum]=&_field;
            fieldColors[meshNum]=C;
            if (C.rows()==0)
                fieldColors[meshNum]=default_glyph_color();

            Eigen::MatrixXd VField, CField;
            Eigen::MatrixXi FField;
            /*const Eigen::MatrixXd& sources,
             const Eigen::MatrixXi& normals,
             const Eigen::MatrixXi& adjSpaces,
             const Eigen::MatrixXd& extField,
             const Eigen::MatrixXd &glyphColors,
             const double sizeRatio,
             const double avgScale,
             Eigen::MatrixXd &fieldV,
             Eigen::MatrixXi &fieldF,
             Eigen::MatrixXd &fieldC,

             const int sparsity=0,
             const double offsetRatio = 0.2)*/
            directional::glyph_lines_mesh(fieldList[meshNum]->tb->sources, fieldList[meshNum]->tb->normals, fieldList[meshNum]->tb->adjSpaces, fieldList[meshNum]->extField, fieldColors[meshNum], sizeRatio, meshList[meshNum]->avgEdgeLength, VField, FField, CField, sparsity, offsetRatio);
            vmList[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].clear();
            vmList[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_mesh(VField,FField);
            vmList[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_colors(CField);

            set_singularities(fieldList[meshNum]->singLocalCycles,
                              fieldList[meshNum]->singIndices,
                              meshNum);
        }

        void inline set_field_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
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

            vmList[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_mesh(VField,FField);
            vmList[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_colors(CField);
        }


        void inline set_singularities(const Eigen::VectorXi& singElements,
                                          const Eigen::VectorXi& singIndices,
                                          const int meshNum=0,
                                          const double radiusRatio=1.25)
        {
            Eigen::MatrixXd VSings, CSings;
            Eigen::MatrixXi FSings;
            directional::singularity_spheres(fieldList[meshNum]->tb->cycleSources, fieldList[meshNum]->tb->cycleNormals, fieldList[meshNum]->N, meshList[meshNum]->avgEdgeLength, singElements, singIndices, default_singularity_colors(fieldList[meshNum]->N), VSings, FSings, CSings, radiusRatio);
            vmList[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].clear();
            vmList[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].set_mesh(VSings,FSings);
            vmList[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].set_colors(CSings);
            vmList[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].show_lines=false;

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
            directional::halfedge_highlights(meshList[meshNum]->V, meshList[meshNum]->F, hlHalfedges, default_seam_color(),VSeams1,FSeams1,CSeams1, widthRatio, 1e-4);
            directional::vertex_highlights(meshList[meshNum]->V, meshList[meshNum]->F, seamVertices, default_seam_color().replicate(seamVertexSet.size(),1),VSeams2,FSeams2,CSeams2, widthRatio, 1e-4);

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

            vmList[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].clear();
            vmList[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].set_mesh(VSeams, FSeams);
            vmList[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].set_colors(CSeams);
        }



        void inline init_streamlines(const int meshNum=0,
                                         const Eigen::VectorXi& seedLocations=Eigen::VectorXi(),
                                         const double distRatio=3.0)
        {
            if (slData.size()<meshNum+1){
                slData.resize(meshNum+1);
                slState.resize(meshNum+1);
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
            Eigen::MatrixXd slColors(slState[meshNum].segStart.size(),3);

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
            }

            Eigen::MatrixXd VStream, CStream;
            Eigen::MatrixXi FStream;
            Eigen::MatrixXd P1(slState[meshNum].segStart.size(),3), P2(slState[meshNum].segEnd.size(),3);
            for (int i=0;i<slState[meshNum].segStart.size();i++){
                P1.row(i)=slState[meshNum].segStart[i];
                P2.row(i)=slState[meshNum].segEnd[i];
            }
            directional::line_cylinders(P1,P2, width, slColors, 4, VStream, FStream, CStream);
            vmList[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].clear();
            vmList[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].set_mesh(VStream, FStream);
            vmList[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].set_colors(CStream);

        }

        void inline set_isolines(const directional::TriMesh& cutMesh,
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

            vmList[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].clear();
            vmList[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].set_mesh(VIso, FIso);
            vmList[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].set_colors(CIso);
            vmList[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_lines = false;
        }

        void inline set_uv(const Eigen::MatrixXd UV,
                               const int meshNum=0)
        {
            vmList[NUMBER_OF_SUBMESHES*meshNum].set_uv(UV);
        }

        void inline set_texture(const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                                    const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
                                    const int meshNum=0)
        {
            vmList[NUMBER_OF_SUBMESHES*meshNum].set_texture(R,G,B);
        }

        void inline set_active(const bool active, const int meshNum=0){
            for (int i=NUMBER_OF_SUBMESHES*meshNum;i<NUMBER_OF_SUBMESHES*meshNum+NUMBER_OF_SUBMESHES;i++)
                vmList[i].show_faces=active;
        }

        void inline toggle_mesh(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum].show_faces=active;
        }

        void inline toggle_mesh_edges(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum].show_lines=active;
        }

        void inline toggle_field(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].show_faces=active;
        }

        void inline toggle_singularities(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].show_faces=active;
        }

        void inline toggle_seams(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].show_faces=active;
        }

        void inline toggle_streamlines(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].show_faces=active;
        }

        void inline toggle_isolines(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_faces=active;
        }

        void inline toggle_texture(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum].show_texture=active;
        }

        //disabling the original mesh
        void inline toggle_edge_data(const bool active, const int meshNum=0){
            vmList[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_faces=active;
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
            jet(-NList,true,fullColors);
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
        }

        static Eigen::MatrixXd inline Viridis()
        {
            static double viridis_color_map[256][3] = {
                    { 0.267004, 0.004874, 0.329415 },
                    { 0.268510, 0.009605, 0.335427 },
                    { 0.269944, 0.014625, 0.341379 },
                    { 0.271305, 0.019942, 0.347269 },
                    { 0.272594, 0.025563, 0.353093 },
                    { 0.273809, 0.031497, 0.358853 },
                    { 0.274952, 0.037752, 0.364543 },
                    { 0.276022, 0.044167, 0.370164 },
                    { 0.277018, 0.050344, 0.375715 },
                    { 0.277941, 0.056324, 0.381191 },
                    { 0.278791, 0.062145, 0.386592 },
                    { 0.279566, 0.067836, 0.391917 },
                    { 0.280267, 0.073417, 0.397163 },
                    { 0.280894, 0.078907, 0.402329 },
                    { 0.281446, 0.084320, 0.407414 },
                    { 0.281924, 0.089666, 0.412415 },
                    { 0.282327, 0.094955, 0.417331 },
                    { 0.282656, 0.100196, 0.422160 },
                    { 0.282910, 0.105393, 0.426902 },
                    { 0.283091, 0.110553, 0.431554 },
                    { 0.283197, 0.115680, 0.436115 },
                    { 0.283229, 0.120777, 0.440584 },
                    { 0.283187, 0.125848, 0.444960 },
                    { 0.283072, 0.130895, 0.449241 },
                    { 0.282884, 0.135920, 0.453427 },
                    { 0.282623, 0.140926, 0.457517 },
                    { 0.282290, 0.145912, 0.461510 },
                    { 0.281887, 0.150881, 0.465405 },
                    { 0.281412, 0.155834, 0.469201 },
                    { 0.280868, 0.160771, 0.472899 },
                    { 0.280255, 0.165693, 0.476498 },
                    { 0.279574, 0.170599, 0.479997 },
                    { 0.278826, 0.175490, 0.483397 },
                    { 0.278012, 0.180367, 0.486697 },
                    { 0.277134, 0.185228, 0.489898 },
                    { 0.276194, 0.190074, 0.493001 },
                    { 0.275191, 0.194905, 0.496005 },
                    { 0.274128, 0.199721, 0.498911 },
                    { 0.273006, 0.204520, 0.501721 },
                    { 0.271828, 0.209303, 0.504434 },
                    { 0.270595, 0.214069, 0.507052 },
                    { 0.269308, 0.218818, 0.509577 },
                    { 0.267968, 0.223549, 0.512008 },
                    { 0.266580, 0.228262, 0.514349 },
                    { 0.265145, 0.232956, 0.516599 },
                    { 0.263663, 0.237631, 0.518762 },
                    { 0.262138, 0.242286, 0.520837 },
                    { 0.260571, 0.246922, 0.522828 },
                    { 0.258965, 0.251537, 0.524736 },
                    { 0.257322, 0.256130, 0.526563 },
                    { 0.255645, 0.260703, 0.528312 },
                    { 0.253935, 0.265254, 0.529983 },
                    { 0.252194, 0.269783, 0.531579 },
                    { 0.250425, 0.274290, 0.533103 },
                    { 0.248629, 0.278775, 0.534556 },
                    { 0.246811, 0.283237, 0.535941 },
                    { 0.244972, 0.287675, 0.537260 },
                    { 0.243113, 0.292092, 0.538516 },
                    { 0.241237, 0.296485, 0.539709 },
                    { 0.239346, 0.300855, 0.540844 },
                    { 0.237441, 0.305202, 0.541921 },
                    { 0.235526, 0.309527, 0.542944 },
                    { 0.233603, 0.313828, 0.543914 },
                    { 0.231674, 0.318106, 0.544834 },
                    { 0.229739, 0.322361, 0.545706 },
                    { 0.227802, 0.326594, 0.546532 },
                    { 0.225863, 0.330805, 0.547314 },
                    { 0.223925, 0.334994, 0.548053 },
                    { 0.221989, 0.339161, 0.548752 },
                    { 0.220057, 0.343307, 0.549413 },
                    { 0.218130, 0.347432, 0.550038 },
                    { 0.216210, 0.351535, 0.550627 },
                    { 0.214298, 0.355619, 0.551184 },
                    { 0.212395, 0.359683, 0.551710 },
                    { 0.210503, 0.363727, 0.552206 },
                    { 0.208623, 0.367752, 0.552675 },
                    { 0.206756, 0.371758, 0.553117 },
                    { 0.204903, 0.375746, 0.553533 },
                    { 0.203063, 0.379716, 0.553925 },
                    { 0.201239, 0.383670, 0.554294 },
                    { 0.199430, 0.387607, 0.554642 },
                    { 0.197636, 0.391528, 0.554969 },
                    { 0.195860, 0.395433, 0.555276 },
                    { 0.194100, 0.399323, 0.555565 },
                    { 0.192357, 0.403199, 0.555836 },
                    { 0.190631, 0.407061, 0.556089 },
                    { 0.188923, 0.410910, 0.556326 },
                    { 0.187231, 0.414746, 0.556547 },
                    { 0.185556, 0.418570, 0.556753 },
                    { 0.183898, 0.422383, 0.556944 },
                    { 0.182256, 0.426184, 0.557120 },
                    { 0.180629, 0.429975, 0.557282 },
                    { 0.179019, 0.433756, 0.557430 },
                    { 0.177423, 0.437527, 0.557565 },
                    { 0.175841, 0.441290, 0.557685 },
                    { 0.174274, 0.445044, 0.557792 },
                    { 0.172719, 0.448791, 0.557885 },
                    { 0.171176, 0.452530, 0.557965 },
                    { 0.169646, 0.456262, 0.558030 },
                    { 0.168126, 0.459988, 0.558082 },
                    { 0.166617, 0.463708, 0.558119 },
                    { 0.165117, 0.467423, 0.558141 },
                    { 0.163625, 0.471133, 0.558148 },
                    { 0.162142, 0.474838, 0.558140 },
                    { 0.160665, 0.478540, 0.558115 },
                    { 0.159194, 0.482237, 0.558073 },
                    { 0.157729, 0.485932, 0.558013 },
                    { 0.156270, 0.489624, 0.557936 },
                    { 0.154815, 0.493313, 0.557840 },
                    { 0.153364, 0.497000, 0.557724 },
                    { 0.151918, 0.500685, 0.557587 },
                    { 0.150476, 0.504369, 0.557430 },
                    { 0.149039, 0.508051, 0.557250 },
                    { 0.147607, 0.511733, 0.557049 },
                    { 0.146180, 0.515413, 0.556823 },
                    { 0.144759, 0.519093, 0.556572 },
                    { 0.143343, 0.522773, 0.556295 },
                    { 0.141935, 0.526453, 0.555991 },
                    { 0.140536, 0.530132, 0.555659 },
                    { 0.139147, 0.533812, 0.555298 },
                    { 0.137770, 0.537492, 0.554906 },
                    { 0.136408, 0.541173, 0.554483 },
                    { 0.135066, 0.544853, 0.554029 },
                    { 0.133743, 0.548535, 0.553541 },
                    { 0.132444, 0.552216, 0.553018 },
                    { 0.131172, 0.555899, 0.552459 },
                    { 0.129933, 0.559582, 0.551864 },
                    { 0.128729, 0.563265, 0.551229 },
                    { 0.127568, 0.566949, 0.550556 },
                    { 0.126453, 0.570633, 0.549841 },
                    { 0.125394, 0.574318, 0.549086 },
                    { 0.124395, 0.578002, 0.548287 },
                    { 0.123463, 0.581687, 0.547445 },
                    { 0.122606, 0.585371, 0.546557 },
                    { 0.121831, 0.589055, 0.545623 },
                    { 0.121148, 0.592739, 0.544641 },
                    { 0.120565, 0.596422, 0.543611 },
                    { 0.120092, 0.600104, 0.542530 },
                    { 0.119738, 0.603785, 0.541400 },
                    { 0.119512, 0.607464, 0.540218 },
                    { 0.119423, 0.611141, 0.538982 },
                    { 0.119483, 0.614817, 0.537692 },
                    { 0.119699, 0.618490, 0.536347 },
                    { 0.120081, 0.622161, 0.534946 },
                    { 0.120638, 0.625828, 0.533488 },
                    { 0.121380, 0.629492, 0.531973 },
                    { 0.122312, 0.633153, 0.530398 },
                    { 0.123444, 0.636809, 0.528763 },
                    { 0.124780, 0.640461, 0.527068 },
                    { 0.126326, 0.644107, 0.525311 },
                    { 0.128087, 0.647749, 0.523491 },
                    { 0.130067, 0.651384, 0.521608 },
                    { 0.132268, 0.655014, 0.519661 },
                    { 0.134692, 0.658636, 0.517649 },
                    { 0.137339, 0.662252, 0.515571 },
                    { 0.140210, 0.665859, 0.513427 },
                    { 0.143303, 0.669459, 0.511215 },
                    { 0.146616, 0.673050, 0.508936 },
                    { 0.150148, 0.676631, 0.506589 },
                    { 0.153894, 0.680203, 0.504172 },
                    { 0.157851, 0.683765, 0.501686 },
                    { 0.162016, 0.687316, 0.499129 },
                    { 0.166383, 0.690856, 0.496502 },
                    { 0.170948, 0.694384, 0.493803 },
                    { 0.175707, 0.697900, 0.491033 },
                    { 0.180653, 0.701402, 0.488189 },
                    { 0.185783, 0.704891, 0.485273 },
                    { 0.191090, 0.708366, 0.482284 },
                    { 0.196571, 0.711827, 0.479221 },
                    { 0.202219, 0.715272, 0.476084 },
                    { 0.208030, 0.718701, 0.472873 },
                    { 0.214000, 0.722114, 0.469588 },
                    { 0.220124, 0.725509, 0.466226 },
                    { 0.226397, 0.728888, 0.462789 },
                    { 0.232815, 0.732247, 0.459277 },
                    { 0.239374, 0.735588, 0.455688 },
                    { 0.246070, 0.738910, 0.452024 },
                    { 0.252899, 0.742211, 0.448284 },
                    { 0.259857, 0.745492, 0.444467 },
                    { 0.266941, 0.748751, 0.440573 },
                    { 0.274149, 0.751988, 0.436601 },
                    { 0.281477, 0.755203, 0.432552 },
                    { 0.288921, 0.758394, 0.428426 },
                    { 0.296479, 0.761561, 0.424223 },
                    { 0.304148, 0.764704, 0.419943 },
                    { 0.311925, 0.767822, 0.415586 },
                    { 0.319809, 0.770914, 0.411152 },
                    { 0.327796, 0.773980, 0.406640 },
                    { 0.335885, 0.777018, 0.402049 },
                    { 0.344074, 0.780029, 0.397381 },
                    { 0.352360, 0.783011, 0.392636 },
                    { 0.360741, 0.785964, 0.387814 },
                    { 0.369214, 0.788888, 0.382914 },
                    { 0.377779, 0.791781, 0.377939 },
                    { 0.386433, 0.794644, 0.372886 },
                    { 0.395174, 0.797475, 0.367757 },
                    { 0.404001, 0.800275, 0.362552 },
                    { 0.412913, 0.803041, 0.357269 },
                    { 0.421908, 0.805774, 0.351910 },
                    { 0.430983, 0.808473, 0.346476 },
                    { 0.440137, 0.811138, 0.340967 },
                    { 0.449368, 0.813768, 0.335384 },
                    { 0.458674, 0.816363, 0.329727 },
                    { 0.468053, 0.818921, 0.323998 },
                    { 0.477504, 0.821444, 0.318195 },
                    { 0.487026, 0.823929, 0.312321 },
                    { 0.496615, 0.826376, 0.306377 },
                    { 0.506271, 0.828786, 0.300362 },
                    { 0.515992, 0.831158, 0.294279 },
                    { 0.525776, 0.833491, 0.288127 },
                    { 0.535621, 0.835785, 0.281908 },
                    { 0.545524, 0.838039, 0.275626 },
                    { 0.555484, 0.840254, 0.269281 },
                    { 0.565498, 0.842430, 0.262877 },
                    { 0.575563, 0.844566, 0.256415 },
                    { 0.585678, 0.846661, 0.249897 },
                    { 0.595839, 0.848717, 0.243329 },
                    { 0.606045, 0.850733, 0.236712 },
                    { 0.616293, 0.852709, 0.230052 },
                    { 0.626579, 0.854645, 0.223353 },
                    { 0.636902, 0.856542, 0.216620 },
                    { 0.647257, 0.858400, 0.209861 },
                    { 0.657642, 0.860219, 0.203082 },
                    { 0.668054, 0.861999, 0.196293 },
                    { 0.678489, 0.863742, 0.189503 },
                    { 0.688944, 0.865448, 0.182725 },
                    { 0.699415, 0.867117, 0.175971 },
                    { 0.709898, 0.868751, 0.169257 },
                    { 0.720391, 0.870350, 0.162603 },
                    { 0.730889, 0.871916, 0.156029 },
                    { 0.741388, 0.873449, 0.149561 },
                    { 0.751884, 0.874951, 0.143228 },
                    { 0.762373, 0.876424, 0.137064 },
                    { 0.772852, 0.877868, 0.131109 },
                    { 0.783315, 0.879285, 0.125405 },
                    { 0.793760, 0.880678, 0.120005 },
                    { 0.804182, 0.882046, 0.114965 },
                    { 0.814576, 0.883393, 0.110347 },
                    { 0.824940, 0.884720, 0.106217 },
                    { 0.835270, 0.886029, 0.102646 },
                    { 0.845561, 0.887322, 0.099702 },
                    { 0.855810, 0.888601, 0.097452 },
                    { 0.866013, 0.889868, 0.095953 },
                    { 0.876168, 0.891125, 0.095250 },
                    { 0.886271, 0.892374, 0.095374 },
                    { 0.896320, 0.893616, 0.096335 },
                    { 0.906311, 0.894855, 0.098125 },
                    { 0.916242, 0.896091, 0.100717 },
                    { 0.926106, 0.897330, 0.104071 },
                    { 0.935904, 0.898570, 0.108131 },
                    { 0.945636, 0.899815, 0.112838 },
                    { 0.955300, 0.901065, 0.118128 },
                    { 0.964894, 0.902323, 0.123941 },
                    { 0.974417, 0.903590, 0.130215 },
                    { 0.983868, 0.904867, 0.136897 },
                    { 0.993248, 0.906157, 0.143936 }
            };

            /*template <typename T>
            IGL_INLINE void igl::colormap(
                    const double palette[256][3], const T x_in, T & r, T & g, T & b)
            {
                static const unsigned int pal = 256;
                const T zero = 0.0;
                const T one = 1.0;
                T x_in_clamped = static_cast<T>(std::max(zero, std::min(one, x_in)));

                // simple rgb lerp from palette
                unsigned int least = std::floor(x_in_clamped * static_cast<T>(pal - 1));
                unsigned int most = std::ceil(x_in_clamped * static_cast<T>(pal - 1));

                T _r[2] = { static_cast<T>(palette[least][0]), static_cast<T>(palette[most][0]) };
                T _g[2] = { static_cast<T>(palette[least][1]), static_cast<T>(palette[most][1]) };
                T _b[2] = { static_cast<T>(palette[least][2]), static_cast<T>(palette[most][2]) };

                T t = std::max(zero, std::min(one, static_cast<T>(fmod(x_in_clamped * static_cast<T>(pal), one))));

                r = std::max(zero, std::min(one, (one - t) * _r[0] + t * _r[1]));
                g = std::max(zero, std::min(one, (one - t) * _g[0] + t * _g[1]));
                b = std::max(zero, std::min(one, (one - t) * _b[0] + t * _b[1]));
            }

            template <typename DerivedZ, typename DerivedC>
            IGL_INLINE void igl::colormap(
                    const ColorMapType cm,
                    const Eigen::MatrixBase<DerivedZ> & Z,
                    const bool normalize,
                    Eigen::PlainObjectBase<DerivedC> & C)
            {
                const double min_z = normalize ? Z.minCoeff() : 0;
                const double max_z = normalize ? Z.maxCoeff() : 1;
                return colormap(cm, Z, min_z, max_z, C);
            }

            template <typename DerivedZ, typename DerivedC>
            IGL_INLINE void igl::colormap(
                    const ColorMapType cm,
                    const Eigen::MatrixBase<DerivedZ> & Z,
                    const double min_z,
                    const double max_z,
                    Eigen::PlainObjectBase<DerivedC> & C)
            {
                C.resize(Z.rows(),3);
                double denom = (max_z - min_z);
                denom = (denom == 0) ? 1 : denom;
                for(int r = 0; r < Z.rows(); ++r) {
                    colormap(
                            cm,
                            (typename DerivedC::Scalar)((-min_z + Z(r,0)) / denom),
                            C(r,0),
                            C(r,1),
                            C(r,2));
                }*/

            
        }

    };  //of DirectionalViewer class

}


#endif
