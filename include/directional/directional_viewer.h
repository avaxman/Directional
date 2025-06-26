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
#include <directional/bar_chain.h>
#include <directional/branched_isolines.h>

/***
 This class implements the Directional viewer, as an extension of the libigl viewer (as a wrapper). This
 viewer providers specialized functionality for outputting directional fields and their combinatorial and geometic properties.
 The numerous tutorial examples highlight its functionality.
 ***/


namespace directional
{

class DirectionalViewer{
public:
    
    //Directional quantities
    std::vector<const TriMesh*> surfaceMeshList;  //meshes that are being viewed
    std::vector<const CartesianField*> fieldList;
    std::vector<directional::StreamlineData> slData;
    std::vector<directional::StreamlineState> slState;
    
    //Polyscope Quantities
    std::vector<polyscope::SurfaceMesh*> psSurfaceMeshList;
    std::vector<polyscope::PointCloud*> psGlyphSourceList;
    std::vector<std::vector<polyscope::PointCloudVectorQuantity*>> psGlyphList;
    std::vector<polyscope::PointCloud*> psSingList;
    std::vector<polyscope::CurveNetwork*> psStreamlineList;
    std::vector<polyscope::CurveNetwork*> psEdgeHighlightList;
    std::vector<polyscope::CurveNetwork*> psIsolineList;
    
    
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
    
    void inline set_surface_mesh(const TriMesh& mesh,
                                 const int meshNum=0,
                                 const std::string name="")
    {
        if (psSurfaceMeshList.size()<meshNum+1) {
            surfaceMeshList.resize(meshNum + 1);
            psSurfaceMeshList.resize(meshNum + 1);
            //avgEdgeLengths.resize(meshNum+1);
        }
        
        surfaceMeshList[meshNum]=&mesh;
        std::string meshName;
        if (name.empty())
            meshName = "Mesh " + std::to_string(meshNum);
        else
            meshName = name;
        psSurfaceMeshList[meshNum]=polyscope::registerSurfaceMesh(meshName, mesh.V, mesh.F)->setSurfaceColor(default_face_color());
        
        std::vector<int> permArr(mesh.EV.rows());
        for (int i=0;i<permArr.size();i++)
            permArr[i]=i;
        psSurfaceMeshList[meshNum]->setEdgePermutation(permArr);
        
        //avgEdgeLengths[meshNum] = mesh.avgEdgeLength;
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
                                const std::string name="",
                                const int meshNum=0)
    {
        std::string dataName;
        if (name.empty())
            dataName = "Vertex data " + std::to_string(meshNum);
        else
            dataName = name;
        polyscope::SurfaceVertexScalarQuantity* vertexQuantity = psSurfaceMeshList[meshNum]->addVertexScalarQuantity(dataName, vertexData);
        vertexQuantity->setMapRange(std::pair<double, double>(minRange, maxRange));
    }
    
    void inline set_face_data(const Eigen::VectorXd& faceData,
                              const double minRange,
                              const double maxRange,
                              const std::string name="",
                              const int meshNum=0)
    {
        std::string dataName;
        if (name.empty())
            dataName = "Face data " + std::to_string(meshNum);
        else
            dataName = name;
        polyscope::SurfaceFaceScalarQuantity* faceQuantity = psSurfaceMeshList[meshNum]->addFaceScalarQuantity(dataName, faceData);
        faceQuantity->setMapRange(std::pair<double, double>(minRange, maxRange));
    }
    
    void inline set_edge_data(const Eigen::VectorXd& edgeData,
                              const double minRange,
                              const double maxRange,
                              const std::string name="",
                              const int meshNum=0)
    {
        
        std::string dataName;
        if (name.empty())
            dataName = "Edge data " + std::to_string(meshNum);
        else
            dataName = name;
        polyscope::SurfaceEdgeScalarQuantity* edgeQuantity = psSurfaceMeshList[meshNum]->addEdgeScalarQuantity(dataName, edgeData);
        edgeQuantity->setMapRange(std::pair<double, double>(minRange, maxRange));
    }
    
    void inline highlight_faces(const Eigen::VectorXi& selectedFaces,
                                const std::string name="",
                                const int meshNum=0){
        
        glm::vec3 defaultColorglm = psSurfaceMeshList[meshNum]->getSurfaceColor();
        Eigen::RowVector3d defaultColor; for (int i=0;i<3;i++) defaultColor(i)=defaultColorglm[i];
        glm::vec3 highlightColor = highlight_face_color();
        Eigen::MatrixXd faceColors(surfaceMeshList[meshNum]->F.rows(),3);
        faceColors.rowwise() = defaultColor;
        for (int i=0;i<selectedFaces.size();i++)
            faceColors.row(selectedFaces(i))<<highlightColor.x,highlightColor.y,highlightColor.z;
        
        std::string highName;
        if (name.empty())
            highName = "Face highlights " + std::to_string(meshNum);
        else
            highName = name;
        psSurfaceMeshList[meshNum]->addFaceColorQuantity(highName, faceColors)->setEnabled(true);
    }
    
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
    
    void inline set_cartesian_field(const CartesianField& _field,
                                    const std::string name="",
                                    const int fieldNum = 0,
                                    const int sparsity = 0,
                                    double maxLengthRatio = 0.3)
    
    {
        //const bool combedColors = false;
        if (psGlyphSourceList.size()<fieldNum+1) {
            fieldList.resize(fieldNum + 1);
            psGlyphSourceList.resize(fieldNum + 1);
            psGlyphList.resize(fieldNum+1);
            psSingList.resize(fieldNum+1);
        }
        
        fieldList[fieldNum] = &_field;
        
        Eigen::VectorXi sampledSpaces = samples_tangent_bundle(_field.tb->sources,_field.tb->adjSpaces, sparsity);
        Eigen::MatrixXd sampledSources, sampledField;
        Eigen::VectorXi three(3); three<<0,1,2;
        Eigen::VectorXi threeN(3*_field.N);
        for (int i=0;i<_field.N;i++)
            threeN.segment(3*i,3)<<3*i,3*i+1,3*i+2;
        directional::slice(_field.tb->sources, sampledSpaces, three, sampledSources);
        directional::slice(_field.extField, sampledSpaces, threeN, sampledField);
        std::string fieldName;
        if (name.empty())
            fieldName = "Field " + std::to_string(fieldNum);
        else
            fieldName = name;
        psGlyphSourceList[fieldNum] = polyscope::registerPointCloud(fieldName, sampledSources);
        psGlyphList[fieldNum].resize(_field.N);
        psGlyphSourceList[fieldNum]->setPointRadius(10e-6);
        psGlyphSourceList[fieldNum]->setPointRenderMode(polyscope::PointRenderMode::Quad);
        for (int i=0;i<_field.N;i++) {
            psGlyphList[fieldNum][i] = psGlyphSourceList[fieldNum]->addVectorQuantity("Vector " + std::to_string(i),
                                                                                      sampledField.block(0, 3 * i, sampledField.rows(),
                                                                                                         3));
            psGlyphList[fieldNum][i]->setVectorLengthScale(maxLengthRatio*_field.tb->avgAdjLength, false);
            psGlyphList[fieldNum][i]->setEnabled(true);
            //if (!combedColors)
            psGlyphList[fieldNum][i]->setVectorColor(default_glyph_color());
        }
        
        //if (addSingularities)
        set_singularities(_field,  _field.singLocalCycles, _field.singIndices, name, fieldNum);
    }
    
    void inline set_raw_field(const Eigen::MatrixXd& sources,
                              const Eigen::MatrixXd& rawField,
                              const double vectorScale,
                              const std::string name="",
                              const int fieldNum = 0)
    
    {
        const double sizeRatio = 0.3;
        //const bool combedColors = false;
        if (psGlyphSourceList.size()<fieldNum+1) {
            //fieldList.resize(fieldNum + 1);
            psGlyphSourceList.resize(fieldNum + 1);
            psGlyphList.resize(fieldNum+1);
            psSingList.resize(fieldNum+1);
        }
        
        //fieldList[fieldNum] = &_field;
        
        std::string fieldName;
        if (name.empty())
            fieldName = "Field " + std::to_string(fieldNum);
        else
            fieldName = name;
        psGlyphSourceList[fieldNum] = polyscope::registerPointCloud(fieldName, sources);
        psGlyphList[fieldNum].resize(rawField.cols()/3);
        psGlyphSourceList[fieldNum]->setPointRadius(10e-6);
        psGlyphSourceList[fieldNum]->setPointRenderMode(polyscope::PointRenderMode::Quad);
        for (int i=0;i<rawField.cols()/3;i++) {
            psGlyphList[fieldNum][i] = psGlyphSourceList[fieldNum]->addVectorQuantity("Vector " + std::to_string(i),
                                                                                      rawField.block(0, 3 * i, rawField.rows(),
                                                                                                     3));
            psGlyphList[fieldNum][i]->setVectorLengthScale(sizeRatio*vectorScale, false);
            psGlyphList[fieldNum][i]->setEnabled(true);
            //if (!combedColors)
            psGlyphList[fieldNum][i]->setVectorColor(default_glyph_color());
        }
        
        //if (addSingularities)
        //set_singularities(_field,  _field.singLocalCycles, _field.singIndices, name, fieldNum);
    }
    
    Eigen::MatrixXd inline set_1form(const Eigen::VectorXd& oneForm,
                                     const std::string formName = "1-form",
                                     const int meshNum=0,
                                     const int fieldNum=0,
                                     const int samplingRate = 2,
                                     const double baryOffset = 0.2,
                                     const double sizeRatio = 0.1)
    
    {
        if (psGlyphList.size()<fieldNum+1) {
            psGlyphSourceList.resize(fieldNum + 1);
            psGlyphList.resize(fieldNum+1);
        }
        
        //creating barcentric template
        Eigen::VectorXd baryRange = Eigen::VectorXd::LinSpaced(samplingRate, 0.0, 1.0);
        //creating all combinations
        std::vector<Eigen::RowVector3d> baryCoords;
        for (int i=0;i<baryRange.size();i++)
            for (int j=0;j<baryRange.size();j++)
                for (int k=0;k<baryRange.size();k++)
                    if (std::abs(baryRange(i)+baryRange(j)+baryRange(k)-1.0)<10e-7)
                        baryCoords.push_back(Eigen::RowVector3d(baryRange(i), baryRange(j), baryRange(k)));
        
        for (int i=0;i<baryCoords.size();i++){
            baryCoords[i] = Eigen::RowVector3d::Ones()*baryOffset+baryCoords[i]*(1.0-3.0*baryOffset);
            //std::cout<<baryCoords[i]<<std::endl;
        }
        
        Eigen::MatrixXd sources(surfaceMeshList[meshNum]->F.rows()*baryCoords.size(),3);
        Eigen::MatrixXd field(surfaceMeshList[meshNum]->F.rows()*baryCoords.size(),3);
        for (int f=0;f<surfaceMeshList[meshNum]->F.rows();f++) {
            Eigen::RowVector3d n = surfaceMeshList[meshNum]->faceNormals.row(f);
            Eigen::RowVector3d e01 = surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f,1)) - surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f,0));
            Eigen::RowVector3d e12 = surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f,2)) - surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f,1));
            Eigen::RowVector3d e20 = surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f,0)) - surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f,2));
            
            Eigen::RowVector3d gB0 = n.cross(e12)/(2*surfaceMeshList[meshNum]->faceAreas(f));
            Eigen::RowVector3d gB1 = n.cross(e20)/(2*surfaceMeshList[meshNum]->faceAreas(f));
            Eigen::RowVector3d gB2 = n.cross(e01)/(2*surfaceMeshList[meshNum]->faceAreas(f));
            for (int i = 0; i < baryCoords.size(); i++) {
                sources.row(baryCoords.size()*f+i) =surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f, 0)) * baryCoords[i](0) +
                surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f, 1)) *baryCoords[i](1) +
                surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->F(f, 2)) * baryCoords[i](2);
                
                Eigen::RowVector3d whit01 = baryCoords[i](0)*gB1 - baryCoords[i](1)*gB0;
                Eigen::RowVector3d whit12 = baryCoords[i](1)*gB2 - baryCoords[i](2)*gB1;
                Eigen::RowVector3d whit20 = baryCoords[i](2)*gB0 - baryCoords[i](0)*gB2;
                
                field.row(baryCoords.size()*f+i) = oneForm(surfaceMeshList[meshNum]->FE(f,0))*surfaceMeshList[meshNum]->FEs(f,0)*whit01 +
                oneForm(surfaceMeshList[meshNum]->FE(f,1))*surfaceMeshList[meshNum]->FEs(f,1)*whit12 +
                oneForm(surfaceMeshList[meshNum]->FE(f,2))*surfaceMeshList[meshNum]->FEs(f,2)*whit20;
            }
        }
        
        psGlyphSourceList[fieldNum] = polyscope::registerPointCloud("sources" + std::to_string(fieldNum), sources);
        psGlyphList[fieldNum].resize(1);
        psGlyphSourceList[fieldNum]->setPointRadius(10e-6);
        psGlyphSourceList[fieldNum]->setPointRenderMode(polyscope::PointRenderMode::Quad);
        psGlyphList[fieldNum][0] = psGlyphSourceList[fieldNum]->addVectorQuantity(std::string(formName),field);
        psGlyphList[fieldNum][0]->setVectorLengthScale(sizeRatio*surfaceMeshList[meshNum]->avgEdgeLength, false);
        return field;
    }
    
    void inline set_singularities(const CartesianField& field,
                                  const Eigen::VectorXi& singElements,
                                  const Eigen::VectorXi& singIndices,
                                  const std::string fieldName="",
                                  const int fieldNum=0,
                                  const double radiusRatio=0.3)
    {
        assert("set_singularities(): singElements.size()!=singIndices.size()" && singElements.size()==singIndices.size());
        Eigen::MatrixXd singSources(singElements.rows(),3);
        for (int i=0;i<singElements.size();i++)
            singSources.row(i) = field.tb->cycleSources.row(singElements(i));
        
        std::string singName;
        if (fieldName.empty())
            singName = "Singularities " + std::to_string(fieldNum);
        else
            singName = fieldName + " singularities";
        
        psSingList[fieldNum] = polyscope::registerPointCloud(singName, singSources)->setPointRadius(radiusRatio*field.tb->avgAdjLength, false);
        std::vector<glm::vec3> singColors(singIndices.size());
        for (int i=0;i<singIndices.size();i++)
            singColors[i] = default_index_color(singIndices[i]);
        
        psSingList[fieldNum]->addColorQuantity("Indices", singColors)->setEnabled(true);
        //psSingList[fieldNum]->addScalarQuantity("Indices", singIndices.cast<double>())->setColorMap("coolwarm")->setMapRange(std::pair<double, double>(-indexRange,indexRange))->setEnabled(true);
    }
    
    void inline highlight_edges(const Eigen::VectorXi& highlightEdges,
                                const std::string name = "",
                                const int meshNum=0,
                                const double widthRatio = 0.05)
    {
        if (psEdgeHighlightList.size()<meshNum+1)
            psEdgeHighlightList.resize(meshNum+1);
        
        Eigen::MatrixXi seamEdges(highlightEdges.size(), 2);
        Eigen::MatrixXd seamNodes(2*highlightEdges.size(), 3);
        for (int i=0;i<highlightEdges.size();i++){
            seamEdges.row(i) << 2 * i, 2 * i+ 1;
            seamNodes.row(2*i)<<surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->EV(highlightEdges[i],0));
            seamNodes.row(2*i+1)<<surfaceMeshList[meshNum]->V.row(surfaceMeshList[meshNum]->EV(highlightEdges[i],1));
        }
        
        std::string highlightName;
        if (name.empty())
            highlightName = "Edge Highlights " + std::to_string(meshNum);
        else
            highlightName = name;
        psEdgeHighlightList[meshNum] = polyscope::registerCurveNetwork(highlightName, seamNodes, seamEdges);
        psEdgeHighlightList[meshNum]->setColor(glm::vec3());
        psEdgeHighlightList[meshNum]->setRadius(widthRatio*surfaceMeshList[meshNum]->avgEdgeLength, false);
    }
    
    void inline highlight_vertices(const Eigen::VectorXi& highlightVertices,
                                   const std::string name = "",
                                   const int meshNum = 0){
        //TODO: boundary vertices
        std::set<int> highlightFacesList;
        //std::cout<<"surfaceMeshList[meshNum]->vertexValence: "<<surfaceMeshList[meshNum]->vertexValence<<std::endl;
        for (int i=0;i<highlightVertices.size();i++)
            for (int j=0;j<surfaceMeshList[meshNum]->vertexValence(highlightVertices(i));j++)
                highlightFacesList.insert(surfaceMeshList[meshNum]->VF(highlightVertices[i], j));
        
        Eigen::VectorXi highlightFaces(highlightFacesList.size());
        int i = 0;
        for (int val : highlightFacesList)
            highlightFaces(i++) = val;
        
        std::string highlightName;
        if (name.empty())
            highlightName = "Vertex Highlights " + std::to_string(meshNum);
        else
            highlightName = name;
        
        highlight_faces(highlightFaces, highlightName, meshNum);
    }
    
    void inline init_streamlines(const CartesianField& field,
                                 const int fieldNum = 0,
                                 const Eigen::VectorXi& seedLocations=Eigen::VectorXi(),
                                 const double distRatio=3.0)
    {
        if (slData.size()<fieldNum+1){
            slData.resize(fieldNum+1);
            slState.resize(fieldNum+1);
            psStreamlineList.resize(fieldNum+1);
        }
        //assert(fieldList[meshNum]->tb->discTangType()==discTangTypeEnum::FACE_SPACES);
        directional::streamlines_init(field, seedLocations,distRatio,slData[fieldNum], slState[fieldNum]);
    }
    
    void inline advance_streamlines(const double dTimeRatio,
                                    const int fieldNum=0,
                                    const double widthRatio=0.05,
                                    const double colorAttenuationRate = 0.9){
        
        double dTime = dTimeRatio*fieldList[fieldNum]->tb->avgAdjLength;
        directional::streamlines_next(slData[fieldNum], slState[fieldNum],dTime);
        double width = widthRatio*fieldList[fieldNum]->tb->avgAdjLength;
        //TODO: attentuation
        
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
        Eigen::MatrixXd P1(slState[fieldNum].segStart.size(),3), P2(slState[fieldNum].segEnd.size(),3);
        for (int i=0;i<slState[fieldNum].segStart.size();i++){
            P1.row(i)=slState[fieldNum].segStart[i];
            P2.row(i)=slState[fieldNum].segEnd[i];
        }
        Eigen::MatrixXd nodes(P1.rows()+P2.rows(),3);
        nodes<<P1, P2;
        Eigen::MatrixXi edges(P1.rows(),2);
        edges.col(0) = Eigen::VectorXi::LinSpaced(P1.rows(), 0, P1.rows()-1);
        edges.col(1) = Eigen::VectorXi::LinSpaced(P2.rows(), P1.rows(), P1.rows()+P2.rows()-1);
        psStreamlineList[fieldNum] = polyscope::registerCurveNetwork("Streamlines " + std::to_string(fieldNum), nodes, edges)->setRadius(width, false);
    }
    
    void inline set_isolines(const directional::TriMesh& cutMesh,
                             const Eigen::MatrixXd& vertexFunction,
                             const std::string name = "",
                             const int meshNum=0,
                             const int fieldNum=0,
                             const double sizeRatio=0.05)
    {
        
        if (psIsolineList.size()<fieldNum+1)
            psIsolineList.resize(fieldNum+1);
        
        Eigen::MatrixXd isoV, isoN;
        Eigen::MatrixXi isoE, isoOrigE;
        Eigen::VectorXi funcNum;
        
        directional::branched_isolines(cutMesh.V, cutMesh.F, vertexFunction, isoV, isoE, isoOrigE, isoN, funcNum);
        
       // double l = sizeRatio*cutMesh.avgEdgeLength;
        double width = sizeRatio*fieldList[fieldNum]->tb->avgAdjLength;
        
        Eigen::MatrixXd VIso, CIso;
        Eigen::MatrixXi FIso;
        
        Eigen::MatrixXd funcColors = isoline_colors();
        Eigen::MatrixXd CFunc;
        CFunc.resize(funcNum.size(),3);
        for (int i=0;i<funcNum.size();i++)
            CFunc.row(i)=funcColors.row(funcNum(i));
        
        std::string isolineName;
        if (name.empty())
            isolineName = "Isolines m" + std::to_string(meshNum) + "-f" + std::to_string(fieldNum);
        else
            isolineName = name;
        
        psIsolineList[fieldNum] = polyscope::registerCurveNetwork(isolineName, isoV, isoE)->setRadius(width, false);
        psIsolineList[fieldNum]->addEdgeColorQuantity("Function number", CFunc)->setEnabled(true);
        
        //directional::bar_chains(cutMesh, isoV,isoE,isoOrigE, isoN,l,(funcNum.template cast<double>().array()+1.0)*l/1000.0,CFunc, VIso, FIso, CIso);
        //psIsolineList[fieldNum]=polyscope::registerSurfaceMesh("isolines "+std::to_string(fieldNum), VIso, FIso);
        //psIsolineList[fieldNum]->addFaceColorQuantity("branches "+std::to_string(fieldNum), CIso);
    }
    
    void inline set_uv(const Eigen::MatrixXd UV,
                       const std::string name="",
                       const int meshNum=0)
    {
        std::string UVName;
        if (name.empty())
            UVName = "UV " + std::to_string(meshNum);
        else
            UVName = name;
        std::pair<glm::vec3, glm::vec3> gridColors(glm::vec3({1.0,1.0,1.0}), glm::vec3({0.0,0.0,0.0}));
        psSurfaceMeshList[meshNum]->addVertexParameterizationQuantity(UVName, UV, polyscope::ParamCoordsType::UNIT)->setStyle(polyscope::ParamVizStyle::GRID)->setGridColors(gridColors)->setCheckerSize(1.0);
    }
    
    /*void inline set_texture(const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
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
    
    
    void inline toggle_cartesian_field(const bool active, const int fieldNum=0){
        if (fieldNum+1>psGlyphList.size())
            return;  //just ignore the command
        psGlyphSourceList[fieldNum]->setEnabled(active);
        for (int i=0;i<psGlyphList[fieldNum].size();i++)
            psGlyphList[fieldNum][i]->setEnabled(active);
    }
    
    void inline toggle_singularities(const bool active, const int fieldNum=0){
        psSingList[fieldNum]->setEnabled(active);
    }
    
    void inline toggle_edge_highlights(const bool active, const int fieldNum=0){
        psEdgeHighlightList[fieldNum]->setEnabled(active);
    }
    
    void inline toggle_streamlines(const bool active, const int meshNum=0){
        psStreamlineList[meshNum]->setEnabled(active);
    }
    
    
    
    /*void inline toggle_isolines(const bool active, const int meshNum=0){
     data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_faces=active;
     }
     
     void inline toggle_texture(const bool active, const int meshNum=0){
     data_list[NUMBER_OF_SUBMESHES*meshNum].show_texture=active;
     }*/
    
    //disabling the original mesh
    void inline toggle_edge_data(const bool active, const int meshNum=0){
        //TODO: save edge data somewhere
    }
    
    //static functions for default values
    /*//Mesh colors
     static Eigen::RowVector3d inline default_mesh_color(){
     return Eigen::RowVector3d::Constant(1.0);
     }*/
    
    void inline set_glyph_length(double lengthRatio, int fieldNum = 0){
        for (int i=0;i<psGlyphList[fieldNum].size();i++)
            psGlyphList[fieldNum][i]->setVectorLengthScale(lengthRatio*fieldList[fieldNum]->tb->avgAdjLength, false);
    }
    
    void inline toggle_combed_colors(const bool isCombed, const bool signSymmetry=true, int fieldNum=0){
        for (int i=0;i<psGlyphList[fieldNum].size();i++)
            if (!isCombed)
                psGlyphList[fieldNum][i]->setVectorColor(default_glyph_color());
            else
                psGlyphList[fieldNum][i]->setVectorColor(indexed_glyph_colors(i, psGlyphList[fieldNum].size(), signSymmetry));
    }
    
    void inline toggle_field_highlight(const bool isHighlighted, int fieldNum = 0){
        for (int i=0;i<psGlyphList[fieldNum].size();i++)
            if (!isHighlighted)
                psGlyphList[fieldNum][i]->setVectorColor(default_glyph_color());
            else
                psGlyphList[fieldNum][i]->setVectorColor(default_vector_constraint_color());
    }
    
    void inline set_field_color(const glm::vec3 fieldColor, int fieldNum = 0){
        for (int i=0;i<psGlyphList[fieldNum].size();i++)
            psGlyphList[fieldNum][i]->setVectorColor(fieldColor);
    }
    
    static glm::vec3 inline default_face_color(){
        return glm::vec3(243.0/255.0,241.0/255.0,216.0/255.0);
    }
    
    //Color for faces that are selected for editing and constraints
    static glm::vec3 inline highlight_face_color(){
        return glm::vec3(0.7,0.2,0.2);
    }
    
    //Glyph colors
    static glm::vec3 inline default_glyph_color(){
        return glm::vec3(0.05,0.05,1.0);
    }
    
    static glm::vec3 inline default_vector_constraint_color(){
        return glm::vec3(1.0,1.0,0.2);
    }
    
    //a dynamic range of 8, independent of N (so singularity indeices differentiated would be -4/N..4/N
    static glm::vec3 inline default_index_color(const int index){
        if (index<=-4)
            return glm::vec3(0.0,0.0,0.0);
        if (index>=4)
            return glm::vec3(1.0,1.0,1.0);
        switch(index){
            case -3: return glm::vec3(1.0,1.0,0.0); break;
            case -2: return glm::vec3(1.0,0.5,0.0); break;
            case -1: return glm::vec3(1.0,0.0,0.0); break;
            case 0: return default_face_color(); break;
            case 1: return glm::vec3(0.0,1.0,0.0); break;
            case 2: return glm::vec3(0.0,1.0,0.5); break;
            case 3: return glm::vec3(0.0,1.0,1.0); break;
        }
    }
    
    //Glyphs in selected faces
    /*static Eigen::RowVector3d inline selected_face_glyph_color(){
     return Eigen::RowVector3d(223.0/255.0, 210.0/255.0, 16.0/255.0);
     }
     
     //The selected glyph currently edited from a selected face
     static Eigen::RowVector3d inline selected_vector_glyph_color(){
     return Eigen::RowVector3d(0.0,1.0,0.5);
     }*/
    
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
    static glm::vec3 inline indexed_glyph_colors(int vectorNum, int N, bool signSymmetry=true){
        
        assert("Vector number has to be between 0 and 14" && vectorNum>=0 && vectorNum<15);
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
        
        int row = (signSymmetry && (N%2==0) ? vectorNum % (N/2) : vectorNum);
        
        return(glm::vec3(glyphPrincipalColors(row,0),glyphPrincipalColors(row,1), glyphPrincipalColors(row,2)));
        
        
        /*Eigen::MatrixXd fullGlyphColors(field.rows(),field.cols());
         int N = field.cols()/3;
         for (int i=0;i<field.rows();i++)
         for (int j=0;j<N;j++)
         fullGlyphColors.block(i,3*j,1,3)<<
         return fullGlyphColors;*/
    }
    
    //Jet-based singularity colors
    /*static Eigen::MatrixXd inline default_singularity_colors(const int N){
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
    
    Eigen::VectorXi samples_tangent_bundle(const Eigen::MatrixXd& sources,
                                           const Eigen::MatrixXi& adjSpaces,
                                           const int sparsity)
    {
        using namespace Eigen;
        using namespace std;
        VectorXi sampledSpaces;
        if (sparsity!=0){
            //creating adjacency matrix
            vector<Triplet<int>> adjTris;
            for (int i=0;i<adjSpaces.rows();i++)
                if ((adjSpaces(i,0)!=-1)&&(adjSpaces(i,1)!=-1)){
                    adjTris.push_back(Triplet<int>(adjSpaces(i,0), adjSpaces(i,1),1));
                    adjTris.push_back(Triplet<int>(adjSpaces(i,1), adjSpaces(i,0),1));
                }
            
            SparseMatrix<int> adjMat(sources.rows(),sources.rows());
            adjMat.setFromTriplets(adjTris.begin(), adjTris.end());
            SparseMatrix<int> newAdjMat(sources.rows(),sources.rows()),matMult;
            directional::speye(sources.rows(), sources.rows(), matMult);
            for (int i=0;i<sparsity;i++){
                matMult=matMult*adjMat;
                newAdjMat+=matMult;
            }
            
            //cout<<"newAdjMat: "<<newAdjMat<<endl;
            
            adjMat=newAdjMat;
            
            vector<set<int>> ringAdjacencies(sources.rows());
            for (int k=0; k<adjMat.outerSize(); ++k){
                for (SparseMatrix<int>::InnerIterator it(adjMat,k); it; ++it){
                    ringAdjacencies[it.row()].insert(it.col());
                    ringAdjacencies[it.col()].insert(it.row());
                }
            }
            
            VectorXi sampleMask=VectorXi::Zero(sources.rows());
            for (int i=0;i<sources.rows();i++){
                if (sampleMask(i)!=0) //occupied face
                    continue;
                
                sampleMask(i)=2;
                //clearing out all other faces
                for (set<int>::iterator si=ringAdjacencies[i].begin();si!=ringAdjacencies[i].end();si++){
                    if (sampleMask(*si)==0)
                        sampleMask(*si)=1;
                    
                }
            }
            
            vector<int> samplesList;
            for (int i=0;i<sampleMask.size();i++)
                if (sampleMask(i)==2)
                    samplesList.push_back(i);
            
            sampledSpaces = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(samplesList.data(), samplesList.size());
            return sampledSpaces;
        } else return Eigen::VectorXi::LinSpaced(sources.rows(),0,sources.rows()-1);
    }
    
};  //of DirectionalViewer class

}


#endif
