// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_VIEWER_H
#define DIRECTIONAL_VIEWER_H

#include <Eigen/Core>
#include <igl/jet.h>
#include <igl/parula.h>
#include <igl/opengl/glfw/Viewer.h>
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
#include <directional/FaceField.h>
#include <igl/edge_topology.h>


//This file contains the default libdirectional visualization paradigms
namespace directional
  {
  
#define NUMBER_OF_SUBMESHES 7
#define FIELD_MESH 1
#define SING_MESH 2
#define SEAMS_MESH 3
#define STREAMLINE_MESH 4
#define EDGE_DIAMOND_MESH 5
#define ISOLINES_MESH 6
  
  class DirectionalViewer: public igl::opengl::glfw::Viewer{
  private:
    std::vector<const TriMesh*> meshList;  //meshes that are being viewed

    std::vector<const FaceField*> fieldList;
    std::vector<Eigen::MatrixXd> fieldColors;
    std::vector<directional::StreamlineData> slData;
    std::vector<directional::StreamlineState> slState;
    
    std::vector<Eigen::MatrixXd> edgeVList;  //edge-diamond vertices list
    std::vector<Eigen::MatrixXi> edgeFList;  //edge-diamond faces list
    std::vector<Eigen::VectorXi> edgeFEList;  //edge-diamond faces->original mesh edges list
    
    std::vector<Eigen::MatrixXd> fieldVList;
    std::vector<Eigen::MatrixXi> fieldFList;
    
  public:
    DirectionalViewer(){}
    ~DirectionalViewer(){}
    
    void IGL_INLINE set_mesh(const TriMesh& mesh,
                             const int meshNum=0)
    {
      Eigen::MatrixXd meshColors;
      meshColors=default_mesh_color();

      if (NUMBER_OF_SUBMESHES*(meshNum+1)>data_list.size()){  //allocating until there
        int currDLSize=data_list.size();
        for (int i=currDLSize;i<NUMBER_OF_SUBMESHES*(meshNum+1);i++)
          append_mesh();
      }
        
      selected_data_index=NUMBER_OF_SUBMESHES*meshNum;  //the last triangle mesh
      data_list[NUMBER_OF_SUBMESHES*meshNum].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_mesh(mesh.V,mesh.F);
      if ((mesh.V.rows()==mesh.F.rows())||(meshColors.rows()!=mesh.V.rows()))  //assume face_based was meant
        data_list[NUMBER_OF_SUBMESHES*meshNum].set_face_based(true);
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_colors(meshColors);
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_lines=false;
      
      if (meshList.size()<meshNum+1){
        meshList.resize(meshNum+1);
        fieldList.resize(meshNum+1);
        edgeVList.resize(meshNum+1);
        edgeFList.resize(meshNum+1);
        edgeFEList.resize(meshNum+1);
      }
      meshList[meshNum]=&mesh;
    }
    
    void IGL_INLINE set_mesh_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
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
    }
    
    void IGL_INLINE set_vertex_data(const Eigen::VectorXd& vertexData,
                                    const double minRange,
                                    const double maxRange,
                                    const int meshNum=0)
    {
      Eigen::MatrixXd C;
      igl::parula(vertexData, minRange,maxRange, C);
      set_mesh_colors(C, meshNum);
    }
    
    void IGL_INLINE set_face_data(const Eigen::VectorXd& faceData,
                                  const double minRange,
                                  const double maxRange,
                                  const int meshNum=0)
    {
      Eigen::MatrixXd C;
      igl::parula(faceData, minRange,maxRange, C);
      set_mesh_colors(C, meshNum);
    }
    
    void IGL_INLINE set_edge_data(const Eigen::VectorXd& edgeData,
                                  const double minRange,
                                  const double maxRange,
                                  const int meshNum=0)
    {
    
      if (edgeVList[meshNum].size()==0){  //allocate
        edge_diamond_mesh(meshList[meshNum]->V,meshList[meshNum]->F,meshList[meshNum]->EV,meshList[meshNum]->EF,edgeVList[meshNum],edgeFList[meshNum],edgeFEList[meshNum]);
        data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].clear();
        data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].set_mesh(edgeVList[meshNum],edgeFList[meshNum]);
      }
      
      Eigen::VectorXd edgeFData(edgeFList[meshNum].rows());
      for (int i=0;i<edgeFList[meshNum].rows();i++)
        edgeFData(i)=edgeData(edgeFEList[meshNum](i));
      
      Eigen::MatrixXd C;
      igl::parula(edgeFData, minRange,maxRange, C);
      data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].set_colors(C);
      
      data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_faces=false;
      data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_lines=false;
      
      selected_data_index=NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH;
    }
    
    void IGL_INLINE set_selected_faces(const Eigen::VectorXi& selectedFaces, const int meshNum=0){
      Eigen::MatrixXd CMesh=directional::DirectionalViewer::default_mesh_color().replicate(meshList[meshNum]->F.rows(),1);
      for (int i=0;i<selectedFaces.size();i++)
        CMesh.row(selectedFaces(i))=selected_face_color();
      set_mesh_colors(CMesh,meshNum);
      
      //coloring field
      Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->F.rows(),fieldList[meshNum]->N);
      for (int i=0;i<selectedFaces.rows();i++)
        glyphColors.row(selectedFaces(i))=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);
      
      set_field_colors(glyphColors,meshNum);
    }
    
    void IGL_INLINE set_selected_vector(const int selectedFace, const int selectedVector, const int meshNum=0)
    {
      Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(meshList[meshNum]->F.rows(),fieldList[meshNum]->N);
      glyphColors.row(selectedFace)=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,fieldList[meshNum]->N);
      glyphColors.block(selectedFace,3*selectedVector,1,3)=directional::DirectionalViewer::selected_vector_glyph_color();
      set_field_colors(glyphColors, meshNum);
    }
    
    void IGL_INLINE set_field(const FaceField& _field,
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
      directional::glyph_lines_mesh(meshList[meshNum]->V, meshList[meshNum]->F, meshList[meshNum]->EF, fieldList[meshNum]->extField, fieldColors[meshNum], VField, FField, CField,sizeRatio, sparsity, offsetRatio);
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_mesh(VField,FField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_colors(CField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].show_lines=false;
      
      set_singularities(fieldList[meshNum]->singVertices,
                        fieldList[meshNum]->singIndices,
                        meshNum);
    }
    
    void IGL_INLINE set_field_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
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
      directional::glyph_lines_mesh(meshList[meshNum]->V, meshList[meshNum]->F, meshList[meshNum]->EF, fieldList[meshNum]->extField, fieldColors[meshNum], VField, FField, CField,sizeRatio, sparsity);
      
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_mesh(VField,FField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].set_colors(CField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].show_lines=false;
    }
    
    
    void IGL_INLINE set_singularities(const Eigen::VectorXi& singVertices,
                                      const Eigen::VectorXi& singIndices,
                                      const int meshNum=0,
                                      const double radiusRatio=1.25)
    {
      Eigen::MatrixXd VSings, CSings;
      Eigen::MatrixXi FSings;
      directional::singularity_spheres(meshList[meshNum]->V, meshList[meshNum]->F, fieldList[meshNum]->N, singVertices, singIndices, default_singularity_colors(fieldList[meshNum]->N), VSings, FSings, CSings, radiusRatio);
      data_list[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].set_mesh(VSings,FSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].set_colors(CSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].show_lines=false;
      
    }
    
    void IGL_INLINE set_seams(const Eigen::VectorXi& combedMatching,
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

      data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].set_mesh(VSeams, FSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].set_colors(CSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].show_lines = false;
    }
    
    
    
    void IGL_INLINE init_streamlines(const int meshNum=0,
                                     const Eigen::VectorXi& seedLocations=Eigen::VectorXi(),
                                     const int sparsity=3)
    {
      if (slData.size()<meshNum+1){
        slData.resize(meshNum+1);
        slState.resize(meshNum+1);
      }
      Eigen::MatrixXd stam;
      directional::streamlines_init(meshList[meshNum]->V, meshList[meshNum]->F, fieldList[meshNum]->extField, seedLocations,sparsity,slData[meshNum], slState[meshNum]);
      
    }
    
    void IGL_INLINE advance_streamlines(const int meshNum=0,
                                        const double widthRatio=0.05,
                                        const double colorAttenuation = 0.9){
      
      directional::streamlines_next(meshList[meshNum]->V, meshList[meshNum]->F, slData[meshNum], slState[meshNum]);
      double width = widthRatio*igl::avg_edge_length(meshList[meshNum]->V, meshList[meshNum]->F);
      
      //generating colors according to original elements and their time signature
      Eigen::MatrixXd slColors(slState[meshNum].P1.rows(),3);
      
      for (int i=0;i<slState[meshNum].origFace.size();i++){
        if (fieldColors[meshNum].rows()==1)
          slColors.row(i)=fieldColors[meshNum];
        else{
          double blendFactor = pow(colorAttenuation,(double)slState[meshNum].timeSignature(i)-1.0);
          //std::cout<<"slState[meshNum].origVector(i): "<<slState[meshNum].origVector(i)<<std::endl;
          slColors.row(i)=fieldColors[meshNum].block(slState[meshNum].origFace(i), 3*slState[meshNum].origVector(i), 1,3);
          slColors.row(i).array()=slColors.row(i).array()*blendFactor+default_mesh_color().array()*(1.0-blendFactor);
        }
      }
      
          
      Eigen::MatrixXd VStream, CStream;
      Eigen::MatrixXi FStream;
      directional::line_cylinders(slState[meshNum].P1,slState[meshNum].P2, width, slColors, 4, VStream, FStream, CStream);
      data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].set_mesh(VStream, FStream);
      data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].set_colors(CStream);
      data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].show_lines = false;
      
      
      
    }
    
    void IGL_INLINE set_isolines(const Eigen::MatrixXd& cutV,
                                 const Eigen::MatrixXi& cutF,
                                 const Eigen::MatrixXd& vertexFunction,
                                 const int meshNum=0,
                                 const double sizeRatio=0.1)
    {
      
      
      Eigen::MatrixXd isoV, isoN;
      Eigen::MatrixXi isoE, isoOrigE;
      Eigen::VectorXi funcNum;
      
      directional::branched_isolines(cutV, cutF, vertexFunction, isoV, isoE, isoOrigE, isoN, funcNum);
      
      double l = sizeRatio*igl::avg_edge_length(cutV, cutF);
      
      Eigen::MatrixXd VIso, CIso;
      Eigen::MatrixXi FIso;
      
      Eigen::MatrixXd funcColors = isoline_colors();
      Eigen::MatrixXd CFunc;
      CFunc.resize(funcNum.size(),3);
      for (int i=0;i<funcNum.size();i++)
        CFunc.row(i)=funcColors.row(funcNum(i));
      
      directional::bar_chains(cutV, cutF, isoV,isoE,isoOrigE, isoN,l,(funcNum.template cast<double>().array()+1.0)*l/1000.0,CFunc, VIso, FIso, CIso);

      data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].set_mesh(VIso, FIso);
      data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].set_colors(CIso);
      data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_lines = false;
    }
    
    void IGL_INLINE set_uv(const Eigen::MatrixXd UV,
                           const int meshNum=0)
    {
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_uv(UV);
    }
    
    void IGL_INLINE set_texture(const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
                                const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
                                const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
                                const int meshNum=0)
    {
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_texture(R,G,B);
    }
    
    void IGL_INLINE set_active(const bool active, const int meshNum=0){
      for (int i=NUMBER_OF_SUBMESHES*meshNum;i<NUMBER_OF_SUBMESHES*meshNum+NUMBER_OF_SUBMESHES;i++)
        data_list[i].show_faces=active;
    }
    
    void IGL_INLINE toggle_mesh(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_faces=active;
    }
    
    void IGL_INLINE toggle_mesh_edges(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_lines=active;
    }
    
    void IGL_INLINE toggle_field(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+FIELD_MESH].show_faces=active;
    }
    
    void IGL_INLINE toggle_singularities(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+SING_MESH].show_faces=active;
    }
    
    void IGL_INLINE toggle_seams(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+SEAMS_MESH].show_faces=active;
    }
    
    void IGL_INLINE toggle_streamlines(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+STREAMLINE_MESH].show_faces=active;
    }
    
    void IGL_INLINE toggle_isolines(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+ISOLINES_MESH].show_faces=active;
    }
    
    void IGL_INLINE toggle_texture(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_texture=active;
    }
    
    //disabling the original mesh
    void IGL_INLINE toggle_edge_data(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+EDGE_DIAMOND_MESH].show_faces=active;
      //data_list[NUMBER_OF_SUBMESHES*meshNum].show_faces=!active;
    }
    
    //static functions for default values
    //Mesh colors
    static Eigen::RowVector3d IGL_INLINE default_mesh_color(){
      return Eigen::RowVector3d::Constant(1.0);
    }
    
    //Color for faces that are selected for editing and constraints
    static Eigen::RowVector3d IGL_INLINE selected_face_color(){
      return Eigen::RowVector3d(0.7,0.2,0.2);
    }
    
    //Glyph colors
    static Eigen::RowVector3d IGL_INLINE default_glyph_color(){
      return Eigen::RowVector3d(0.0,0.2,1.0);
    }
    
    //Glyphs in selected faces
    static Eigen::RowVector3d IGL_INLINE selected_face_glyph_color(){
      return Eigen::RowVector3d(223.0/255.0, 210.0/255.0, 16.0/255.0);
    }
    
    //The selected glyph currently edited from a selected face
    static Eigen::RowVector3d IGL_INLINE selected_vector_glyph_color(){
      return Eigen::RowVector3d(0.0,1.0,0.5);
    }
    
    //Colors by indices in each directional object.
    static Eigen::MatrixXd IGL_INLINE isoline_colors(){
      
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
    static Eigen::MatrixXd IGL_INLINE indexed_glyph_colors(const Eigen::MatrixXd& field, bool signSymmetry=true){
      
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
    static Eigen::MatrixXd IGL_INLINE default_singularity_colors(const int N){
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
    static Eigen::RowVector3d IGL_INLINE default_seam_color(){
      return Eigen::RowVector3d(0.0,0.0,0.0);
    }
    
    static Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> IGL_INLINE default_texture(){
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
    
  };  //of DirectionalViewer class
  
  }


#endif
