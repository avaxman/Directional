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
#include <directional/glyph_lines_raw.h>
#include <directional/singularity_spheres.h>
#include <directional/seam_lines.h>
#include <directional/edge_diamond_mesh.h>




//This file contains the default libdirectional visualization paradigms
namespace directional
  {
  
#define NUMBER_OF_SUBMESHES 6  //triangle mesh, field, singularities, seams, streamlines, edge-diamond mesh
  
  class DirectionalViewer: public igl::opengl::glfw::Viewer{
  private:
    std::vector<Eigen::MatrixXd> VList;  //vertices of mesh list
    std::vector<Eigen::MatrixXi> FList;  //faces of mesh list
    
    std::vector<Eigen::MatrixXd> edgeVList;  //edge-diamond vertices list
    std::vector<Eigen::MatrixXi> edgeFList;  //edge-diamond faces list
    std::vector<Eigen::VectorXi> edgeFEList;  //edge-diamond faces->original mesh edges list
    
  public:
    DirectionalViewer(){}
    ~DirectionalViewer(){}
    
    void IGL_INLINE set_mesh(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                             const int meshNum=0)
    {
      Eigen::MatrixXd meshColors;
      if (C.rows()==0)
        meshColors=default_mesh_color();
      else
        meshColors=C;
      
      if (NUMBER_OF_SUBMESHES*(meshNum+1)>data_list.size()){  //allocating until there
        int currDLSize=data_list.size();
        for (int i=currDLSize;i<NUMBER_OF_SUBMESHES*(meshNum+1);i++)
          append_mesh();
        
        selected_data_index=NUMBER_OF_SUBMESHES*meshNum;  //the last triangle mesh
        
        data_list[NUMBER_OF_SUBMESHES*meshNum].clear();
        data_list[NUMBER_OF_SUBMESHES*meshNum].set_mesh(V,F);
        data_list[NUMBER_OF_SUBMESHES*meshNum].set_colors(meshColors);
        
        if (VList.size()<meshNum+1){
          VList.resize(meshNum+1);
          FList.resize(meshNum+1);
          edgeVList.resize(meshNum+1);
          edgeFList.resize(meshNum+1);
          edgeFEList.resize(meshNum+1);
          //CList.resize(meshNum+1);
          
          VList[meshNum]=V;
          FList[meshNum]=F;
          //CList[meshNum]=C;
        }
        
      }
    }
    
    void IGL_INLINE set_mesh_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                                    const int meshNum=0)
    {
      Eigen::MatrixXd meshColors;
      if (C.rows()==0)
        meshColors=default_mesh_color();
      else
        meshColors=C;
      
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_colors(meshColors);
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_faces=true;
      if (edgeVList[meshNum].size()!=0){
        data_list[NUMBER_OF_SUBMESHES*meshNum+5].show_faces=false;
        data_list[NUMBER_OF_SUBMESHES*meshNum+5].show_lines=false;
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
                                  const Eigen::MatrixXi& EV,
                                  const Eigen::MatrixXi& FE,
                                  const Eigen::MatrixXi& EF,
                                  const int meshNum=0)
    {
      
      if (edgeVList[meshNum].size()==0){  //allocate
        edge_diamond_mesh(VList[meshNum],FList[meshNum],EV,EF,edgeVList[meshNum],edgeFList[meshNum],edgeFEList[meshNum]);
        data_list[NUMBER_OF_SUBMESHES*meshNum+5].clear();
        data_list[NUMBER_OF_SUBMESHES*meshNum+5].set_mesh(edgeVList[meshNum],edgeFList[meshNum]);
      }
      
      Eigen::VectorXd edgeFData(edgeFList[meshNum].rows());
      for (int i=0;i<edgeFList[meshNum].rows();i++)
        edgeFData(i)=edgeData(edgeFEList[meshNum](i));
      
      Eigen::MatrixXd C;
      igl::parula(edgeFData, minRange,maxRange, C);
      data_list[NUMBER_OF_SUBMESHES*meshNum+5].set_colors(C);
      
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_faces=false;
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_lines=false;
      data_list[NUMBER_OF_SUBMESHES*meshNum+5].show_faces=true;
      data_list[NUMBER_OF_SUBMESHES*meshNum+5].show_lines=false;
      
      selected_data_index=NUMBER_OF_SUBMESHES*meshNum+5;
    }
    
    void IGL_INLINE set_selected_faces(const Eigen::VectorXi& selectedFaces, const int meshNum=0){
      Eigen::MatrixXd CMesh=directional::DirectionalViewer::default_mesh_color().replicate(FList[meshNum].rows(),1);
      for (int i=0;i<selectedFaces.size();i++)
        CMesh.row(selectedFaces(i))=selected_face_color();
      set_mesh_colors(CMesh,meshNum);
    }
    
    void IGL_INLINE set_selected_vector(const Eigen::MatrixXd& rawField, const int selectedFace, const int selectedVector, const int meshNum=0)
    {
      int N = rawField.cols()/3;
      Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(FList[meshNum].rows(),N);
      glyphColors.row(selectedFace)=directional::DirectionalViewer::selected_face_glyph_color().replicate(1,N);
      glyphColors.block(selectedFace,3*selectedVector,1,3)=directional::DirectionalViewer::selected_vector_glyph_color();
      
      set_field(rawField, glyphColors);
    }
    
    
    
    
    void IGL_INLINE set_field(const Eigen::MatrixXd& rawField,
                              const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                              const int meshNum=0)
    
    {
      Eigen::MatrixXd fieldColors=C;
      if (C.rows()==0)
        fieldColors=default_glyph_color();
      
      Eigen::MatrixXd VField, CField;
      Eigen::MatrixXi FField;
      directional::glyph_lines_raw(VList[meshNum], FList[meshNum], rawField, fieldColors, VField, FField, CField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].set_mesh(VField,FField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].set_colors(CField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].show_lines=false;
    }
    
    void IGL_INLINE set_singularities(const int N,
                                      const Eigen::VectorXi& singVertices,
                                      const Eigen::VectorXi& singIndices,
                                      const int meshNum=0)
    {
      Eigen::MatrixXd VSings, CSings;
      Eigen::MatrixXi FSings;
      directional::singularity_spheres(VList[meshNum], FList[meshNum], N, singVertices, singIndices, default_singularity_colors(N), VSings, FSings, CSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].set_mesh(VSings,FSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].set_colors(CSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].show_lines=false;
      
    }
    
    void IGL_INLINE set_seams(const Eigen::MatrixXi& EV,
                              const Eigen::VectorXi& combedMatching,
                              const int meshNum=0)
    {
      Eigen::MatrixXd VSeams, CSeams;
      Eigen::MatrixXi FSeams;
      directional::seam_lines(VList[meshNum], FList[meshNum],EV,combedMatching, default_seam_color(), VSeams,FSeams,CSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].set_mesh(VSeams, FSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].set_colors(CSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].show_lines = false;
    }
    
    void IGL_INLINE set_streamlines(const Eigen::MatrixXd& P1,
                                    const Eigen::MatrixXd& P2,
                                    const Eigen::MatrixXd& C,
                                    const int meshNum=0,
                                    const double width=0.0005)
    {
      Eigen::MatrixXd VStream, CStream;
      Eigen::MatrixXi FStream;
      directional::line_cylinders(P1,P2, width, C, 4, VStream, FStream, CStream);
      
      data_list[NUMBER_OF_SUBMESHES*meshNum+4].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+4].set_mesh(VStream, FStream);
      data_list[NUMBER_OF_SUBMESHES*meshNum+4].set_colors(CStream);
      data_list[NUMBER_OF_SUBMESHES*meshNum+4].show_lines = false;
    }
    
    void IGL_INLINE set_uv(const Eigen::MatrixXd UV,
                           const int meshNum=0)
    {
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_uv(UV);
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_texture=true;
    }
    
    void IGL_INLINE set_texture(const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& textureMat, const int meshNum=0)
    {
      int textureSize=textureMat.rows()/3;
      
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_texture(textureMat.block(0,0,textureSize,textureMat.cols()),
                                                         textureMat.block(textureSize,0,textureSize,textureMat.cols()),
                                                         textureMat.block(2*textureSize,0,textureSize,textureMat.cols()));
    }
    
    void IGL_INLINE set_active(const bool active, const int meshNum=0){
      for (int i=NUMBER_OF_SUBMESHES*meshNum;i<NUMBER_OF_SUBMESHES*meshNum+NUMBER_OF_SUBMESHES;i++)
        data_list[i].show_faces=active;
    }
    
    void IGL_INLINE toggle_mesh_edges(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_lines=active;
    }
    
    void IGL_INLINE toggle_field(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].show_faces=active;
    }
    
    void IGL_INLINE toggle_singularities(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].show_faces=active;
    }
    
    void IGL_INLINE toggle_seams(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].show_faces=active;
    }
    
    void IGL_INLINE toggle_streamlines(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+4].show_faces=active;
    }
    
    void IGL_INLINE toggle_uv(const bool active, const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum].show_texture=active;
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
