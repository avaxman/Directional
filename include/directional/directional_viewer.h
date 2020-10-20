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
#include <igl/opengl/glfw/Viewer.h>
#include <directional/glyph_lines_raw.h>
#include <directional/singularity_spheres.h>
#include <directional/seam_lines.h>




//This file contains the default libdirectional visualization paradigms
namespace directional
  {
  
#define NUMBER_OF_SUBMESHES 5  //triangle mesh, field, singularities, seams, streamlines
  
  class DirectionalViewer: public igl::opengl::glfw::Viewer{
  private:
    std::vector<Eigen::MatrixXd> VMesh;  //vertices of mesh list
    std::vector<Eigen::MatrixXi> FMesh;  //faces of mesh list
    std::vector<Eigen::MatrixXd> CMesh;  //colors of mesh list
    
    //TODO: do I really need to save rawfields?
  public:
    DirectionalViewer(){}
    ~DirectionalViewer(){}
    
    void set_mesh(const Eigen::MatrixXd& V,
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
        
        if (VMesh.size()<meshNum+1){
          VMesh.resize(meshNum+1);
          FMesh.resize(meshNum+1);
          CMesh.resize(meshNum+1);
          
          VMesh[meshNum]=V;
          FMesh[meshNum]=F;
          CMesh[meshNum]=C;
        }
        
      }
    }
    
    void set_mesh_colors(const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                         const int meshNum=0)
    {
      Eigen::MatrixXd meshColors;
      if (C.rows()==0)
        meshColors=default_mesh_color();
      else
        meshColors=C;
      
      data_list[NUMBER_OF_SUBMESHES*meshNum].set_colors(meshColors);
      CMesh[meshNum]=C;
    }
    
    void set_field(const Eigen::MatrixXd& rawField,
                   const Eigen::MatrixXd& C=Eigen::MatrixXd(),
                   const int meshNum=0)
    
    {
      Eigen::MatrixXd fieldColors=C;
      if (C.rows()==0)
        fieldColors=default_glyph_color();
      
      Eigen::MatrixXd VField, CField;
      Eigen::MatrixXi FField;
      directional::glyph_lines_raw(VMesh[meshNum], FMesh[meshNum], rawField, fieldColors, VField, FField, CField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].set_mesh(VField,FField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].set_colors(CField);
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].show_lines=false;
    }
    
    void set_singularities(const int N,
                           const Eigen::VectorXi& singVertices,
                           const Eigen::VectorXi& singIndices,
                           const int meshNum=0)
    {
      Eigen::MatrixXd VSings, CSings;
      Eigen::MatrixXi FSings;
      directional::singularity_spheres(VMesh[meshNum], FMesh[meshNum], N, singVertices, singIndices, VSings, FSings, CSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].set_mesh(VSings,FSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].set_colors(CSings);
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].show_lines=false;
      
    }
    
    void set_seams(const Eigen::MatrixXi& EV,
                   const Eigen::VectorXi& combedMatching,
                   const int meshNum=0)
    {
      Eigen::MatrixXd VSeams, CSeams;
      Eigen::MatrixXi FSeams;
      directional::seam_lines(VMesh[meshNum],FMesh[meshNum],EV,combedMatching, VSeams,FSeams,CSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].clear();
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].set_mesh(VSeams, FSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].set_colors(CSeams);
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].show_lines = false;
    }
    
    void set_streamlines(const Eigen::MatrixXd& P1,
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
    
    void toggle_field(const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+1].show_faces=!data_list[NUMBER_OF_SUBMESHES*meshNum+1].show_faces;
    }
    
    void toggle_singularities(const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+2].show_faces=!data_list[NUMBER_OF_SUBMESHES*meshNum+2].show_faces;
    }
    
    void toggle_seams(const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+3].show_faces=!data_list[NUMBER_OF_SUBMESHES*meshNum+3].show_faces;
    }
    
    void toggle_streamlines(const int meshNum=0){
      data_list[NUMBER_OF_SUBMESHES*meshNum+4].show_faces=!data_list[NUMBER_OF_SUBMESHES*meshNum+4].show_faces;
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
    static Eigen::MatrixXd IGL_INLINE indexed_glyph_colors(const Eigen::MatrixXd& field){
      
      Eigen::Matrix<double, 15,3> glyphPrincipalColors;
      glyphPrincipalColors<<1.0,0.0,0.5,
      1.0,0.5,0.0,
      0.0,1.0,0.5,
      0.0,0.5,1.0,
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
          fullGlyphColors.block(i,3*j,1,3)<<glyphPrincipalColors.row(j);
      
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
    
  };  //of DirectionalViewer class
  
  }


#endif
