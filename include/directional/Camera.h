#ifndef DIRECTIONAL_CAMERA_H
#define DIRECTIONAL_CAMERA_H

#include <Eigen/Eigen>

namespace directional {

    class Camera{

        typedef enum ProjectionType{ORTHOGRAPHIC, PERSPECTIVE};
#define TRANSLATION(m,T) m<<0.0,0.0,0.0,T(0),0.0,0.0,0.0,T(1),0.0,0.0,0.0,T(2),0.0,0.0,0.0,1.0;
    private:

        //Camera parameters
        Eigen::Vector3d cameraPos, cameraTarget, cameraUp;
        Eigen::Matrix4d modelViewMat, projectionMat;

        //projection parameters
        double projZoom;
        ProjectionType projType;
        double screenWidth, screenHeight;

        //model parameters
        double sceneRadius;
        Eigen::Vector3d sceneCenter;
    public:

        Camera(){
            reset();
        }

        ~Camera(){}

        void reset() {
            modelViewMat = Eigen::Matrix4d::Identity();
            cameraPos = Eigen::Vector3d(0.0,0.0,-1.0);
            cameraTarget = Eigen::Vector3d();
            cameraUp = Eigen::Vector3d(0.0,1.0,0.0);
            zoom = 0.5;
        }

        void set_model(const Eigen::Vector3d& center, const double radius){
            sceneCenter=center; sceneRadius=radius;
        }

        //reproducing gluLookAt
        void set_camera(const Eigen::Vector3d& pos, const Eigen::Vector3d& target, const Eigen::Vector3d& up){
            cameraPos=pos;
            cameraTarget=target;
            cameraUp=up;
            Eigen::Vector3d f = target - pos;
            f.normalize();
            Eigen::Vector3d s = f.cross(up.normalized());
            s.normalize();
            Eigen::Vector3d u = s.cross(f);
            u.normalize();

            Eigen::Matrix4d rotMat; rotMat<<s, 0.0,u, 0.0,-f, 0.0,Eigen::RowVector3d(), 1.0;
            Eigen::Matrix4d transMat;
            TRANSLATION(transMat,-pos);
            modelViewMat=rotMat*transMat;
        }


        void set_projection(ProjectionType pType, const double zoom, const double fov=55){

            double aspectRatio = (double)screenWidth/screenHeight; // aspect ratio
            projType = pType;
            projZoom=zoom;
            double near = 2.0 * sceneRadius;
            double far = 6.0 * sceneRadius;

            if (projType==ORTHOGRAPHIC){
                double left = -2.0 * sceneRadius * zoom * aspectRatio;
                double right = 2.0 * sceneRadius * zoom * aspectRatio;
                double bottom = -2.0 * sceneRadius * zoom;
                double up = 2.0 * sceneRadius * zoom;

                double tx = -(right+left)/(right-left);
                double ty = -(top+bottom)/(top-bottom);
                double tz = -(far+near)/(far-near);

                projectionMat<<2/(right-left),0.0,0.0,tx,
                0.0,2/(top-bottom),0.0,ty,
                0.0,0.0,-2/(far-near), tz,
                0.0,0.0,0.0,1.0;
            }

            if (projType = PERSPECTIVE){
                double f = cot(((fov*projZoom/180.0)*M_PI/2.0);

                projectionMat<<f/aspectRatio,0.0,0.0,0.0,
                0.0,f,0.0,0.0,
                0.0,0.0,(far+near)/(near-far), 2*far*near/(near-far),
                0.0,0.0,-1.0,0.0;
            }
        }
    };
}

#endif
