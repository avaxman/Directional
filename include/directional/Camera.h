#ifndef DIRECTIONAL_CAMERA_H
#define DIRECTIONAL_CAMERA_H

#include <Eigen/Eigen>




namespace directional {

#define ORTHOGRAPHIC 1
#define PERSPECTIVE 2

    class Camera{

    private:


        //Camera parameters
        Eigen::Vector3d cameraPos, cameraTarget, cameraUp;
        Eigen::Matrix4d model, view, projection;
        double zoom;


        //model parameters
        double sceneRadius;
        Eigen::Vector3d sceneCenter;
    public:

        Camera(){
            reset();
        }

        ~Camera(){}

        void reset() {
            model = view = projection = Eigen::Matrix4d::Identity();
            pos = Eigen::Vector3d(0.0,0.0,-1.0);
            target = Eigen::Vector3d();
            up = Eigen::Vector3d(0.0,1.0,0.0);
            zoom = 0.5;



        }

        void set_model(const Eigen::Vector3d& center, const double radius){
            sceneCenter=center; sceneRadius=radius;
        }

        void set_camera(const Eigen::Vector3d& pos, const Eigen::Vector3d& target, const Eigen::Vector3d& up, const double zoom){

        }
    };
}

#endif
