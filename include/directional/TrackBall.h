#ifndef DIRECTIONAL_TRACKBALL_H
#define DIRECTIONAL_TRACKBALL_H

#include <Eigen/Eigen>


namespace directional {

    class TrackBall{
    private:

        Eigen::Quaterniond currRot, lastRot;
        Eigen::Quaterniond startPos, currPos;
    public:

        TrackBall(){
            reset();
        }

        ~TrackBall(){}

        void reset(){lastRot=currRot=Eigen::Quaterniond(1.0,0.0,0.0,0.0);}
        void start(const double px, const double py){
            startPos=Eigen::Quaterniond(0.0,px,py, (px*px+py*py<1 ? sqrt(1-px*px+py*py) : 0.0));
            startPos.normalize();

        }
        void update(const double px, const double py, Eigen::Matrix4d& fullRotMat){
            currPos = Eigen::Quaterniond(0.0,px,py, (px*px+py*py<1 ? sqrt(1-px*px+py*py) : 0.0));
            currPos.normalize();
            currRot =currPos*startPos.inverse();
            fullRotMat<<(currRot*lastRot).normalized().toRotationMatrix(), Eigen::Vector3d(), Eigen::RowVector4d;
        }

        void stop(){
            lastRot = currRot*lastRot;
            currRot = Eigen::Quaterniond(1.0,0.0,0.0,0.0);
        }
    };
}

#endif
