//
// Created by avaxman on 24/10/2023.
//

#ifndef DIRECTIONAL_TUTORIALS_OPENGLVIEWER_H
#define DIRECTIONAL_TUTORIALS_OPENGLVIEWER_H

#include <Eigen/Core>
#include <GLFW/glfw3.h>
#include <directional/Camera.h>

namespace Directional{

    class OpenGLViewer{

    private:
        //GLFW stuff
        GLFWwindow* glfwWindow;
        GLuint vertexBuffer, vertexShader, fragmentShader, glfwProgram;
        void(* callback_key_down) (GLFWwindow *window, int key, int scancode, int action, int mods);
        Camera* camera;

        std::string vertexShaderText, fragmentShaderText;



    public:
        OpenGLViewer(){}
        ~OpenGLViewer(){}

        void init(const std::string& _vertexShaderText, const std::string& _fragmentShaderText){
            vertexShaderText=_vertexShaderText;
            fragmentShaderText= _fragmentShaderText;
            assert(glfwInit() && "GLFW failed to initialize!");
            glfwWindow=NULL;
            callback_key_down=&default_key_down;

            glGenBuffers(1, &vertexBuffer);
            glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);

        }

        void terminate(){
            if (glfwWindow!=NULL)
                glfwDestroyWindow(glfwWindow);
            glfwTerminate();
        }

        bool launch(){
            glfwWindow = glfwCreateWindow(1920, 1280, "Directional Viewer", NULL, NULL);
            if (!glfwWindow)
                return false;

            glfwSetKeyCallback(glfwWindow, callback_key_down);

            glfwMakeContextCurrent(glfwWindow);
            //gladLoadGL(glfwGetProcAddress);
            glfwSwapInterval(1);

            //Default shaders
            vertexShader = glCreateShader(GL_VERTEX_SHADER);
            glShaderSource(vertexShader, 1, &vertexShaderText, NULL);
            glCompileShader(vertexShader);

            fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
            glShaderSource(fragmentShader, 1, &fragmentShaderText, NULL);
            glCompileShader(fragmentShader);

            glfwProgram = glCreateProgram();
            glAttachShader(glfwProgram, vertexShader);
            glAttachShader(glfwProgram, fragmentShader);
            glLinkProgram(glfwProgram);

            while (!glfwWindowShouldClose(glfwWindow))
            {
                float ratio;
                int width, height;

                glfwGetFramebufferSize(glfwWindow, &width, &height);

                glViewport(0, 0, width, height);
                glClear(GL_COLOR_BUFFER_BIT);

                glUseProgram(glfwProgram);
                //glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (const GLfloat*) mvp);
                glDrawArrays(GL_TRIANGLES, 0, 3);

                glfwSwapBuffers(glfwProgram);
                glfwPollEvents();
            }
        }

        void set_mesh(const Eigen::MatrixXd& V,
                      const Eigen::MatrixXi& F,
                      const Eigen::MatrixXi& C){

            glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

        }


        static void default_key_down(GLFWwindow* window, int key, int scancode, int action, int mods)
        {
            if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
                glfwSetWindowShouldClose(window, GLFW_TRUE);
        }



    };
}


#endif //DIRECTIONAL_TUTORIALS_OPENGLVIEWER_H
