//
// Created by clemence on 2020-12-19. AND JACK
//
#include <GLFW/glfw3.h>
#include <fluidsim.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>

int main(int argc, char **argv) {

    /*  Generalized coordinates is a 2n vector of particle positions:
        [x1, y1, x2, y2, ...] , where n is the number of fluid particles defined
        in simulation. */  
    Eigen::VectorXf q = Eigen::VectorXf::Zero(2*NUM_PARTICLES);

    /*  Generalized velocity of fluid particles */
    Eigen::VectorXf qdot = Eigen::VectorXf::Zero(2*NUM_PARTICLES);

    glfwInit();

    // Create window
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    GLFWwindow* window = glfwCreateWindow(800, 800, "OpenGL", nullptr, nullptr);
    glfwMakeContextCurrent(window);

    // initialize main sim variables here

    while(!glfwWindowShouldClose(window))
    {
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
}