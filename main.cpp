//
// Created by clemence on 2020-12-19.
//
#include <GLFW/glfw3.h>

int main(int argc, char **argv) {
    glfwInit();

    // Create window
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);
    GLFWwindow* window = glfwCreateWindow(800, 800, "OpenGL", nullptr, nullptr);
    glfwMakeContextCurrent(window);

    while(!glfwWindowShouldClose(window))
    {
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwTerminate();
}