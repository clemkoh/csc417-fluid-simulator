cmake_minimum_required(VERSION 3.10)
project(csc417-final-project)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#setup libigl properly

option(LIBIGL_WITH_OPENGL               "Use OpenGL"            ON)
option(LIBIGL_WITH_OPENGL_GLFW          "Use GLFW"              ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI    "Use ImGui"             ON)

cmake_policy(SET CMP0072 NEW)

find_package(LIBIGL REQUIRED QUIET)
find_package(glfw3 REQUIRED)

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${LIBIGL_INCLUDE_DIR})

file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)


add_executable(${PROJECT_NAME} main.cpp  ${H_FILES} ${CPP_FILES})
target_link_libraries(${PROJECT_NAME} glfw igl::core igl::opengl)