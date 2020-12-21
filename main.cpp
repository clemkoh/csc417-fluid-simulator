//
// Created by clemence on 2020-12-19. AND JACK
//
#define GL_GLEXT_PROTOTYPES
#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <fluidsim.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stdio.h>
#include <iostream>

#define RESCALE (1.0/40.0)

int main(int argc, char **argv) {

    /*  Generalized coordinates is a 2n vector of particle positions:
        [x1, y1, x2, y2, ...] , where n is the number of fluid particles defined
        in simulation. */  
    Eigen::VectorXf q = Eigen::VectorXf::Zero(2*NUM_PARTICLES);

    /*  Generalized velocity of fluid particles */
    Eigen::VectorXf qdot = Eigen::VectorXf::Zero(2*NUM_PARTICLES);

    init_fluid_sim(q, qdot);

    run_one_iteration(q, qdot);
    //std::cout << q << "\n";
    //std::cout << qdot << "\n";

    igl::opengl::glfw::Viewer viewer;
    viewer.core().is_animating = true;
    viewer.core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);

    Eigen::MatrixXd points(NUM_PARTICLES,3);
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & )->bool
    {
        run_one_iteration(q, qdot);
        for (int i = 0; i < points.rows(); i++) {
            points(i, 0) = q(2 * i) * RESCALE;
            points(i, 1) = q(2 * i + 1) * RESCALE;
        }
        points = points.eval();
        viewer.data().set_points(points, Eigen::RowVector3d(1,1,1));
        return false;
    };
    viewer.launch();
}