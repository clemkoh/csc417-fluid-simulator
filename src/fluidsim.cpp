#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fluidsim.h>


/* TODO: set this to non zero */
void init_fluid_sim(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot) {
    q = Eigen::VectorXf::Zero(2*NUM_PARTICLES);
    qdot = Eigen::VectorXf::Zero(2*NUM_PARTICLES);    
}



void pressure_project_step(Eigen::Ref<Eigen::VectorXf> density, 
                           Eigen::Ref<const Eigen::VectorXf> q, 
                           Eigen::Ref<const Eigen::VectorXf> qdot, 
                           Eigen::Ref<Eigen::MatrixXf> mac_x, /* staggered grid part of x coords */
                           Eigen::Ref<Eigen::MatrixXf> mac_y, /* staggered grid part of y coords */
                           Eigen::Ref<Eigen::MatrixXf> mac_p, /* staggered grid part of pressure */                           
                           const double dt) {
    
    
    /* collect particles to MAC grid, maybe pass this as an argument */
    Eigen::VectorXf cells[LENGTH/GRID_DX][HEIGHT/GRID_DY]; /* scratch space for collecting particles */



    /* calculate the new pressure gradient */

    /* convert back to particle form with new velocities */
}