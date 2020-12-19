#include <Eigen/Dense>
#include <Eigen/Sparse>

// MAIN fluid variables
#define DENSITY 0.0
#define VISCOSITY 0.0
#define NUM_PARTICLES 8
#define DT 1

// simulation area size
#define LENGTH 100
#define HEIGHT 100

#define GRID_DX 10
#define GRID_DY 10





/* Initialize coordinates of fluid simulation */
void init_fluid_sim(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot);

void advect_step(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot, double dt);

void diffusion_step(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot, double dt);

/* I think we can hardcode the external forces (probably just gravity) inside this function */
void external_forces_step(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot, double dt);

void pressure_project_step(Eigen::Ref<Eigen::VectorXf> density, 
                           Eigen::Ref<const Eigen::VectorXf> q, 
                           Eigen::Ref<const Eigen::VectorXf> qdot, 
                           Eigen::Ref<Eigen::MatrixXf> mac_x, /* staggered grid part of x coords */
                           Eigen::Ref<Eigen::MatrixXf> mac_y, /* staggered grid part of y coords */
                           Eigen::Ref<Eigen::MatrixXf> mac_p, /* staggered grid part of pressure */                           
                           const double dt);