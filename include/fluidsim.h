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



void run_one_iteration(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot);

/* Initialize coordinates of fluid simulation */
void init_fluid_sim(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot);

void advect_step(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot, double dt);

//void diffusion_step(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot, double dt);

/* I think we can hardcode the external forces (probably just gravity) inside this function */
void external_forces_step(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot, double dt);

void pressure_project_step(
    Eigen::Ref<Eigen::VectorXf> qdot, 
    Eigen::Ref<const Eigen::VectorXf> q, 
    const double dt,
    const double density);

void bilinear_weights(
    Eigen::Ref<Eigen::Vector4d> weights,
    double x, double y,
    double x1, double x2, double y1, double y2);

double bilinear_interpolate(
    double x, double y,
    double x1, double x2, double y1, double y2,
    double q11, double q12, double q21, double q22);

void linear_weights(
    Eigen::Ref<Eigen::Vector2d> weights,
    double x,
    double x1, double x2);

double linear_interpolate(
    double x,
    double x1, double x2,
    double q1, double q2);

