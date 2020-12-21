#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <algorithm>
#include <utility>
#include <fluidsim.h>


/* TODO: set this to non zero */
void init_fluid_sim(Eigen::Ref<Eigen::VectorXf> q, Eigen::Ref<Eigen::VectorXf> qdot) {
    q = Eigen::VectorXf::Zero(2*NUM_PARTICLES);
    qdot = Eigen::VectorXf::Zero(2*NUM_PARTICLES);    
}



void pressure_project_step(
    Eigen::Ref<Eigen::VectorXf> qdot, 
    Eigen::Ref<const Eigen::VectorXf> q, 
    const double dt,
    const double density) 
{    
    int num_cells_x = LENGTH/GRID_DX;
    int num_cells_y = HEIGHT/GRID_DY;

    /* collect particles to MAC grid*/
    std::vector<int> cells[2*num_cells_y][2*num_cells_x]; 
    
    Eigen::MatrixXd mac_x = Eigen::MatrixXd::Zero(num_cells_y, num_cells_x+1);
    Eigen::MatrixXd mac_y = Eigen::MatrixXd::Zero(num_cells_y+1, num_cells_x);    

    for (int i = 0; i < NUM_PARTICLES; i++) {
        double qx = q(2*i);
        double qy = q(2*i + 1);
        double qdotx = qdot(2*i);
        double qdoty = qdot(2*i +1);
        cells[((int)(2*qy))/GRID_DY][((int)(2*qx))/GRID_DX].push_back(i);
    }

    /* no idea how the memory layout is for this so we might not be accessing it in the best way */
    for (int yi = 0; yi < 2*num_cells_y; yi++) {
        for (int xi = 0; xi < 2*num_cells_x; xi++) {
            for (int p=0; p < cells[yi][xi].size(); p++) {

                int index = cells[yi][xi][p];
                double qdotx = qdot(2*index);
                double qdoty = qdot(2*index+1);
                double qx = q(2*index);
                double qy = q(2*index+1);

                if (yi == 0) {
                    Eigen::Vector2d w;
                    linear_weights(w, qx, (xi/2)*GRID_DX, (xi/2 +1)*GRID_DX);
                    mac_x(0, xi/2) += w(0)*qdotx;
                    mac_x(0, xi/2 +1) += w(1)*qdotx;
                } else if (yi == 2*num_cells_y -1) {
                    Eigen::Vector2d w;
                    linear_weights(w, qx, (xi/2)*GRID_DX, (xi/2 +1)*GRID_DX);
                    mac_x((yi-1)/2, xi/2) += w(0)*qdotx;
                    mac_x((yi-1)/2, xi/2 +1) += w(1)*qdotx;
                } else { /* i think this is right? */
                    Eigen::Vector4d w;
                    bilinear_weights(w, qx, qy, (xi/2)*GRID_DX, (xi/2 +1)*GRID_DX, ((yi-1)/2)*GRID_DY, ((yi-1)/2 +1)*GRID_DY);
                    mac_x((yi-1)/2, xi/2) += w(0)*qdotx;
                    mac_x((yi-1)/2 + 1, xi/2) += w(1)*qdotx;
                    mac_x((yi-1)/2, xi/2 + 1) += w(2)*qdotx;                    
                    mac_x((yi-1)/2 + 1, xi/2 +1) += w(3)*qdotx;
                }

                if (xi == 0) {
                    Eigen::Vector2d w;
                    linear_weights(w, qy, (yi/2)*GRID_DY, (yi/2 +1)*GRID_DY);
                    mac_y(yi/2, 0) += w(0)*qdoty;
                    mac_y(yi/2 +1, 0) += w(1)*qdoty;
                } else if (xi == 2*num_cells_x -1) {
                    Eigen::Vector2d w;
                    linear_weights(w, qy, (yi/2)*GRID_DY, (yi/2 +1)*GRID_DY);
                    mac_x(yi/2, (xi-1)/2) += w(0)*qdoty;
                    mac_x(yi/2 +1, (xi-1)/2) += w(1)*qdoty;
                } else {
                    Eigen::Vector4d w;
                    bilinear_weights(w, qx, qy, ((xi-1)/2)*GRID_DX, ((xi-1)/2 +1)*GRID_DX, (yi/2)*GRID_DY, (yi/2 +1)*GRID_DY);
                    mac_y(yi/2, (xi-1)/2) += w(0)*qdoty;
                    mac_y(yi/2 + 1, (xi-1)/2) += w(1)*qdoty;
                    mac_y(yi/2, (xi-1)/2 + 1) += w(2)*qdoty;
                    mac_y(yi/2 + 1, (xi-1)/2 +1) += w(3)*qdoty;
                }
            }
        }
    }
 
    Eigen::SparseMatrix<double, Eigen::RowMajor> A(num_cells_x*num_cells_y, num_cells_x*num_cells_y);
    Eigen::VectorXd b(num_cells_x*num_cells_y);


    for (int i=0; i < num_cells_x*num_cells_y; i++) {

            int num_neighbours = 0;
            int yval = i/num_cells_x;
            int xval = i%num_cells_x;
            // assumes that GRID_DX === GRID_DY
            b(i) = (mac_x(yval, xval+1) - mac_x(yval, xval) + mac_y(yval+1,xval) - mac_y(yval, xval))/GRID_DX;

            if (xval != 0) { // insert to the left
                A.insert(i, i-1) = -1;
                num_neighbours++;
            }

            if (xval != num_cells_x-1) { // insert to the right
                A.insert(i, i+1) = -1;
                num_neighbours++;
            }
        
            if (yval != 0){  // insert below
                A.insert(i, i-num_cells_x) = -1;
                num_neighbours++;
            }

            if (yval != num_cells_y -1) { // insert above
                A.insert(i, i+num_cells_x) = -1;
                num_neighbours++;
            }

            A.insert(i,i) = num_neighbours;
    }
    Eigen::VectorXd p = Eigen::VectorXd::Zero(num_cells_x*num_cells_y);
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double>,
                              Eigen::Lower|Eigen::Upper                              
                              > cg;
    
    cg.setMaxIterations(64);
    cg.compute((dt/(density * GRID_DX)) * A);
    
    /* TODO: should we negate the divergence vectors? */
    p = cg.solve(b);

    /* update grid velocities with computed pressures */
    for (int pi = 0; pi < num_cells_x*num_cells_y; pi++) {
        int py = pi/num_cells_x;
        int px = pi%num_cells_x;

        if (py != 0) {
            double dv = p(py*num_cells_x+px) - p((py-1)*num_cells_x+px);
            dv *= dt/(density*GRID_DX);
            mac_y(py,px) -= dv;
        }

        if (px != 0) {
            double du = p(py*num_cells_x+px) - p(py*num_cells_x +px -1);
            du *= dt/(density*GRID_DX);
            mac_x(py, px) -= du;
        }
    }

    /* convert back to particle form */
    for (int yi = 0; yi < 2*num_cells_y; yi++) {
        for (int xi = 0; xi < 2*num_cells_x; xi++) {
            for (int p=0; p < cells[yi][xi].size(); p++) {

                double x_v_interpolation, y_v_interpolation;
                
                int index = cells[xi][yi][p];                
                double qdotx = qdot(2*index);
                double qdoty = qdot(2*index+1); 
                double qx = q(2*index);
                double qy = q(2*index+1);

                // PIC transfer from grid to particle
                if (yi == 0 || yi == (2*num_cells_y)-1) { // only do linear interpolation? is this right?
                    x_v_interpolation = linear_interpolate(
                        qx, 
                        (xi/2)*GRID_DX, (xi/2 +1)*GRID_DX, 
                        mac_x(yi/2,xi/2), mac_x(yi/2, xi/2 +1)
                    );

                } else {
                    x_v_interpolation = bilinear_interpolate(
                        qx, qy,
                        (xi/2)*GRID_DX, (xi/2 +1)*GRID_DX,
                        (yi/2)*GRID_DY, (yi/2 +1)*GRID_DY,
                        mac_x(yi/2, xi/2), mac_x(yi/2, xi/2 +1),
                        mac_x(yi/2 +1, xi/2), mac_x(yi/2+1, xi/2 +1)
                    );
                }

                if (xi == 0 || xi == (2*num_cells_x)-1) {
                    y_v_interpolation = linear_interpolate(
                        qy,
                        (yi/2)*GRID_DY, (yi/2 +1)*GRID_DY,
                        mac_y(yi/2,xi/2), mac_y(yi/2 + 1, xi/2)
                    );
                } else {
                    y_v_interpolation = bilinear_interpolate(
                        qx, qy,
                        (xi/2)*GRID_DX, (xi/2 +1)*GRID_DX,
                        (yi/2)*GRID_DY, (yi/2 +1)*GRID_DY,
                        mac_y(yi/2, xi/2), mac_y(yi/2, xi/2+1),
                        mac_y(yi/2 +1, xi/2), mac_y(yi/2+1, xi/2 +1)
                    );
                }

                qdot(2*index) = x_v_interpolation;
                qdot(2*index+1) = y_v_interpolation;
            }
        }
    }
    
}


/* Helpers */
void bilinear_weights(
    Eigen::Ref<Eigen::Vector4d> weights,
    double x, double y,
    double x1, double x2, double y1, double y2) 
{
    double wx1 = (x2 - x) / (x2 - x1);
    double wx2 = (x - x1) / (x2 - x1);
    double wy1 = (y2 - y) / (y2 - y1);
    double wy2 = (y - y1) / (y2 - y1);

    weights(0) = wx1 * wy1;
    weights(1) = wx1 * wy2;
    weights(2) = wx2 * wy1;
    weights(3) = wx2 * wy2;
}
double bilinear_interpolate(
    double x, double y,
    double x1, double x2, double y1, double y2,
    double q11, double q12, double q21, double q22) 
{    
    Eigen::Vector4d w;
    bilinear_weights(w, x, y, x1, x2, y1, y2);

    return w(0)*q11 + w(1)*q12 + w(2)*q21 + w(3)*q22; 
}

void linear_weights(
    Eigen::Ref<Eigen::Vector2d> weights,
    double x,
    double x1, double x2) 
{    
    weights(0) = (x2 - x) / (x2 - x1);
    weights(1) = (x - x1) / (x2 - x1);
}

double linear_interpolate(
    double x,
    double x1, double x2,
    double q1, double q2) 
{
    Eigen::Vector2d w;
    linear_weights(w, x, x1, x2);
    return q1*w(0) + q2*w(1);
}
