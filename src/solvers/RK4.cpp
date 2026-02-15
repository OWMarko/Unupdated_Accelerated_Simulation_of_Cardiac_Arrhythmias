#include "include/solvers/RK4.hpp"

namespace CardiacArrhySim{

Eigen::VectorXd RK4::steps(const AlievPanfilov& model, Eigen::VectorXd& y, double dt) const
{
    Eigen::VectorXd k1 = model.compute_choice(y);
    Eigen::VectorXd k2 = model.compute_choice(y + 0.5 * dt * k1);
    Eigen::VectorXd k3 = model.compute_choice(y + 0.5 * dt * k2);
    Eigen::VectorXd k4 = model.compute_choice(y + dt * k3);

    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

}