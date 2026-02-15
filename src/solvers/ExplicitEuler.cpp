#include "include/solvers/EulerExplicite.hpp"

namespace CardiacArrhySim{

Eigen::VectorXd EE::steps(const AlievPanfilov& model, const Eigen::VectorXd& y, double dt) const
{
    return y + dt*model.compute_choice(y);
}
}