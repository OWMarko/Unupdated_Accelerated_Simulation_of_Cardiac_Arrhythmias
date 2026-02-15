#pragma once
#include "include/solvers/AlievPanfilov.hpp"
#include "include/solvers/Methods.hpp"

namespace CardiacArrhySim{

class RK4 : public Methods{

    public : 
    [[nodiscard]] Eigen::VectorXd steps(const AlievPanfilov& model, const Eigen::VectorXd& y, double dt) const override;

    [[nodiscard]] std::string name() const override { return "Runge Kutta 4"; }

};
}