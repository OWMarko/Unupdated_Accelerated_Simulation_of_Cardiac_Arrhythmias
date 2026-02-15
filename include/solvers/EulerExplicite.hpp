#pragma once
#include "include/solvers/AlievPanfilov.hpp"
#include "include/solvers/Methods.hpp"

namespace CardiacArrhySim{

class EE : public Methods {

    public : 
    [[nodiscard]] Eigen::VectorXd steps(const AlievPanfilov& model, const Eigen::VectorXd& y, double dt) const override; // We want to implement the steps method for our ExpliciteEuler method, we use override to indicate that we are overriding a virtual method from the base class

    [[nodiscard]] std::string name() const override { return "Euler Explicite"; }
};
}