#pragma once
#include "AlievPanfilov.hpp"
#include <string>
#include <Eigen/Dense>

namespace CardiacArrhySim{

class Methods{

    public : 
        virtual ~Methods() = default; // We want to be able to destroy our class without problem so we create a virtual destructor

        [[nodiscard]] virtual Eigen::VectorXd steps(const AlievPanfilov& model, const Eigen::VectorXd& y, double dt) const = 0; // We want to create a method that will be used by all our methods, but we don't know how to do it, so we create a pure virtual method, this is the meaning of = 0, it means that this method is pure virtual and must be implemented by all the classes that inherit from Methods

        [[nodiscard]] virtual std::string name() const = 0; // We want to know the name of our method, so we create a pure virtual method that will return the name of our method
};

}
