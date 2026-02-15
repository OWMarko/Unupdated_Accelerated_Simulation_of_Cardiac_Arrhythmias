#pragma once
#include "Methods.hpp"

namespace CardiacArrhySim{

class PFE : public Methods{

    public : 

    PFE(double M, int k_micro) : M_(M), k_micro_(k_micro) {} // We create a constructor for our class, we want to be able to set M and k_micro when we create an object of our class

    [[nodiscard]] Eigen::VectorXd steps(const AlievPanfilov& model, const Eigen::VectorXd& y, double dt) const override;

    [[nodiscard]] std::string name() const override {return "PFE";}

    private:
    double M_;
    int k_micro_;

};
}