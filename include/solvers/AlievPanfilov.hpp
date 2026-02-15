#pragma once
#include <Eigen/Dense>

// We have to create our model, to be efficient we create a class with our parameters, instructions and so on
// Here we'll use canine myocardium to do the basic experience

namespace CardiacArrhySim{

class AlievPanfilov {

    // User can see and use these parameters
    public :
        // our bag filled with parameters. 
        // Note that we put here our struct in a class. Thanks to this encapsulation we avoid confusion with other param.
        struct Parameters
        {
            double k = 8.0;
            double a = 0.15;
            double eps0 = 0.002;
            double mu1 = 0.2;
            double mu2 = 0.3;

        };
    
    // In our lecture it's : AlievPanfilov() {}
    AlievPanfilov() = default; // If we don't do this, the compiler while create one but it can be deleted if we create a other constructor
    
    // Here we have something new compared with our lecture, let's see the meaning of this line
    // nodiscard says to the compiler that if someone call this function but don't store the output it gives a warning because of compute_result(vector)
    // Indeed we create thanks to eigen a vector of size X which contains doubles and this is the type of our function compute_result
    // our function has a argument which is a vector that cannot be modified (const), only readed, and we use it's memory adress, we avoid to copy the vector 
    [[nodiscard]] Eigen::VectorXd compute_result(const Eigen::VectorXd& y) const; // the const at the end allows to do some calculus without change the internal state of our class.

    // Plus we want to do some OpenMP to optimize our cost, so we create a parallel method :
    [[nodiscard]] Eigen::VectorXd compute_result_parallel(const Eigen::VectorXd& y) const;

    // But if we have several nodes, it seems logic to use in an intelligent way our methods
    [[nodiscard]] Eigen::VectorXd compute_choice(const Eigen::VectorXd& y) const
    {
        if (y.size() < 2000) return compute_result(y);
        else return compute_result_parallel(y);
    }
    

    private :
        Parameters p_; // Only functions of our class can touch p_, for example if we want to set k for a calculation we can do it

};

}