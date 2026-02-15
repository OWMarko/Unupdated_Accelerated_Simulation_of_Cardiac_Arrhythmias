#include "include/solvers/AlievPanfilov.hpp"

// Now that we have our header, we need to create our equation and to modulate it to be 0D,1D and 2D. 

Eigen::VectorXd CardiacArrhySim::AlievPanfilov::compute_result(const Eigen::VectorXd& y) const {
    
    // First we check size of our vector, we have to deal with a 2 var syst so the size is pair
    long total_size = y.size();
    if (total_size % 2 != 0)
    {
        throw std::runtime_error("Our vector muste be U + V (our 2 states, time and space)");
    }
    long N = total_size / 2;

    Eigen::VectorXd dydt(total_size); // same size as y, it's the temp. der.

    // Now we play with our vectors, we know that first terms of u starts at the beginning of the vect. and v start at the end of u (start + x)

    const double* u_ptr = y.data();
    const double* v_ptr = y.data() + N;

    double* du_ptr = dydt.data();
    double* dv_ptr = dydt.data() + N;

    const double k = p_.k;
    const double a = p_.a;
    const double eps0 = p_.eps0;
    const double mu1 = p_.mu1;
    const double mu2 = p_.mu2;

    for (long i = 0; i < N; ++i)
    {
        double u = u_ptr[i];
        double v = v_ptr[i];

        double epsilon = eps0 + (mu1*v) / (u+mu2);

        du_ptr[i] = k *u*(u-a)*(u-1) - u*v;
        dv_ptr[i] = epsilon * (-v-k*u*(u - a - 1.0));

    }

    return dydt;
}

Eigen::VectorXd CardiacArrhySim::AlievPanfilov::compute_result_parallel(const Eigen::VectorXd& y) const{


    long total_size = y.size();
    if(total_size % 2 != 0)
    {
        throw std::runtime_error("Our vector muste be U + V (our 2 states, time and space)");

    }

    long N = total_size / 2;

    Eigen::VectorXd dydt(total_size);

    const double* u_ptr = y.data();
    const double* v_ptr = y.data() + N;

    double* du_ptr = dydt.data();
    double* dv_ptr = dydt.data() + N;

    const double k = p_.k;
    const double a = p_.a;
    const double eps0 = p_.eps0;
    const double mu1 = p_.mu1;
    const double mu2 = p_.mu2;

    #pragma omp parallel for default(none) \
        shared(u_ptr,v_ptr,du_ptr,dv_ptr,k,a,eps0,mu1,mu2) \
        schedule(static)
        for(long i = 0;i < N; ++i)
        {
            double u = u_ptr[i];
            double v = v_ptr[i];

            double epsilon = eps0 + (mu1*v) / (u+mu2);

            du_ptr[i] = k *u*(u-a)*(u-1) - u*v;
            dv_ptr[i] = epsilon * (-v-k*u*(u - a - 1.0));
    }

    return dydt;
}

