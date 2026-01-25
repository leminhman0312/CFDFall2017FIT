#include <cmath>
#include "driver.h"


int main()
{
    int imax = 21; // domain length
    double delta_x = 0.05; //spacing in x 
    double delta_t = 0.05; //spacing in time

    double alpha = 0.1; // thermal diffusivity 
    double t0 = 100.0; // initial temperature
    double tboundary = 300.0; // boundary side temperatures

    double Fourier_num = (alpha*delta_t)/(pow(delta_x,2.0)); // alpha*dt/dx^2 the Fourier number

    sim(delta_x,delta_t,imax,t0,tboundary,Fourier_num);
    plot(delta_t);

    compare_schemes(0.01,0.1);

    compare_dt_all_times("ftcs_explicit", 0.01, 0.05);
    compare_dt_all_times("dufort",        0.01, 0.05);
    compare_dt_all_times("ftcs_implicit", 0.01, 0.05);
    compare_dt_all_times("cn",            0.01, 0.05);

    return 0;
}

