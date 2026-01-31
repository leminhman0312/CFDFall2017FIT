#pragma once
#include <headers.h>

std::vector<std::vector<double>> FTCS_implicit(
    const std::vector<std::vector<double>>& u0, int nmax, double deltax,
    double deltay, double dt, double alpha, double t1, double t2, double t3,
    double t4);
