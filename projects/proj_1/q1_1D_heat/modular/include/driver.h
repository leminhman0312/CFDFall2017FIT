#ifndef DRIVER_H
#define DRIVER_H

#include "headers.h"

// prototypes
void sim(double delta_x, double delta_t, int imax, double t0, double tboundary,
         double F_num);
void plot(double delta_t);
void compare_schemes(double delta_t, double t_target = -1.0);
void compare_dt(const char* scheme, double dt1, double dt2, double t_target);
void compare_dt_all_times(const char* scheme, double dt1, double dt2);

void compare_error_schemes(double delta_t, double t_target);

void convergence_study(double delta_x, int imax, double alpha, double t0,
                       double tboundary, double t_target,
                       const double dt_list[], int ndt);

#endif
