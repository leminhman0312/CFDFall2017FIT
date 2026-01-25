#ifndef METHODS_H
#define METHODS_H

#include <cmath>
#include <cstdio>

// prototypes

void FTCS_explicit(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]);

void DufortFrankel(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]);

void FTCS_implicit(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]);

void CrankNicolson(int nmax, double F_num, double tboundary, double t0,
                   int imax, double out_un[]);

void thomasTriDiagonal(int imax, double a[], double b[], double c[], double d[],
                       double u[]);

void FTCS_return_implicit_interior(int nmax, double F_num, double tboundary,
                                   double t0, int imax, double u_full[]);

void exact_solution_profile(int imax, const double x[], double t, double alpha,
                            double L, double tboundary, double t0, int nterms,
                            double out_T[]);

#endif
