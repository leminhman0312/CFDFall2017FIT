#include "methods.h"
#include <cmath>

// FTCS explicit

void FTCS_explicit(int nmax, double F_num, double tboundary, double t0,int imax, double out_un[]){
    double u[imax+1];
    double un[imax+1];

    for (int i = 2; i <= imax-1; i++){
        u[i]  = t0;
        un[i] = t0;
    }

    u[1] = tboundary;   u[imax] = tboundary;
    un[1] = tboundary;  un[imax] = tboundary;

    for (int n = 1; n <= nmax; n++){
        for (int i = 1; i <= imax; i++) u[i] = un[i];

        for (int i = 2; i <= imax-1; i++){
            un[i] = u[i] + F_num*(u[i+1] - 2.0*u[i] + u[i-1]);
        }

        un[1] = tboundary;
        un[imax] = tboundary;
    }

    for (int i = 1; i <= imax; i++) out_un[i] = un[i];
}

// Dufort Frankel
void DufortFrankel(int nmax, double F_num, double tboundary, double t0,int imax, double out_un[]){
    double d = 2.0*F_num;

    double un_m1[imax+1] = {0.0};
    double un_p1[imax+1] = {0.0};
    double un[imax+1];

    un_m1[1] = tboundary;
    un_m1[imax] = tboundary;
    for (int i = 2; i <= imax-1; i++) un_m1[i] = t0;

    // starter step using your implicit interior helper
    double u1[imax+1];
    FTCS_return_implicit_interior(1, F_num, tboundary, t0, imax, u1);
    for (int i = 1; i <= imax; i++) un[i] = u1[i];

    for (int n = 1; n <= nmax-1; n++){
        for (int i = 2; i <= imax-1; i++){
            un_p1[i] = ((1.0 - d)*un_m1[i] + d*(un[i+1] + un[i-1])) / (1.0 + d);
        }

        un_p1[1] = tboundary;
        un_p1[imax] = tboundary;

        for (int i = 1; i <= imax; i++){
            un_m1[i] = un[i];
            un[i] = un_p1[i];
        }
    }

    for (int i = 1; i <= imax; i++) out_un[i] = un[i];
}

// implicit interior helper
void FTCS_return_implicit_interior(int nmax, double F_num, double tboundary, double t0,int imax, double u_full[]){
    
    int N = imax - 2;

    double a[N+1], b[N+1], c[N+1], d[N+1], u_int[N+1];

    for (int i = 0; i <= N; i++){
        a[i]=b[i]=c[i]=d[i]=u_int[i]=0.0;
    }

    u_full[1] = tboundary;
    u_full[imax] = tboundary;
    for (int i = 2; i <= imax-1; i++) u_full[i] = t0;

    // uniform interior system, then enforce edge structure
    for (int j = 1; j <= N; j++){
        d[j] = 1.0 + 2.0*F_num;
        a[j] = -F_num;
        b[j] = -F_num;
    }
    b[1] = 0.0;
    a[N] = 0.0;

    for (int t = 1; t <= nmax; t++){
        for (int j = 1; j <= N; j++){
            int i = j + 1;
            c[j] = u_full[i];
        }

        c[1] += F_num*tboundary;
        c[N] += F_num*tboundary;

        thomasTriDiagonal(N, a, b, c, d, u_int);

        u_full[1] = tboundary;
        u_full[imax] = tboundary;
        for (int j = 1; j <= N; j++){
            int i = j + 1;
            u_full[i] = u_int[j];
        }
    }
}

// full implicit solver
void FTCS_implicit(int nmax, double F_num, double tboundary, double t0,int imax, double out_un[]){

    double a[imax+1], b[imax+1], c[imax+1], d[imax+1], un[imax+1];

    un[1] = tboundary;
    un[imax] = tboundary;
    for (int i = 2; i <= imax-1; i++) un[i] = t0;

    d[1] = 1.0; a[1]=0.0; b[1]=0.0;
    d[imax] = 1.0; a[imax]=0.0; b[imax]=0.0;

    for (int i = 2; i <= imax-1; i++){
        d[i] = 1.0 + 2.0*F_num;
        a[i] = -F_num;
        b[i] = -F_num;
    }

    for (int t = 1; t <= nmax; t++){
        c[1] = tboundary;
        c[imax] = tboundary;
        for (int i = 2; i <= imax-1; i++) c[i] = un[i];
        thomasTriDiagonal(imax, a, b, c, d, un);
    }

    for (int i = 1; i <= imax; i++) out_un[i] = un[i];
}

// Crank Nicolson (your current logic)
void CrankNicolson(int nmax, double F_num, double tboundary, double t0,int imax, double out_un[]){

    double d_NC = F_num / 2.0;

    double a[imax+1], b[imax+1], c[imax+1], d[imax+1];
    double u0[imax+1];
    double u_half[imax+1];

    u0[1] = tboundary;
    u0[imax] = tboundary;
    for (int i = 2; i <= imax-1; i++) u0[i] = t0;

    u_half[1] = tboundary;
    u_half[imax] = tboundary;

    d[1]=1.0; a[1]=0.0; b[1]=0.0;
    d[imax]=1.0; a[imax]=0.0; b[imax]=0.0;

    for (int i = 2; i <= imax-1; i++){
        d[i] = 1.0 + 2.0*d_NC;
        a[i] = -d_NC;
        b[i] = -d_NC;
    }

    for (int n = 1; n <= nmax; n++){
        u_half[1] = tboundary;
        u_half[imax] = tboundary;

        for (int i = 2; i <= imax-1; i++){
            u_half[i] = u0[i] + d_NC*(u0[i+1] - 2.0*u0[i] + u0[i-1]);
        }

        c[1] = tboundary;
        c[imax] = tboundary;
        for (int i = 2; i <= imax-1; i++) c[i] = u_half[i];

        thomasTriDiagonal(imax, a, b, c, d, u0);

        u0[1] = tboundary;
        u0[imax] = tboundary;
    }

    for (int i = 1; i <= imax; i++) out_un[i] = u0[i];
}

// thomas
void thomasTriDiagonal(int imax, double a[], double b[], double c[], double d[], double u[]){
    double dprime[imax+1];
    double cprime[imax+1];

    dprime[1] = d[1];
    cprime[1] = c[1];

    for (int i = 2; i <= imax; i++){
        dprime[i] = d[i] - (b[i]*a[i-1]) / dprime[i-1];
        cprime[i] = c[i] - (cprime[i-1]*b[i]) / dprime[i-1];
    }

    u[imax] = cprime[imax] / dprime[imax];

    for (int i = imax-1; i >= 1; i--){
        u[i] = (cprime[i] - a[i]*u[i+1]) / dprime[i];
    }
}
