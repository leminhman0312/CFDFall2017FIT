#include "explicit.h"

#include <cstdio>

std::vector<std::vector<double>> FTCS_Explicit(
    const std::vector<std::vector<double>>& u0, int nmax, double deltax,
    double deltay, double dt, double alpha, double t1, double t2, double t3,
    double t4) {
  // create u and pass u0 from init
  std::vector<std::vector<double>> u = u0;

  int imax = (int)u.size() - 1;
  int jmax = (int)u[1].size() - 1;

  std::vector<std::vector<double>> u_new(imax + 1,
                                         std::vector<double>(jmax + 1, 0.0));

  double fx = (alpha * dt) / (deltax * deltax);
  double fy = (alpha * dt) / (deltay * deltay);

  // optional stability warning (2D explicit can blow up)
  if (fx + fy > 0.5) {
    printf("WARNING: FTCS explicit may be unstable: fx+fy = %.6f > 0.5\n",
           fx + fy);
  }

  for (int n = 1; n <= nmax; n++) {
    // carry everything over first
    for (int i = 1; i <= imax; i++) {
      for (int j = 1; j <= jmax; j++) {
        u_new[i][j] = u[i][j];
      }
    }

    // interior update
    for (int i = 2; i <= imax - 1; i++) {
      for (int j = 2; j <= jmax - 1; j++) {
        u_new[i][j] = (1.0 - 2.0 * fx - 2.0 * fy) * u[i][j] +
                      fx * (u[i + 1][j] + u[i - 1][j]) +
                      fy * (u[i][j + 1] + u[i][j - 1]);
      }
    }

    // enforce BCs
    for (int j = 1; j <= jmax; j++) {
      u_new[1][j] = t2;
      u_new[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u_new[i][1] = t1;
      u_new[i][jmax] = t3;
    }

    u.swap(u_new);
  }

  return u;
}