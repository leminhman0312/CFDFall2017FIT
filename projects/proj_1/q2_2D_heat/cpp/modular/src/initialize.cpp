#include "initialize.h"

std::vector<std::vector<double>> initializeField(int imax, int jmax, double t0,
                                                 double t1, double t2,
                                                 double t3, double t4) {
  std::vector<std::vector<double>> u(imax + 1,
                                     std::vector<double>(jmax + 1, t0));

  // boundaries
  for (int j = 1; j <= jmax; j++) {
    u[1][j] = t2;     // left
    u[imax][j] = t4;  // right
  }
  for (int i = 1; i <= imax; i++) {
    u[i][1] = t1;     // bottom
    u[i][jmax] = t3;  // top
  }

  return u;
}
