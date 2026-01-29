#include "thomas.h"

void thomasTriDiagonal(int imax, double a[], double b[], double c[], double d[],
                       double u[]) {
  double dprime[imax + 1];
  double cprime[imax + 1];

  dprime[1] = d[1];
  cprime[1] = c[1];

  for (int i = 2; i <= imax; i++) {
    dprime[i] = d[i] - (b[i] * a[i - 1]) / dprime[i - 1];
    cprime[i] = c[i] - (cprime[i - 1] * b[i]) / dprime[i - 1];
  }

  u[imax] = cprime[imax] / dprime[imax];

  for (int i = imax - 1; i >= 1; i--) {
    u[i] = (cprime[i] - a[i] * u[i + 1]) / dprime[i];
  }
}