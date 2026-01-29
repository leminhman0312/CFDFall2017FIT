#include "implicit.h"

#include "thomas.h"

std::vector<std::vector<double>> FTCS_implicit(
    const std::vector<std::vector<double>>& u0, int nmax, double deltax,
    double deltay, double dt, double alpha, double t0, double t1, double t2,
    double t3, double t4) {
  // create u and pass u0 from init
  std::vector<std::vector<double>> u = u0;

  int imax = (int)u.size() - 1;
  int jmax = (int)u[1].size() - 1;

  std::vector<std::vector<double>> u_dummy(imax + 1,
                                           std::vector<double>(jmax + 1, 0.0));

  double fx = (alpha * dt) / (deltax * deltax);
  double fy = (alpha * dt) / (deltay * deltay);

  // DEFINE VECTORS FOR THOMAS IN X
  double ax[imax + 1];  // above
  double bx[imax + 1];  // below
  double cx[imax + 1];  // rhs
  double dx[imax + 1];  // diagonal

  // Thomas arrays interior points
  for (int i = 2; i <= imax; i++) {
    ax[i] = -fx / 2.;
    bx[i] = -fx / 2.;
    dx[i] = 1.0 + fx;
  }

  // endpoints
  dx[1] = 1.0;
  ax[1] = 0.0;
  bx[1] = 0.0;

  dx[imax] = 1.0;
  ax[imax] = 0.0;
  bx[imax] = 0.0;

  // DEFINE VECTORS FOR THOMAS IN Y
  double ay[jmax + 1];  // above
  double by[jmax + 1];  // below
  double cy[jmax + 1];  // rhs
  double dy[jmax + 1];  // diagonal

  // Thomas arrays interior points
  for (int j = 2; j <= jmax; j++) {
    ay[j] = -fy / 2.;
    by[j] = -fy / 2.;
    dy[j] = 1.0 + fy;
  }

  // endpoints
  dy[1] = 1.0;
  ay[1] = 0.0;
  by[1] = 0.0;

  dy[jmax] = 1.0;
  ay[jmax] = 0.0;
  by[jmax] = 0.0;

  for (int n = 1; n <= nmax; n++) {
    // enforce BCs on u at start of step (good habit)
    for (int j = 1; j <= jmax; j++) {
      u[1][j] = t2;     // top
      u[imax][j] = t4;  // bottom
    }
    for (int i = 1; i <= imax; i++) {
      u[i][1] = t1;     // left
      u[i][jmax] = t3;  // right
    }

    // -------------------------
    // STEP 1: X sweep, solve for u_dummy = T^(n+1/2)
    // For each j, solve tridiagonal in i
    // -------------------------
    for (int j = 2; j <= jmax - 1; j++) {
      // Build RHS cx for interior i
      for (int i = 2; i <= imax - 1; i++) {
        cx[i] = (1.0 - fy) * u[i][j] + (fy / 2.0) * (u[i][j + 1] + u[i][j - 1]);
      }

      // Add boundary contributions in x into RHS
      // At i=2, term involves u[1][j] which is known boundary
      cx[2] += (fx / 2.0) * u[1][j];
      // At i=imax-1, term involves u[imax][j]
      cx[imax - 1] += (fx / 2.0) * u[imax][j];

      // Thomas solve along i direction
      // Solution goes into a temp 1D line, then copy to u_dummy
      // Reuse cx as RHS, and store solution in a line array solx
      double solx[imax + 1] = {0.0};

      // endpoints are boundaries
      solx[1] = u[1][j];
      solx[imax] = u[imax][j];

      cx[1] = u[1][j];
      cx[imax] = u[imax][j];

      thomasTriDiagonal(imax, ax, bx, cx, dx, solx);

      // copy solved line into u_dummy
      for (int i = 1; i <= imax; i++) {
        u_dummy[i][j] = solx[i];
      }
    }

    // enforce BCs on u_dummy too
    for (int j = 1; j <= jmax; j++) {
      u_dummy[1][j] = t2;
      u_dummy[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u_dummy[i][1] = t1;
      u_dummy[i][jmax] = t3;
    }

    // -------------------------
    // STEP 2: Y sweep, solve for u = T^(n+1)
    // For each i, solve tridiagonal in j
    // -------------------------
    for (int i = 2; i <= imax - 1; i++) {
      // Build RHS cy for interior j
      for (int j = 2; j <= jmax - 1; j++) {
        cy[j] = (1.0 - fx) * u_dummy[i][j] +
                (fx / 2.0) * (u_dummy[i + 1][j] + u_dummy[i - 1][j]);
      }

      // Add boundary contributions in y into RHS
      cy[2] += (fy / 2.0) * u_dummy[i][1];
      cy[jmax - 1] += (fy / 2.0) * u_dummy[i][jmax];

      // Thomas solve along j direction
      double soly[jmax + 1] = {0.0};
      soly[1] = u_dummy[i][1];
      soly[jmax] = u_dummy[i][jmax];

      thomasTriDiagonal(jmax, ay, by, cy, dy, soly);

      // copy solved column into u
      for (int j = 1; j <= jmax; j++) {
        u[i][j] = soly[j];
      }
    }

    // enforce final BCs on u
    for (int j = 1; j <= jmax; j++) {
      u[1][j] = t2;
      u[imax][j] = t4;
    }
    for (int i = 1; i <= imax; i++) {
      u[i][1] = t1;
      u[i][jmax] = t3;
    }
  }
  return u;
}