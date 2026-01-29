// MAE 5150: Coding Project 1
// Max Le
// Question 2: 2D Heat Conduction Equation

//----------------------------------------------------------------------------------//

#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

// Functions Prototypes

// Put this near the top (prototype)
std::vector<std::vector<double>> initializeField(int imax, int jmax, double t0,
                                                 double t1, double t2,
                                                 double t3, double t4);

std::vector<std::vector<double>> FTCS_implicit(
    const std::vector<std::vector<double>>& u0, int nmax, double deltax,
    double deltay, double dt, double alpha, double t0, double t1, double t2,
    double t3, double t4);

std::vector<std::vector<double>> FTCS_Explicit(
    const std::vector<std::vector<double>>& u0, int nmax, double deltax,
    double deltay, double dt, double alpha, double t0, double t1, double t2,
    double t3, double t4);

void thomasTriDiagonal(int imax, double a[], double b[], double c[], double d[],
                       double u[]);

void write_field_xyz(const char* filename,
                     const std::vector<std::vector<double>>& u, double deltax,
                     double deltay);

void plotContour2D(const char* infile, const char* outpng, double xmax,
                   double ymax, const char* tlabel);

void plotContourMatlabLike(const char* datafile, const char* outpng);

void ensure_dir(const char* name);

void ensure_project_dirs();

#include <cmath>   // std::lround
#include <cstdio>  // optional

//----------------------------------------------------------------------------------//
// Main code
int main() {
  ensure_project_dirs();

  // time controls
  const double dt = 0.01;    // hr
  const double hours = 0.5;  // hr
  const int nmax = (int)std::lround(hours / dt);

  // grid
  const double deltax = 0.1;  // ft
  const double deltay = 0.1;  // ft

  // physics
  const double alpha = 0.645;  // (units consistent with your equation)

  // domain
  const double minlength_X = 0.0;
  const double maxlength_X = 3.5;
  const double minlength_Y = 0.0;
  const double maxlength_Y = 3.5;
  int imax = (int)ceil(((maxlength_X - minlength_X) / deltax) + 1.0);
  int jmax = (int)ceil(((maxlength_Y - minlength_Y) / deltay) + 1.0);

  // temperatures
  const double t0 = 0.0;    // interior initial
  const double t1 = 200.0;  // bottom boundary
  const double t2 = 200.0;  // left boundary
  const double t3 = 0.0;    // top boundary
  const double t4 = 0.0;    // right boundary

  // calculate dt for explicit scheme: fx + fy <= 0.5
  const double dt_explicit =
      0.5 / (alpha * ((1. / pow(deltax, 2)) + (1. / (pow(deltay, 2)))));

  // pre processing
  auto u0 = initializeField(imax, jmax, t0, t1, t2, t3, t4);

  // simulation
  auto u_explicit = FTCS_Explicit(u0, nmax, deltax, deltay, dt_explicit, alpha,
                                  t0, t1, t2, t3, t4);

  auto u_implicit =
      FTCS_implicit(u0, nmax, deltax, deltay, dt, alpha, t0, t1, t2, t3, t4);

  // post processing
  write_field_xyz("data/initial.dat", u0, deltax, deltay);
  plotContourMatlabLike("data/initial.dat", "data/contour_initial.png");

  write_field_xyz("data/implicit.dat", u_implicit, deltax, deltay);
  plotContourMatlabLike("data/implicit.dat", "data/contour_implicit.png");

  write_field_xyz("data/explicit.dat", u_explicit, deltax, deltay);
  plotContourMatlabLike("data/explicit.dat", "data/contour_explicit.png");
  return 0;
}

//----------------------------------------------------------------------------------//

// Functions

std::vector<std::vector<double>> FTCS_Explicit(
    const std::vector<std::vector<double>>& u0, int nmax, double deltax,
    double deltay, double dt, double alpha, double t0, double t1, double t2,
    double t3, double t4) {
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
      u_new[1][j] = t3;
      u_new[imax][j] = t1;
    }
    for (int i = 1; i <= imax; i++) {
      u_new[i][1] = t2;
      u_new[i][jmax] = t4;
    }

    u.swap(u_new);
  }

  return u;
}

std::vector<std::vector<double>> initializeField(int imax, int jmax, double t0,
                                                 double t1, double t2,
                                                 double t3, double t4) {
  // allocate and fill everything with t0 first
  // initial conditions
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

  // // optional: explicit corner rule (average the two touching sides)
  // u[1][1] = 0.5 * (t2 + t3);        // top-left
  // u[1][jmax] = 0.5 * (t4 + t3);     // top-right
  // u[imax][1] = 0.5 * (t2 + t1);     // bottom-left
  // u[imax][jmax] = 0.5 * (t4 + t1);  // bottom-right

  return u;
}

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
      u[1][j] = t3;     // top
      u[imax][j] = t1;  // bottom
    }
    for (int i = 1; i <= imax; i++) {
      u[i][1] = t2;     // left
      u[i][jmax] = t4;  // right
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
      u_dummy[1][j] = t3;
      u_dummy[imax][j] = t1;
    }
    for (int i = 1; i <= imax; i++) {
      u_dummy[i][1] = t2;
      u_dummy[i][jmax] = t4;
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
      u[1][j] = t3;
      u[imax][j] = t1;
    }
    for (int i = 1; i <= imax; i++) {
      u[i][1] = t2;
      u[i][jmax] = t4;
    }
  }
  return u;
}

// thomas
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

void write_field_xyz(const char* filename,
                     const std::vector<std::vector<double>>& u, double deltax,
                     double deltay) {
  int imax = (int)u.size() - 1;
  int jmax = (int)u[1].size() - 1;

  FILE* f = fopen(filename, "w");
  if (!f) return;

  for (int j = 1; j <= jmax; j++) {
    double y = (j - 1) * deltay;
    for (int i = 1; i <= imax; i++) {
      double x = (i - 1) * deltax;
      fprintf(f, "%.8f %.8f %.10f\n", x, y, u[i][j]);
    }
    fprintf(f, "\n");
  }

  fclose(f);
}

void plotContour2D(const char* infile, const char* outpng, double xmax,
                   double ymax, const char* tlabel) {
  char cmd[1024];

  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \"infile='%s'; outpng='%s'; xmax=%.6f; ymax=%.6f; "
                "tlabel='%s'\" "
                "gnuplot_scripts/plot_contour_2d.gp",
                infile, outpng, xmax, ymax, tlabel);

  system(cmd);
}

void plotContourMatlabLike(const char* datafile, const char* outpng) {
  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \"datafile=\\\"%s\\\"; outpng=\\\"%s\\\"\" "
                "gnuplot_scripts/plot_contour_2d_matlab_like.gp",
                datafile, outpng);
  printf("GNUPLOT CMD:\n%s\n", cmd);
  system(cmd);
}

void ensure_dir(const char* name) {
  struct stat st;
  if (stat(name, &st) != 0) {
    mkdir(name, 0755);
  }
}

void ensure_project_dirs() {
  ensure_dir("data");
  ensure_dir("plot");
  ensure_dir("gnuplot_scripts");
}
