#include <cmath>
#include <vector>

#include "explicit.h"
#include "filesystem.h"
#include "implicit.h"
#include "initialize.h"
#include "io.h"
#include "plotting.h"

int main() {
  ensure_project_dirs();

  // time controls
  const double dt = 0.01;    // hr
  const double hours = 1.0;  // hr
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
  plotContourMatlabLike("data/initial.dat", "plot/contour_initial.png", hours,
                        "Initial conditions");

  write_field_xyz("data/implicit.dat", u_implicit, deltax, deltay);
  plotContourMatlabLike("data/implicit.dat", "plot/contour_implicit.png", hours,
                        "FTCS implicit");

  write_field_xyz("data/explicit.dat", u_explicit, deltax, deltay);
  plotContourMatlabLike("data/explicit.dat", "plot/contour_explicit.png", hours,
                        "FTCS Explicit");
  return 0;
}