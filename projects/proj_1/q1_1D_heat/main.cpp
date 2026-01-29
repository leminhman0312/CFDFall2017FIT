#include <cmath>

#include "driver.h"

int main() {
  int imax = 21;             // domain length
  double delta_x = 0.05;     // spacing in x
  double delta_t = 0.01;     // spacing in time
  double alpha = 0.1;        // thermal diffusivity
  double t0 = 100.0;         // initial temperatures
  double tboundary = 300.0;  // boundary side temperatures

  double Fourier_num = (alpha * delta_t) /
                       (pow(delta_x, 2.0));  // alpha*dt/dx^2 the Fourier number

  sim(delta_x, delta_t, imax, t0, tboundary,
      Fourier_num);  // run the simulation

  compare_schemes(delta_t);

  compare_error_schemes(delta_t, 0.4);

  // For error norm analysis
  double dt_list[] = {0.10, 0.05, 0.02, 0.01, 0.005};
  convergence_study(delta_x, imax, alpha, t0, tboundary, 0.4, dt_list, 5);

  return 0;
}
