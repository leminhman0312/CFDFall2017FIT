#include "io.h"

#include <cstdio>

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