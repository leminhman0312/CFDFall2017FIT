#include "plotting.h"

#include <cstdio>
#include <cstring>

void plotContourMatlabLike(const char* datafile, const char* outpng,
                           double time_hr, const char* scheme) {
  char cmd[1024];
  std::snprintf(cmd, sizeof(cmd),
                "gnuplot -e \""
                "datafile=\\\"%s\\\"; "
                "outpng=\\\"%s\\\"; "
                "tlabel=%.3f; "
                "scheme=\\\"%s\\\"\" "
                "gnuplot_scripts/plot_contour_2d_matlab_like.gp",
                datafile, outpng, time_hr, scheme);
  printf("GNUPLOT CMD:\n%s\n", cmd);
  system(cmd);
}