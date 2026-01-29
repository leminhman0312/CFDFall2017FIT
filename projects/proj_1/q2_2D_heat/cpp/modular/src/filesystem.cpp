#include "filesystem.h"

#include <sys/stat.h>
#include <sys/types.h>

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
