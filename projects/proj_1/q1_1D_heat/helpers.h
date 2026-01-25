#ifndef HELPERS_H
#define HELPERS_H

#include <sys/stat.h>
#include <sys/types.h>

void ensure_plot_dir();
void ensure_compare_dir();
void ensure_data_dir();

void make_dt_tag(double dt, char tag[8]);

bool file_exists(const char* path);
bool scheme_data_exists(const char* scheme, double dt);

#endif
