#ifndef HELPERS_H
#define HELPERS_H

#include "headers.h"

void ensure_plot_dir();
void ensure_compare_dir();
void ensure_data_dir();

void make_dt_tag(double dt, char tag[8]);

bool file_exists(const char* path);
bool scheme_data_exists(const char* scheme, double dt);

int last_block_index_and_time(const char* path, double* last_time_out);

double latest_time_in_file(const char* path);

double err_L2(int imax, const double a[], const double b[]);
double err_Linf(int imax, const double a[], const double b[]);

#endif
