#include "helpers.h"
#include <cstdio>
#include <cmath>

void ensure_plot_dir()
{
    struct stat st;
    if (stat("plot", &st) != 0) mkdir("plot", 0755);
}

void ensure_compare_dir()
{
    struct stat st;
    if (stat("compare", &st) != 0) mkdir("compare", 0755);
}

void ensure_data_dir()
{
    struct stat st;
    if (stat("data", &st) != 0) mkdir("data", 0755);
}

void make_dt_tag(double dt, char tag[8])
{
    int code = (int)lround(dt * 100.0);
    std::snprintf(tag, 8, "%03d", code);
}

bool file_exists(const char* path)
{
    struct stat st;
    return (stat(path, &st) == 0);
}

bool scheme_data_exists(const char* scheme, double dt)
{
    char tag[8];
    make_dt_tag(dt, tag);

    char fname[128];
    std::snprintf(fname, sizeof(fname), "data/%s_%s.txt", scheme, tag);
    return file_exists(fname);
}
