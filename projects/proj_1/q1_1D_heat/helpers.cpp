#include "helpers.h"
#include <cstdio>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstring>

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

int last_block_index_and_time(const char* path, double* last_time_out)
{
    std::ifstream in(path);
    if (!in) {
        if (last_time_out) *last_time_out = 0.0;
        return 0;
    }

    std::string line;
    int idx = -1;
    double last_t = 0.0;

    while (std::getline(in, line)) {
        if (line.rfind("# t =", 0) == 0) {
            idx += 1;
            const char* s = line.c_str();
            const char* eq = std::strchr(s, '=');
            if (eq) {
                last_t = std::atof(eq + 1);
            }
        }
    }

    if (idx < 0) idx = 0;
    if (last_time_out) *last_time_out = last_t;
    return idx;
}

double latest_time_in_file(const char* path)
{
    FILE* f = std::fopen(path, "r");
    if (!f) return 0.0;

    char line[256];
    double last_t = 0.0;

    while (std::fgets(line, sizeof(line), f)) {
        if (std::strncmp(line, "# t =", 5) == 0) {
            double t = 0.0;
            if (std::sscanf(line, "# t = %lf", &t) == 1) {
                last_t = t;
            }
        }
    }

    std::fclose(f);
    return last_t;
}