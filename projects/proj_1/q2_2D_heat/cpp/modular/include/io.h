#pragma once
#include <headers.h>

void write_field_xyz(const char* filename,
                     const std::vector<std::vector<double>>& u, double deltax,
                     double deltay);
