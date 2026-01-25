#include "driver.h"

#include <cstdio>
#include <cmath>
#include <cstdlib>

#include "helpers.h"
#include "methods.h"


// compare dt single time
void compare_dt(const char* scheme, double dt1, double dt2, double t_target)
{
    ensure_compare_dir();

    if (!scheme_data_exists(scheme, dt1) || !scheme_data_exists(scheme, dt2)) return;

    char tag1[8], tag2[8];
    make_dt_tag(dt1, tag1);
    make_dt_tag(dt2, tag2);

    int idx = (int)lround(t_target / 0.1);

    char cmd[512];
    std::snprintf(cmd, sizeof(cmd),
        "gnuplot -e \"scheme='%s'; tag1='%s'; dt1='%.2f'; tag2='%s'; dt2='%.2f'; idx=%d; tlabel='%.1f'\" gnuplot_scripts/compare_dt.gp",
        scheme, tag1, dt1, tag2, dt2, idx, t_target);

    system(cmd);
}

void compare_dt_all_times(const char* scheme, double dt1, double dt2)
{
    double times[4] = {0.1, 0.2, 0.3, 0.4};
    for (int k = 0; k < 4; k++) compare_dt(scheme, dt1, dt2, times[k]);
}

// compare schemes at one time
void compare_schemes(double delta_t, double t_target)
{
    ensure_compare_dir();

    char tag[8];
    make_dt_tag(delta_t, tag);

    int idx = (int)lround(t_target / 0.1);

    char cmd[256];
    std::snprintf(cmd, 256,
        "gnuplot -e \"tag='%s'; dt='%.2f'; idx=%d; tlabel='%.1f'\" gnuplot_scripts/compare_schemes.gp",
        tag, delta_t, idx, t_target);

    system(cmd);
}

// run one dt and write all scheme data files
void sim(double delta_x,double delta_t, int imax, double t0, double tboundary, double F_num)
{
    ensure_data_dir();

    char tag[8];
    make_dt_tag(delta_t, tag);

    char f_explicit[64], f_dufort[64], f_implicit[64], f_cn[64];
    std::snprintf(f_explicit, 64, "data/ftcs_explicit_%s.txt", tag);
    std::snprintf(f_dufort,   64, "data/dufort_%s.txt",        tag);
    std::snprintf(f_implicit, 64, "data/ftcs_implicit_%s.txt", tag);
    std::snprintf(f_cn,       64, "data/cn_%s.txt",            tag);

    FILE* outfile_explicit = fopen(f_explicit, "w");
    FILE* outfile_dufort   = fopen(f_dufort,   "w");
    FILE* outfile_implicit = fopen(f_implicit, "w");
    FILE* outfile_cn       = fopen(f_cn,       "w");

    double x_vector[imax+1];
    x_vector[1] = 0.0;
    for (int i = 1; i <= imax-1; i++) x_vector[i+1] = x_vector[i] + delta_x;

    double T_explicit[imax+1], T_dufort[imax+1], T_implicit[imax+1], T_cn[imax+1];

    for (int k = 0; k < 5; k++){
        double t_target = 0.1 * k;
        int nmax = (int)lround(t_target / delta_t);

        FTCS_explicit(nmax, F_num, tboundary, t0, imax, T_explicit);
        DufortFrankel(nmax, F_num, tboundary, t0, imax, T_dufort);
        FTCS_implicit(nmax, F_num, tboundary, t0, imax, T_implicit);
        CrankNicolson(nmax, F_num, tboundary, t0, imax, T_cn);

        fprintf(outfile_explicit, "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_dufort,   "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_implicit, "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_cn,       "# t = %.2f hr\n# x\tT\n", t_target);

        for (int i = 1; i <= imax; i++){
            fprintf(outfile_explicit, "%.6f\t%.6f\n", x_vector[i], T_explicit[i]);
            fprintf(outfile_dufort,   "%.6f\t%.6f\n", x_vector[i], T_dufort[i]);
            fprintf(outfile_implicit, "%.6f\t%.6f\n", x_vector[i], T_implicit[i]);
            fprintf(outfile_cn,       "%.6f\t%.6f\n", x_vector[i], T_cn[i]);
        }

        fprintf(outfile_explicit, "\n\n");
        fprintf(outfile_dufort,   "\n\n");
        fprintf(outfile_implicit, "\n\n");
        fprintf(outfile_cn,       "\n\n");
    }

    fclose(outfile_explicit);
    fclose(outfile_dufort);
    fclose(outfile_implicit);
    fclose(outfile_cn);
}

// plot all schemes for one dt
void plot(double delta_t)
{
    ensure_plot_dir();

    char tag[8];
    make_dt_tag(delta_t, tag);

    char cmd[256];

    std::snprintf(cmd, 256,
        "gnuplot -e \"scheme='ftcs_explicit'; tag='%s'; dt='%.2f'; outpng='plot/ftcs_explicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
        tag, delta_t, tag);
    system(cmd);

    std::snprintf(cmd, 256,
        "gnuplot -e \"scheme='dufort'; tag='%s'; dt='%.2f'; outpng='plot/dufort_%s.png'\" gnuplot_scripts/plot_scheme.gp",
        tag, delta_t, tag);
    system(cmd);

    std::snprintf(cmd, 256,
        "gnuplot -e \"scheme='ftcs_implicit'; tag='%s'; dt='%.2f'; outpng='plot/ftcs_implicit_%s.png'\" gnuplot_scripts/plot_scheme.gp",
        tag, delta_t, tag);
    system(cmd);

    std::snprintf(cmd, 256,
        "gnuplot -e \"scheme='cn'; tag='%s'; dt='%.2f'; outpng='plot/cn_%s.png'\" gnuplot_scripts/plot_scheme.gp",
        tag, delta_t, tag);
    system(cmd);
}
