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

    char exact_path[128];
    std::snprintf(exact_path, sizeof(exact_path), "data/exact_%s.txt", tag);

    int idx = 0;
    double tlabel = 0.0;

    if (t_target < 0.0) {
        idx = last_block_index_and_time(exact_path, &tlabel);
    } else {
        tlabel = t_target;
        idx = (int)lround(t_target / 0.1);
    }

    char cmd[256];
    std::snprintf(cmd, sizeof(cmd),
        "gnuplot -e \"tag='%s'; dt='%.2f'; idx=%d; tlabel='%.1f'\" gnuplot_scripts/compare_schemes.gp",
        tag, delta_t, idx, tlabel);

    system(cmd);
}



// run one dt and write all scheme data files
void sim(double delta_x,double delta_t, int imax, double t0, double tboundary, double F_num){
    ensure_data_dir();

    char tag[8];
    make_dt_tag(delta_t, tag);    

    char f_explicit[64], f_dufort[64], f_implicit[64], f_cn[64],f_exact[64];
    std::snprintf(f_explicit, 64, "data/ftcs_explicit_%s.txt", tag);
    std::snprintf(f_dufort,   64, "data/dufort_%s.txt",        tag);
    std::snprintf(f_implicit, 64, "data/ftcs_implicit_%s.txt", tag);
    std::snprintf(f_cn,       64, "data/cn_%s.txt",            tag);
    std::snprintf(f_exact,    64, "data/exact_%s.txt",         tag);


    FILE* outfile_explicit = fopen(f_explicit, "w");
    FILE* outfile_dufort   = fopen(f_dufort,   "w");
    FILE* outfile_implicit = fopen(f_implicit, "w");
    FILE* outfile_cn       = fopen(f_cn,       "w");
    FILE* outfile_exact    = fopen(f_exact,    "w");


    char f_err_exp[64], f_err_imp[64], f_err_df[64], f_err_cn[64];
    std::snprintf(f_err_exp, 64, "data/error_ftcs_explicit_%s.txt", tag);
    std::snprintf(f_err_imp, 64, "data/error_ftcs_implicit_%s.txt", tag);
    std::snprintf(f_err_df,  64, "data/error_dufort_%s.txt",        tag);
    std::snprintf(f_err_cn,  64, "data/error_cn_%s.txt",            tag);

    FILE* outfile_err_exp = fopen(f_err_exp, "w");
    FILE* outfile_err_imp = fopen(f_err_imp, "w");
    FILE* outfile_err_df  = fopen(f_err_df,  "w");
    FILE* outfile_err_cn  = fopen(f_err_cn,  "w");


    double x_vector[imax+1];
    x_vector[1] = 0.0;
    for (int i = 1; i <= imax-1; i++) x_vector[i+1] = x_vector[i] + delta_x;

    double L = (imax - 1) * delta_x;
    double alpha = F_num * (delta_x * delta_x) / delta_t;
    double T_exact[imax+1];

    double T_explicit[imax+1], T_dufort[imax+1], T_implicit[imax+1], T_cn[imax+1];

    for (int k = 0; k < 5; k++){
        double t_target = 0.1 * k;
        int nmax = (int)lround(t_target / delta_t);

        // compute the schemes
        FTCS_explicit(nmax, F_num, tboundary, t0, imax, T_explicit);
        DufortFrankel(nmax, F_num, tboundary, t0, imax, T_dufort);
        FTCS_implicit(nmax, F_num, tboundary, t0, imax, T_implicit);
        CrankNicolson(nmax, F_num, tboundary, t0, imax, T_cn);

        // compute the exact solution
        exact_solution_profile(imax, x_vector, t_target, alpha, L, tboundary, t0, 200, T_exact);

        // write to files
        fprintf(outfile_explicit, "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_dufort,   "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_implicit, "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_cn,       "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_exact,    "# t = %.2f hr\n# x\tT\n", t_target);

        fprintf(outfile_err_exp, "# t = %.2f hr\n# x\terr\n", t_target);
        fprintf(outfile_err_imp, "# t = %.2f hr\n# x\terr\n", t_target);
        fprintf(outfile_err_df,  "# t = %.2f hr\n# x\terr\n", t_target);
        fprintf(outfile_err_cn,  "# t = %.2f hr\n# x\terr\n", t_target);



        for (int i = 1; i <= imax; i++){
            fprintf(outfile_explicit, "%.6f\t%.6f\n", x_vector[i], T_explicit[i]);
            fprintf(outfile_dufort,   "%.6f\t%.6f\n", x_vector[i], T_dufort[i]);
            fprintf(outfile_implicit, "%.6f\t%.6f\n", x_vector[i], T_implicit[i]);
            fprintf(outfile_cn,       "%.6f\t%.6f\n", x_vector[i], T_cn[i]);
            fprintf(outfile_exact, "%.6f\t%.6f\n", x_vector[i], T_exact[i]);

            fprintf(outfile_err_exp, "%.6f\t%.6f\n", x_vector[i], T_explicit[i] - T_exact[i]);
            fprintf(outfile_err_imp, "%.6f\t%.6f\n", x_vector[i], T_implicit[i] - T_exact[i]);
            fprintf(outfile_err_df,  "%.6f\t%.6f\n", x_vector[i], T_dufort[i]   - T_exact[i]);
            fprintf(outfile_err_cn,  "%.6f\t%.6f\n", x_vector[i], T_cn[i]       - T_exact[i]);

        }

        fprintf(outfile_explicit, "\n\n");
        fprintf(outfile_dufort,   "\n\n");
        fprintf(outfile_implicit, "\n\n");
        fprintf(outfile_cn,       "\n\n");
        fprintf(outfile_exact,    "\n\n");

        fprintf(outfile_err_exp, "\n\n");
        fprintf(outfile_err_imp, "\n\n");
        fprintf(outfile_err_df,  "\n\n");
        fprintf(outfile_err_cn,  "\n\n");

    }

    fclose(outfile_explicit);
    fclose(outfile_dufort);
    fclose(outfile_implicit);
    fclose(outfile_cn);
    fclose(outfile_exact);

    fclose(outfile_err_exp);
    fclose(outfile_err_imp);
    fclose(outfile_err_df);
    fclose(outfile_err_cn);

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

void compare_error_schemes(double delta_t)
{
    char tag[8];
    make_dt_tag(delta_t, tag);

    char exact_path[128];
    std::snprintf(exact_path, sizeof(exact_path), "data/exact_%s.txt", tag);

    double t_latest = latest_time_in_file(exact_path);
    if (t_latest <= 0.0) return;

    compare_error_schemes(delta_t, t_latest);
}

void compare_error_schemes(double delta_t, double t_target)
{
    ensure_compare_dir();

    char tag[8];
    make_dt_tag(delta_t, tag);

    int idx = (int)lround(t_target / 0.1);

    char cmd[512];
    std::snprintf(cmd, sizeof(cmd),
    "gnuplot -e \"tag='%s'; dt=%.2f; idx=%d; tlabel='%.1f'; outpng='compare/error_schemes_%s_t%.1f.png'\" gnuplot_scripts/compare_error_schemes.gp",
    tag, delta_t, idx, t_target, tag, t_target);
    system(cmd);


}