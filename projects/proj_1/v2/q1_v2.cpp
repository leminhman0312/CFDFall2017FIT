//MAE 5150 CFD: Coding Project 1 
//Max Le
//Question 1: 1D Heat Conduction Equation
// du/dt + alpha* du/dx = 0


#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <sys/stat.h>
#include <sys/types.h>
//Functions' Prototypes

//Calculations
void thomasTriDiagonal(int imax, double a[],double b[],double c[], double d[],double u[]);
void FTCS_return_implicit(int nmax, double diffusion, double tboundary, double t0, int imax,double x[],double u[]);
void FTCS_return_implicit_interior(int nmax, double diffusion, double tboundary, double t0, int imax, double u_full[]);



//Finite Difference Methods
void FTCS_explicit(double t0, double tboundary, int imax, int nmax, double dt, double diffusion, double x[],FILE* outfile1, double out_un[]);
void DufortFrankel(double t0, double tboundary , int imax, int nmax, double dt, double diffusion, double x[], FILE* outfile2, double out_un[]);
void FTCS_implicit(int nmax, double diffusion, double tboundary, double t0, int imax,double x[],FILE *outfile3, double out_un[]);
void CrankNicolson_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[],double out_un[]);

//Simulations
void sim(double delta_x,double delta_t, int imax, double t0, double tboundary, double diffusion);
void plot(double delta_t);
void compare_schemes(double delta_t, double t_target);
void compare_dt(const char* scheme, double dt1, double dt2, double t_target);
void compare_dt_all_times(const char* scheme, double dt1, double dt2);

// helpers
static void make_dt_tag(double dt, char tag[8]);
static void ensure_compare_dir();
static void ensure_plot_dir();
static void ensure_data_dir();
static bool file_exists(const char* path);
static bool scheme_data_exists(const char* scheme, double dt);













//MAIN CODE

int main(){
	
	//BASIC PARAMETERS

	int imax = 21; //max space
	double delta_x = 0.05; 
	double delta_t = 0.05;
	
	double alpha = 0.1;//thermal diffusivity
	double t0 = 100.00;//initial temperature 
	double tboundary = 300.00;//boundary temperature
	double diffusion = (alpha*delta_t)/(pow(delta_x,2.)); //diffusion number
	double x_vector[imax+1]; //spatial vector X
	//------------------------------------------------------------------------------------------------------------------
	// Fill up first value of X-vector
	x_vector[1] = 0.00;


	//Fill up X-vector for space
	
	for (int i = 1; i<=imax-1; i++){
		x_vector[i+1] = x_vector[i] + delta_x;
	}	


	// simulations
	sim(delta_x,delta_t,imax,t0,tboundary,diffusion);

	//plotting
	plot(delta_t);
	compare_schemes(0.01,0.1);

	compare_dt_all_times("ftcs_explicit", 0.01, 0.05);
	compare_dt_all_times("dufort",        0.01, 0.05);
	compare_dt_all_times("ftcs_implicit", 0.01, 0.05);
	compare_dt_all_times("cn",            0.01, 0.05);


	
	return 0;

}


static bool scheme_data_exists(const char* scheme, double dt)
{
    char tag[8];
    make_dt_tag(dt, tag);

    char fname[128];
    std::snprintf(fname, sizeof(fname), "data/%s_%s.txt", scheme, tag);
    return file_exists(fname);
}



static bool file_exists(const char* path)
{
    struct stat st;
    return (stat(path, &st) == 0);
}


static void ensure_plot_dir()
{
    struct stat st;
    if (stat("plot", &st) != 0) {
        mkdir("plot", 0755);
    }
}

static void ensure_compare_dir()
{
    struct stat st;
    if (stat("compare", &st) != 0) {
        mkdir("compare", 0755);
    }
}

void compare_dt(const char* scheme, double dt1, double dt2, double t_target)
{
    ensure_compare_dir();

    if (!scheme_data_exists(scheme, dt1) || !scheme_data_exists(scheme, dt2)) {
        return; // do nothing if either dataset missing
    }

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

    for (int k = 0; k < 4; k++){
        compare_dt(scheme, dt1, dt2, times[k]);
    }
}


static void make_dt_tag(double dt, char tag[8])
{
    // dt in hours, like 0.01 or 0.05
    // convert to integer thousandths of an hour: 0.01 -> 10 -> "001" style is your old convention
    // Your files use 001 for 0.01 and 005 for 0.05 so tag is dt*100 rounded
    int code = (int)lround(dt * 100.0);   // 0.01 -> 1, 0.05 -> 5
    std::snprintf(tag, 8, "%03d", code);  // 1 -> "001", 5 -> "005"
}

void sim(double delta_x, double delta_t, int imax,
         double t0, double tboundary, double diffusion)
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
    FILE* outfile_CN       = fopen(f_cn,       "w");

    double x_vector[imax+1];
    x_vector[1] = 0.0;
    for (int i = 1; i <= imax-1; i++) x_vector[i+1] = x_vector[i] + delta_x;

    double T_explicit[imax+1];
    double T_dufort[imax+1];
    double T_implicit[imax+1];
    double T_CN[imax+1];

    for (int k = 0; k < 5; k++){
        double t_target = 0.1 * k;
        int nmax = (int)lround(t_target / delta_t);

        FTCS_explicit(t0, tboundary, imax, nmax, delta_t, diffusion, x_vector, outfile_explicit, T_explicit);
        DufortFrankel(t0, tboundary, imax, nmax, delta_t, diffusion, x_vector, outfile_dufort, T_dufort);
        FTCS_implicit(nmax, diffusion, tboundary, t0, imax, x_vector, outfile_implicit, T_implicit);
        CrankNicolson_implicit(nmax, diffusion, tboundary, t0, imax, x_vector, T_CN);

        fprintf(outfile_explicit, "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_dufort,   "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_implicit, "# t = %.2f hr\n# x\tT\n", t_target);
        fprintf(outfile_CN,       "# t = %.2f hr\n# x\tT\n", t_target);

        for (int i = 1; i <= imax; i++){
            fprintf(outfile_explicit, "%.6f\t%.6f\n", x_vector[i], T_explicit[i]);
            fprintf(outfile_dufort,   "%.6f\t%.6f\n", x_vector[i], T_dufort[i]);
            fprintf(outfile_implicit, "%.6f\t%.6f\n", x_vector[i], T_implicit[i]);
            fprintf(outfile_CN,       "%.6f\t%.6f\n", x_vector[i], T_CN[i]);
        }

        fprintf(outfile_explicit, "\n\n");
        fprintf(outfile_dufort,   "\n\n");
        fprintf(outfile_implicit, "\n\n");
        fprintf(outfile_CN,       "\n\n");
    }

    fclose(outfile_explicit);
    fclose(outfile_dufort);
    fclose(outfile_implicit);
    fclose(outfile_CN);
}



void compare_schemes(double delta_t, double t_target)
{
	ensure_compare_dir();
    char tag[8];
    make_dt_tag(delta_t, tag);

    int idx = (int)lround(t_target / 0.1);   // because your blocks are 0.0,0.1,0.2,... and index 0 is t=0
    // so t=0.1 -> idx=1, t=0.2 -> idx=2, etc

    char cmd[256];
    std::snprintf(cmd, 256,
        "gnuplot -e \"tag='%s'; dt='%.2f'; idx=%d; tlabel='%.1f'\" gnuplot_scripts/compare_schemes.gp",
        tag, delta_t, idx, t_target);
    system(cmd);
}


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


static void ensure_data_dir()
{
    struct stat st;
    if (stat("data", &st) != 0) {
        mkdir("data", 0755);
    }
}



//END MAIN CODE

//FUNCTIONS

//FTCS EXPLICIT
void FTCS_explicit(double t0, double tboundary, int imax, int nmax, double dt, double diffusion, double x[],FILE* outfile1, double out_un[]){
	
	//set up initial conditions, T = 100 of interior nodes 
	double u[imax];
	double un[imax];
	for (int i = 2; i<=imax-1;i++){
		u[i] = t0; 
		un[i] = t0;
	}

	//set up boundary condition = 300 for outside node, 1 and imax
	u[1] = tboundary;
	u[imax] = tboundary;
	un[1] = tboundary;
	un[imax] = tboundary;

	//time marching
	for (int n = 1; n<=nmax; n++){	
		
		//copy un to u 
		for (int i = 1; i<=imax;i++){
			u[i] = un[i];
		}

		//update interior
		for (int i = 2; i<=imax-1;i++){
			un[i] = u[i] + diffusion*(u[i+1]-(2*u[i])+u[i-1]);
		}

		//enforce BC everystep
		un[1] = tboundary;
		un[imax] = tboundary;
	}
	
	
	//copy un to output array

	for (int i = 1; i <= imax; i++){
		out_un[i] = un[i];
	}

}


//DUFORT FRANKEL IMPLICIT 1st step
void DufortFrankel(double t0, double tboundary , int imax, int nmax, double dt, double diffusion, double x[], FILE* outfile2, double out_un[]){
	
    double d = 2.*diffusion;
	//set up initial conditions, T=t0  of interior nodes 
	double un_m1[imax+1] ={0.0};
	double un_p1[imax+1] = {0.0};
	double un[imax+1];

    un_m1[1] = tboundary;
    un_m1[imax] = tboundary;

	for (int i = 2; i<=imax-1;i++){
		un_m1[i] = t0;
	}	


	//USE FTCS EXPLICIT TO SOLVE FOR 2ND TIME STEP
    double u1[imax+1];
    // FTCS_return_implicit(1,diffusion,tboundary,t0,imax,x,u1);
	FTCS_return_implicit_interior(1,diffusion,tboundary,t0,imax,u1);
    for (int i = 1; i <= imax; i++) un[i] = u1[i]; //update

    //Finite Difference Loop
	//TIME LOOP
	for (int n = 1; n<=nmax-1;n++){

		//FDE using the intial FTCS implict as starters
		for (int i = 2; i<=imax-1;i++){
            un_p1[i] = ((1.0-d)*un_m1[i] + d*(un[i+1]+un[i-1]))/(1.0+d);
        }

        // enforce BC everystep
        un_p1[1] = tboundary;
        un_p1[imax] = tboundary;
        
        // update the time levels
        for (int i = 1; i <= imax; i++){
            un_m1[i] = un[i];
            un[i] = un_p1[i];
        }
	}

    // copy un to output array
    for (int i = 1; i <= imax; i++){
        out_un[i] = un[i];
    }
}


void FTCS_return_implicit_interior(int nmax, double diffusion, double tboundary, double t0, int imax, double u_full[]){
	
	int N = imax -2; //only deal with interior node so less size
	
	//Define arrays
    double a[N+1];
    double b[N+1];
    double c[N+1];
	double d[N+1];
	double u_int[N+1];

	// zeroes out
	for (int i = 0; i < N+1; i++) {
		a[i] = 0.0;
		b[i] = 0.0;
		c[i] = 0.0;
		d[i] = 0.0;
        u_int[i] = 0.0;
    }

	//initialize full field at t = 0
	u_full[1] = tboundary;
	u_full[imax] = tboundary;
	for (int i = 2; i <= imax-1; i++) {
		u_full[i] = t0;
	}

	// tridiagonal coefficients
	// first points
	b[1] = 0.0;
	d[1] = 1.0 + 2.0*diffusion;
	a[1] = -diffusion;

	// inner points
	for (int i = 1; i<=N; i++) {
		d[i] = 1.0 + 2.0*diffusion;
		a[i] = -diffusion;
		b[i] = -diffusion;
	}

	// end points
	b[N] = -diffusion;
	d[N] = 1.0 + 2.0*diffusion;
	a[N] = 0.0;


	

	for (int t = 1; t <= nmax; t++) {
		
		//RHS from previous time level PLUS boundary conditions
		//map j=1 => node i = 2 
		//map j=N => node i = imax-1
		for (int j = 1; j <= N; j++) {
			int i = j + 1;
			c[j] = u_full[i];
		}

		// add boundary conditions to RHS
		c[1] += diffusion*tboundary; // + d*U1
		c[N] += diffusion*tboundary; // + d*Uimax

		//solve for interior at new time
		thomasTriDiagonal(N,a,b,c,d,u_int);

		//write back to full array and reapply boundary conditions
		u_full[1] = tboundary;
		u_full[imax] = tboundary;
		for (int j = 1; j <= N; j++) {
			int i = j+1;
			u_full[i] = u_int[j];
		}
	}
}


void FTCS_return_implicit(int nmax, double diffusion, double tboundary, double t0, int imax,double x[], double u[]){

	//Define arrays
    double a[imax+1];
    double b[imax+1];
    double c[imax+1];
	double d[imax+1];


    // initial conditions
    u[1] = tboundary;
    u[imax] = tboundary;
    for (int i = 2; i <= imax-1 ; i++) u[i] = t0;


    // end points for thomas
    d[1] = 1.0;
    a[1] = 0.0;
    b[1] = 0.0;

    d[imax] = 1.0;
    a[imax] = 0.0;
    b[imax] = 0.0;


    // interior thomas coefficients
	for (int i = 2;i<=imax-1; i++){
		d[i] = (1+(2*diffusion));
		a[i] = -diffusion;
		b[i] = -diffusion;
	}


	//Call thomas and time loop
	for (int t = 1; t<=nmax;t++){
        c[1] = tboundary;
        c[imax] = tboundary;
        for (int i = 2; i <= imax - 1; i++){
            c[i] = u[i];
        }
        thomasTriDiagonal(imax,a ,b ,c ,d ,u );
	}

}



void thomasTriDiagonal(int imax, double a[],double b[],double c[], double d[],double u[]){
			
	double dprime[imax+1];
	double cprime[imax+1];
	dprime[1] = d[1];
	cprime[1] = c[1];

	//FORWARD LOOP

	for (int i = 2; i<=imax;i++){
		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));		
	}
    
    u[imax] = cprime[imax] / dprime[imax];
    
	//BACKWARD LOOP
	for (int i = imax-1;i>=1;i--){
		u[i] = (cprime[i]-(a[i]*u[i+1]))/(dprime[i]);
	}	
}


void FTCS_implicit(int nmax, double diffusion, double tboundary, double t0, int imax,double x[],FILE *outfile3, double out_un[]){
	

	//Define arrays
    double a[imax+1];
    double b[imax+1];
    double c[imax+1];
	double d[imax+1];
	double un[imax+1];

	// initial conditions
    un[1] = tboundary;
    un[imax] = tboundary;
    for (int i = 2; i <= imax-1 ; i++) un[i] = t0;


    // end points for thomas
    d[1] = 1.0;
    a[1] = 0.0;
    b[1] = 0.0;

    d[imax] = 1.0;
    a[imax] = 0.0;
    b[imax] = 0.0;

	 // interior thomas coefficients
	for (int i = 2;i<=imax-1; i++){
		d[i] = (1+(2*diffusion));
		a[i] = -diffusion;
		b[i] = -diffusion;
	}

  
	//Call thomas and time loop
	for (int t = 1; t<=nmax;t++){
        c[1] = tboundary;
        c[imax] = tboundary;
        for (int i = 2; i <= imax - 1; i++){
            c[i] = un[i];
        }
        thomasTriDiagonal(imax,a ,b ,c ,d ,un );
	}



	//copy un to output array

	for (int i = 1; i <= imax; i++){
		out_un[i] = un[i];
	}
}

void CrankNicolson_implicit(int nmax, double diffusion, double tboundary, double t0, int imax,double x[],double out_un[]){


	double d_NC = diffusion/2.;
	

	//Define arrays
    double a[imax+1];
    double b[imax+1];
    double c[imax+1];
	double d[imax+1];
	double un[imax+1];
	double u0[imax+1]; // to store full step size


	//set up boundary condition = 300 for outside node, 1 and imax
	u0[1] = tboundary;
	u0[imax] = tboundary;

	// initial conditions
	for (int i = 2; i<=imax-1;i++){
		u0[i] = t0;
	}

	double u_half[imax+1]; //to store half step size
	u_half[1] = tboundary;
	u_half[imax] = tboundary;

	// end points for thomas
	d[1] = 1.0;
	a[1] = 0.0;
	b[1] = 0.0;

	d[imax] = 1.0;	
	a[imax] = 0.0;
	b[imax] = 0.0;

	// interior thomas coefficients
	for (int i = 2;i<=imax-1; i++){
		d[i] = (1+2.0*d_NC);
		a[i] = -d_NC;
		b[i] = -d_NC;
	}
		
	//TIME LOOP
	for (int n = 1; n<=nmax;n++){

		// STEP 1: Explicit half step
		//enforce boundary at every point
		u_half[1] = tboundary;
		u_half[imax] = tboundary;

		//interior points
		for (int i = 2; i<=imax-1;i++){
			u_half[i] = u0[i] + d_NC*(u0[i+1]-(2*u0[i])+u0[i-1]);
		}

		// STEP 2: Implicit half step
		// update the rhs vector c
		// boundary for c
		c[1] = tboundary;
		c[imax] = tboundary;

		//interior points for c
		for (int i = 2; i <= imax-1; i++){
			c[i] = u_half[i];
		}
		
		thomasTriDiagonal(imax,a,b,c,d,u0);

		u0[1] = tboundary;
		u0[imax] = tboundary;
	}	

	//copy un to output array

	for (int i = 1; i <= imax; i++){
		out_un[i] = u0[i];
	}

}
//END FUNCTIONS
