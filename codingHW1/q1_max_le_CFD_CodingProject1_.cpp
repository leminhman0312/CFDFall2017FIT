//MAE 5150 CFD: Coding Project 1 
//Max Le
//Question 1: 1D Heat Conduction Equation

#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cstdlib>

//Functions' Prototypes

//Calculations
double thomasTriDiagonal(int imax, double a[],double b[],double c[], double d[],double u[],double x[]);
double * FTCS_return_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[]);

//Finite Difference Methods
double FTCS_explicit(double t0, double tboundary, int imax, int nmax, double dt, double dx, double x[],FILE* outfile1);
double DufortFrankelEXPLICIT(double t0, double tboundary , int imax, int nmax, double dt, double dx, double x[], FILE* outfile2);
double FTCS_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[],FILE *outfile3);
double CrankNicolson_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[],FILE *outfile4);

//MAIN CODE

int main(){
	
	//BASIC PARAMETERS

	int imax = 21; //max space
	//int nmax= 2;//max time, 1 = 1 time step, 2 = 2 time steps....
	double deltax = 0.05; //delta x
	double dt = 0.05;
	
	double alpha = 0.1;//thermal diffusivity
	double t0 = 100.00;//initial 
	double tboundary = 300.00;//boundary
	double dx = (alpha*dt)/(pow(deltax,2.)); //diffusion number
	double x_vector[imax+1]; //spatial vector X
	

	//------------------------------------------------------------------------------------------------------------------
	//set up initial condition for X vector (0 at not i=1)
	x_vector[1] = 0.00;

	// CrankNicolson_implicit(nmax_vector[i],dx,tboundary,t0,imax,x_vector,outfile4);

	

	//Fill up X-vector for space
	
	for (int i = 1; i<=imax; i++){
		x_vector[i+1] = x_vector[i] + deltax;
	}	

	// double ftcs_explicit_value[5][imax];
	// double DufortFrankel_value[5][imax];
	double ftcs_implicit_value[5][imax];
	// double CrankNicolson_value[5][imax];


	//FILES to output.  Change dt value above then comment out the appropriate sections


	//--------------------------------------------------------------------------------------------------------------------

	//FILES FOR dt = 0.01

	// FILE *outfile1;
	// outfile1 = fopen("FTCS_explicit001.txt","w");
	// FILE *outfile2;
	// outfile2 = fopen("DufortFrankel001.txt","w");
	FILE *outfile3;
	outfile3 = fopen("FTCS_implicit001.txt","w");
	// FILE* outfile4;
	// outfile4 = fopen("CrankNicolson001.txt","w");
	int nmax_vector[5] = {0,10,20,30,40}; //This is to store the step size at 0.1, 0.2,0.3 and 0.4 hr

	for (int i = 1; i<=4;i++){
	// ftcs_explicit_value[i][imax] = FTCS_explicit(t0,tboundary,imax,nmax_vector[i],dt,dx,x_vector,outfile1);
	// DufortFrankel_value[i][imax] = DufortFrankelEXPLICIT(t0,tboundary,imax,nmax_vector[i],dt,dx,x_vector,outfile2);
	ftcs_implicit_value[i][imax] = FTCS_implicit(nmax_vector[i],dx,tboundary,t0,imax,x_vector,outfile3);
	// CrankNicolson_value[i][imax] = CrankNicolson_implicit(nmax_vector[i],dx,tboundary,t0,imax,x_vector,outfile4);
	}


	//-----------------------------------------------------------------------------------------------------------------------

	
	// //FILES FOR dt = 0.05

	// FILE *outfile1;
	// outfile1 = fopen("FTCS_explicit005.txt","w");
	// FILE *outfile2;
	// outfile2 = fopen("DufortFrankel005.txt","w");
	// FILE *outfile3;
	// outfile3 = fopen("FTCS_implicit005.txt","w");
	// FILE* outfile4;
	// outfile4 = fopen("CrankNicolson005.txt","w");

	// int nmax_vector[5] = {0,2,4,6,8}; //This is to store the step size at 0.1, 0.2,0.3 and 0.4 hr


	// for (int i = 1; i<=4;i++){
	// ftcs_explicit_value[i][imax] = FTCS_explicit(t0,tboundary,imax,nmax_vector[i],dt,dx,x_vector,outfile1);
	// DufortFrankel_value[i][imax] = DufortFrankelEXPLICIT(t0,tboundary,imax,nmax_vector[i],dt,dx,x_vector,outfile2);
	// // ftcs_implicit_value[i][imax] = FTCS_implicit(nmax_vector[i],dx,tboundary,t0,imax,x_vector,outfile3);
	// CrankNicolson_value[i][imax] = CrankNicolson_implicit(nmax_vector[i],dx,tboundary,t0,imax,x_vector,outfile4);
	// }


	return 0;

}


//END MAIN CODE


//-------------------------------------------------------------------------------------------------------------------------


//FUNCTIONS

//FTCS EXPLICIT
double FTCS_explicit(double t0, double tboundary, int imax, int nmax, double dt, double dx, double x[],FILE* outfile1){
	
	//set up initial conditions = 100 of interior nodes 
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
	//Print X
	
	for (int i = 1; i<=imax;i++){
		fprintf(outfile1,"%f\t",x[i]);
	}
	fprintf(outfile1,"\n");

	//Finite Difference Loop
	//TIME LOOP
	for (int n = 1; n<=nmax; n++){	
		
		//Copy loop 

		for (int i = 1; i<=imax;i++){
			u[i] = un[i];
		}
		//FDE looping

		for (int i = 2; i<=imax-1;i++){
			// u0[i] = u[i] + dx*(u[i+1]-(2*u[i])+u[i-1]);
			un[i] = u[i] + dx*(u[i+1]-(2*u[i])+u[i-1]);
		}
	}	
	//Printing to file
	for (int k = 1; k<=imax;k++){
			fprintf(outfile1,"%f\t",un[k]);
			printf("%f\n",un[k]);
	}

	fprintf(outfile1,"\n");
	printf("\n");	
}


double thomasTriDiagonal(int imax, double a[],double b[],double c[], double d[],double u[],double x[]){
			
	double dprime[imax];
	double cprime[imax];
	dprime[1] = d[1];
	cprime[1] = c[1];
	//FORWARD LOOP

	for (int i = 2; i<=imax;i++){
		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));		
	}	
	//BACKWARD LOOP
	for (int i = imax-1;i>=2;i--){
		u[i] = (cprime[i]-(a[i]*u[i+1]))/(dprime[i]);
	}	
}


double FTCS_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[],FILE *outfile3){
	//Print X
	for (int i = 1; i<=imax;i++){
		fprintf(outfile3,"%f\t",x[i]);
	}
	fprintf(outfile3,"\n");
	
	double minus1plus2dx = -(1+(2*dx));

	//Define solution vector for FTCS implicit
	double u_FTCS_implicit[imax];
	u_FTCS_implicit[1] = tboundary; //boundary
	u_FTCS_implicit[imax] = tboundary; //boundary

	//Define special vectors for Thomas Algorithm
	// B = below, A = above, D = main diagonal
	double b_FTCS_implicit[imax]={0};
	double a_FTCS_implicit[imax]={0};
	double d_FTCS_implicit[imax]={0};

	for (int i = 2;i<=imax-1; i++){
		d_FTCS_implicit[i] = minus1plus2dx;
		a_FTCS_implicit[i] = dx;
		b_FTCS_implicit[i] = dx;
	}

	d_FTCS_implicit[1] = -1.0;
	d_FTCS_implicit[imax] = -1.0;

	//Set up special conditions for B, and C

	a_FTCS_implicit[imax-1] = dx;
	a_FTCS_implicit[1] = 0.0;

	b_FTCS_implicit[2] = dx;
	b_FTCS_implicit[imax] = 0.0;


	//Define right handed vectors

	double c_FTCS_implicit[imax]={t0};
	c_FTCS_implicit[1] = -1*tboundary;
	c_FTCS_implicit[imax] = -1*tboundary;
	for (int k = 2; k<=imax-1; k++){
		c_FTCS_implicit[k] = -1*t0;
	}
	//Call thomas and time loop
	for (int t =1; t<=nmax;t++){
		for (int i =1; i<=imax;i++){
			thomasTriDiagonal(imax,a_FTCS_implicit,b_FTCS_implicit,c_FTCS_implicit,d_FTCS_implicit,u_FTCS_implicit,x);
			c_FTCS_implicit[i] = -1*u_FTCS_implicit[i];
		}
	}
	//Printing
	for (int k =1; k<=imax;k++){
		fprintf(outfile3,"%f\t",u_FTCS_implicit[k]);
		printf("%f\n", u_FTCS_implicit[k]);
	}
	fprintf(outfile3,"\n");
}

double CrankNicolson_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[],FILE *outfile4){
	
	double u[imax];
	for (int i = 2; i<=imax-1;i++){
		u[i] = t0;
	}

	//set up boundary condition = 300 for outside node, 1 and imax
	u[1] = tboundary;
	u[imax] = tboundary;




	//PART 1: EXPLICIT

	double u_dumb[imax]={0.0}; //to store half step size
	u_dumb[1] = tboundary;
	u_dumb[imax] = tboundary;
	//Finite Difference Loop
	//TIME LOOP
	for (int n = 1; n<=nmax;n++){
		//assign values at previous time
		//use u0 as placeholder
		for (int j = 2; j<=imax-1; j++){
			u_dumb[j] = u[j];
		}
		//FDE looping
		for (int i = 2; i<=imax-1;i++){
			u[i] = u_dumb[i] + (dx/4.)*(u[i+1]-(2*u[i])+u[i-1]);
		}		
	}	

	//-----------------------------------------------------------------------//

	//Now Uhalf is updated = U0
	
	//PART 2: IMPLICIT 
	//Define diagonal vectors to input to Thomas Crank Nicolson
	//same as before in FTCS implicit, but dx=dx/4, -(1+dx) becomes 
	// (1+ dx/2)

	//Define solution vector for Crank Nicolson implicit
	double plusdxby2 = (1+((dx)/(2.)));
	double u_NC_implicit[imax];
	u_NC_implicit[1] = tboundary; //boundary
	u_NC_implicit[imax] = tboundary; //boundary

	//Define special vectors for Thomas Algorithm
	// B = below, A = above, D = main diagonal
	double b_NC_implicit[imax]={0};
	double a_NC_implicit[imax]={0};
	double d_NC_implicit[imax]={0};

	for (int i = 2;i<=imax-1; i++){
		d_NC_implicit[i] = -1.*plusdxby2;
		a_NC_implicit[i] = (dx)/(4.);
		b_NC_implicit[i] = (dx)/(4.);
	}
	
	d_NC_implicit[1] = -1.0;
	d_NC_implicit[imax] = -1.0;

	//Set up special conditions for B, and C

	a_NC_implicit[imax-1] = (dx)/(4.);
	a_NC_implicit[1] = 0.0;

	b_NC_implicit[2] = (dx)/(4.);
	b_NC_implicit[imax] = 0.0;
	

	//Define right handed vectors, basically negative
	// of Udumb

	double c_NC_implicit[imax];
	for (int i = 1;i<=imax;i++){
		c_NC_implicit[i] = -1.*u_dumb[i];
	}


	//Feed A,B,C,D into Thomas and Time loop

	for (int t = 1; t<=nmax;t++){
		for (int i = 1;i<=imax;i++){
			thomasTriDiagonal(imax,a_NC_implicit,b_NC_implicit,c_NC_implicit,d_NC_implicit,u_NC_implicit,x);
			u_NC_implicit[i]=u_dumb[i];
			c_NC_implicit[i] = -1.*u_dumb[i];
		}
	}

	//Print out final solution for NC explicit

	for (int i = 1; i<=imax;i++){
		fprintf(outfile4,"%f\t",u_NC_implicit[i]);
		printf("%f\t%f\n",x[i],u_NC_implicit[i]);
	}
	fprintf(outfile4,"\n");
	printf("\n");
}

double * FTCS_return_implicit(int nmax, double dx, double tboundary, double t0, int imax,double x[]){
	double minus1plus2dx = -(1+(2*dx));

	//Define solution vector for FTCS implicit
	static double u_FTCS_implicit_return[21];
	u_FTCS_implicit_return[1] = tboundary; //boundary
	u_FTCS_implicit_return[imax] = tboundary; //boundary

	//Define special vectors for Thomas Algorithm
	// B = below, A = above, D = main diagonal
	double b_FTCS_implicit_return[imax]={0};
	double a_FTCS_implicit_return[imax]={0};
	double d_FTCS_implicit_return[imax]={0};

	for (int i = 2;i<=imax-1; i++){
		d_FTCS_implicit_return[i] = minus1plus2dx;
		a_FTCS_implicit_return[i] = dx;
		b_FTCS_implicit_return[i] = dx;
	}

	d_FTCS_implicit_return[1] = -1.0;
	d_FTCS_implicit_return[imax] = -1.0;

	//Set up special conditions for B, and C

	a_FTCS_implicit_return[imax-1] = dx;
	a_FTCS_implicit_return[1] = 0.0;

	b_FTCS_implicit_return[2] = dx;
	b_FTCS_implicit_return[imax] = 0.0;

	//Define right handed vectors

	double c_FTCS_implicit_return[imax]={t0};
	c_FTCS_implicit_return[1] = -1*tboundary;
	c_FTCS_implicit_return[imax] = -1*tboundary;
	for (int k = 2; k<=imax-1; k++){
		c_FTCS_implicit_return[k] = -1*t0;
	}

	//Call thomas and time loop
	for (int t =1; t<=nmax;t++){
		for (int i =1; i<=imax;i++){
			thomasTriDiagonal(imax,a_FTCS_implicit_return,b_FTCS_implicit_return,c_FTCS_implicit_return,d_FTCS_implicit_return,u_FTCS_implicit_return,x);
			c_FTCS_implicit_return[i] = -1*u_FTCS_implicit_return[i];
		}
	}

	return u_FTCS_implicit_return; //return the pointer array to feed into DufortFrankel

}


//DUFORT FRANKEL IMPLICIT 1st step
double DufortFrankelEXPLICIT(double t0, double tboundary , int imax, int nmax, double dt, double dx, double x[], FILE* outfile2){
	
	//set up initial conditions = 100 of interior nodes 
	//DECLARE DUMMY U FOR AFTER FIRST STEP

	double temp[imax] ={0.0};
	double u_plus_one[imax] = {0.0};
	// double u_previous[imax] = {0.0};
	double u[imax];
	for (int i = 2; i<=imax-1;i++){
		u[i] = t0;
		u_plus_one[i] = t0;
		temp[i] = t0;	
		// u_previous[i]= t0;
	}

	//set up boundary condition = 300 for outside node, 1 and imax	
	
	u[1] = tboundary;
	u[imax] = tboundary;
	u_plus_one[imax] = tboundary;
	u_plus_one[1] = tboundary;	
	temp[imax] = tboundary;
	temp[1] = tboundary;
	// Print X

	for (int i = 1; i<=imax;i++){
		fprintf(outfile2,"%f\t",x[i]);
	}
	fprintf(outfile2,"\n");


	//Let un be results from previous time step
	double un[imax] = {0.0};

	//USE FTCS EXPLICIT TO SOLVE FOR 2ND TIME STEP

	//Declare pointer to get FTCS implicit to calculate at 1 time step
	double * u_previous;

	//Finite Difference Loop
	//TIME LOOP
	for (int n = 1; n<=nmax;n++){
		//Copy loop

		

		//Calculate at 1 time step
		if (n == 1){
			for (int i = 2; i<=imax-1;i++){
				u_previous = FTCS_return_implicit(2,dx,tboundary,t0,imax,x);
			}
		}

		for (int c = 1; c<=imax;c++){
			u[c] = u_previous[c];
		}
		//FDE using the intial FTCS implict as starters
			for (int j = 2; j<=imax-1;j++){
				u_previous[j] = (((1-2.*dx)/(1+2.*dx))*u[j]) + ((2.*dx)/(1+2.*dx))*(u[j+1]+u[j-1]);
			}			
		
	}


	//Testing to print out U_plus_one

	for (int k = 1; k<=imax;k++){
		printf("%f\n",u[k]);
		fprintf(outfile2, "%f\t",u[k]);
	}
	printf("\n");
	fprintf(outfile2,"\n");
}

//END FUNCTIONS
