#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include <cmath>
#include <fstream>




//SOLVING THE 1D LINEAR DIFFUSION EQUATION

// dU/dt = Nu * d^2(U)/dx^2


int main(){

	
	///////////////////////////////////////////////////
	////////////////Start of C++ /////////////////////

	//Declare basic variables

	FILE* outfile;
	outfile = fopen("data.txt","w");

	int xmax = 2;
	int xmin = 0;
	int nx = 42; //real nx = 41, to start at 1, add 1 
	double dx = double((xmax-xmin)) /(nx-2);


	int nt = 21;
	double nu = 0.3; //viscosity
	double sigma = 0.2; // a parameter
	double dt = pow(dx,2.0)*sigma/nu;

	//cout << dt << endl; 
	
	
	double x[nx]; //space vectors
	x[1] = xmin;
	double u0[nx]; //U at time = 0
	double u[nx]; //Real U[time][space], u[i][j]
	

	//Initial Condition.  Step function

	for (int j = 1; j <nx; j++){
		
		x[j+1] = x[j] + dx;
		if (x[j]<0.5 || x[j] > 1.0){
			u0[j] = 1.0;
		}
		else{
			u0[j] = 2.0;
		}
	}

	




	//Computations

	for (int n = 0; n<=nt; n++){
		//assigning value at previous time
		//using u0 as a holder
		for (int j = 1; j<nx;j++){
			u[j] = u0[j];
		}
		for (int i =2; i< nx-1; i++){
			u0[i]=u[i]+nu*dt/(pow(dx,2.0))*(u[i+1]-2*u[i]+u[i-1]);
		}
	}

	//testing

	for (int k = 1; k<nx; k++){
		printf("%f\n",u0[k]);
	}


	//Print to file

	for (int i = 1; i<nx; i++){

			fprintf(outfile,"%f\t%f\n", x[i],u0[i]);
	}	


	



}


















