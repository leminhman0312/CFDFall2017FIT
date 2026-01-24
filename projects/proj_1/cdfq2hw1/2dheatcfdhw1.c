#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;



double thomasTriDiagonal(int imax, double a[],double b[],double c[], double d[],double u[]){
	double dprime[imax+1];
	double cprime[imax+1];
	dprime[2] = d[2];
	cprime[2] = c[2];

	//FORWARD LOOP
	for (int i = 3; i<=imax-1;i++){

		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));		
	}

	u[imax-1] = cprime[imax-1]/dprime[imax-1];
	
	//BACKWARD LOOP
	for (int i = imax-2;i>=2;i--){
		u[i] = (cprime[i]-(a[i]*u[i+1]))/(dprime[i]);
	}

	for (int i =2;i<=imax-1;i++){
		printf("%f\t", u[i]);
	}
	printf("\n");

}






int main(){
	//Define basic parameters
	double dt = 0.01;
	int nmax = 10;
	double deltax = 0.1; 
	double deltay = deltax;
	double alpha = 0.645; //thermal diffusivity
	double maxlength_X = 3.5; //3.5 by 3.5
	double minlength_X = 0.0; //assume start at 0
	double maxlength_Y = 3.5; //
	double minlength_Y = 0.0; //assume start at 0
	double numberofpointsX = (maxlength_X-minlength_X)/(deltax);
	double numberofpointsY = (maxlength_Y-minlength_Y)/(deltay);
	int numberofpointsX_int = int(numberofpointsX);
	int numberofpointsY_int = int(numberofpointsY);
	int imax = numberofpointsX_int;
	int jmax = numberofpointsY_int;

	double t0 = 0.0; //initial temperature
	double t1 = 200.00; //boundary temperature
	double t2 = 200.00; //boundary temperature
	double t3 = 0.0;
	double t4 = 0.0;

	//Define U
	//Set Initial U at t = 0 (or n = 1) to 0 everywhere
	double u[imax+1][jmax+1]; //solution at zero time

	for (int i = 1; i<=imax;i++){
		for (int j = 1; j<=jmax;j++){
			u[i][j] = 0.0;
			u[1][j] = t1;
			u[i][1] = t2; 
		}
	}

	//U is good, 35 by 35
	

	//Define constants for tri-diagonal matrices

	double dx = (alpha*dt)/(pow(deltax,2.)); //diffusion number in X
	double dy = (alpha*dt)/(pow(deltay,2.)); //diffusion number in Y
	double d1 = 0.5*dx;
	double d2 = 0.5*dy;

	double plus2d2 = 1+(2.*d2);
	double plus2d1 = 1+(2.*d1);
	double minus2d2 = 1-(2.*d2);
	double minus2d1 = 1-(2.*d1);

	//THOMAS PARAMETERS

	//Define vectors for Thomas in X 

	double ax[imax+1]={0.0}; //above
	double bx[imax+1]={0.0}; //below
	double cx[imax+1];//rhs
	double diagonalX[imax+1]={0.0};//diagonal


	//Fill out values
	//For ax,bx
	for (int i = 2; i<=imax-1;i++){
		ax[i] = -1.*d1;
		ax[imax-1] = 0.0;
		bx[i] = -1.*d1;
		bx[2] = 0.0;
		diagonalX[i] = plus2d1;
	}



	//Create dummy matrices for X sweep 
	//Initialized to 0 for now, will be rewritten

	double u_dummy_matrix[imax+1][jmax+1] = {0.0};
	double sweep_X[imax+1]; //vectorto store individual thomas in X

	//Put values into Udummy = U
	for (int i = 1; i<=imax;i++){
		for (int j = 1; j<=jmax;j++){
			u_dummy_matrix[i][j] = u[i][j];
		}
	}

	//Define vectors for Thomas in Y

	double ay[imax+1]={0.0}; //above
	double by[imax+1]={0.0}; //below
	double cy[imax+1];//rhs
	double diagonalY[imax+1]={0.0};//diagonal

	//Fill out values
	//For ax,bx
	for (int i = 2; i<=imax-1;i++){
		ay[i] = -1*d2;
		ay[imax-1] = 0.0;
		by[i] = -1*d2;
		by[2] = 0.0;
		diagonalY[i] = plus2d2;
	}

	//Create final matrix after Y sweep = real solution at n+1
	//Initialized to 0 for now, will be rewritten

	double u_real[imax+1][jmax+1] = {0.0};
	double sweep_Y[imax+1]; //vector to store individual thomas in Y

	

	//for (int n = 1; n<=nmax;n++){

				//Except for Cx[2], and Cx[imax-1] (BC)
				//Else it is this (in case t0 is not zero)
	for (int j = 2;j<=jmax-1;j++){
		for (int i = 2;i<=imax-1;i++){
			if (i == 2){
				cx[i] = d2*u[i][j+1] + (minus2d2*u[i][j])+(d2*u[i][j-1])+(d1*t2); //BC
			} 
			else if (i==imax-1){
				cx[i] = d2*u[i][j+1] + (minus2d2*u[i][j])+(d2*u[i][j-1])+(d1*t4);
			}
			else {
				cx[i] = d2*u[i][j+1] + (minus2d2*u[i][j])+(d2*u[i][j-1]);
			}
		}
	}	

	//Checked ax,bx,cx,diagonalX = All OK


	//X sweep 


	//thomasTriDiagonal(imax,ax,bx,cx,diagonalX,sweep_X);


	for (int n = 1; n<=nmax;n++){

		for (int j = 2; j<=jmax;j++){
			//Call thomas, thomas returns a vector (row)
			sweep_X[j] = thomasTriDiagonal(imax,ax,bx,cx,diagonalX,sweep_X);
			printf("%f\n", sweep_X[j]);
			//Replace row from 2 to jmax
			//for (i = )
			std::copy(sweep_X,sweep_X+jmax,u_dummy_matrix[j]);				
		}
	}
}

