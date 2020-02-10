#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;



struct Elliptic
{
	double dx ;
	double dy;
	double imax;
	double jmax;
	double ERROR_MAX;
	vector<vector <double> > psi;	
	vector<vector <double> > psi_updated;
	double ERROR; 
	int maxiter;
	int iter;
};

int main(){
	//basic data

	Elliptic laplace;
	laplace.dx = 0.2;
	laplace.dy = 0.2;
	laplace.imax = 21;
	laplace.jmax = 31;
	laplace.ERROR_MAX = 0.01; //max allowable error
	laplace.maxiter = 500;
	laplace.iter = 0;
	double w = 1.9;


	// vector<vector<double> > psi;
	laplace.psi.resize(laplace.imax+1);
	laplace.psi_updated.resize(laplace.imax+1);

	// psi.resize(imax+1);
	for (int i = 1; i<=laplace.imax;i++){
		laplace.psi[i].resize(laplace.jmax+1);
		laplace.psi_updated[i].resize(laplace.jmax+1);
	}

	//bottom BC,initially ZERO everywhere
		
		for (int j = 7; j<=laplace.jmax;j++){
			laplace.psi[1][j] = 100.00;
			laplace.psi_updated[1][j] = 100.00;
		}
	
	
		
	//POINT GAUSS-SEIDEL

	double beta = (laplace.dx)/(laplace.dy);
	double Omegaby2BetaSquare = (w/(2*(1+pow(beta,2.))));
	double betaSquare = pow(beta,2.);

	// //Iteration loop 
	
	do{
		laplace.ERROR = 0.0;			
		for (int i = 2; i<=laplace.imax-1;i++){
			for (int j = 2; j<=laplace.jmax-1;j++){	
					//Finite Difference, using psi_updated as Latest Data

					laplace.psi_updated[i][j] = ((1-w)*(laplace.psi[i][j]))+(Omegaby2BetaSquare)*(laplace.psi[i+1][j]+laplace.psi_updated[i-1][j]+(betaSquare)*(laplace.psi[i][j+1]+laplace.psi_updated[i][j-1]));					

					//Calculate error, keep doing this until satisfy ERROR MAX
					laplace.ERROR += abs((laplace.psi_updated[i][j] - laplace.psi[i][j]));					
			}
		}

		//Updating Pupdated with P

		for (int i = 2; i<=laplace.imax-1;i++){
			for (int j =2; j<=laplace.jmax-1;j++){
				laplace.psi[i][j] = laplace.psi_updated[i][j];
			}
		}


		// Make sure BC satisfied, dpsi/dx = 0
		for (int i = 1; i<=laplace.imax;i++){
				laplace.psi[i][laplace.jmax] = laplace.psi[i][laplace.jmax-1];
		}

		//Update the iteration counter
		laplace.iter = laplace.iter + 1;

		}while(laplace.ERROR > laplace.ERROR_MAX);
		



		printf("----------------Final result------------------------\n");
		// //Testing to print out results 

		FILE* outfile;
		outfile = fopen("testing.dat","w");



		for (int i = 1; i<=laplace.imax;i++){
			for (int j =1; j<=laplace.jmax;j++){
				fprintf(outfile,"%6.2f\t", laplace.psi[i][j]);
			}
			fprintf(outfile,"\n");
		}


		printf("ERROR: %f\n", laplace.ERROR);

	return 0;
}
		










	

