#include <iostream>
#include <stdio.h>
#include <cmath>
#include <vector>
using namespace std;

//Declaring struct, vector just like the Point GS

struct Elliptic
{
	double dx ;
	double dy;
	int imax;
	int jmax;
	double ERROR_MAX;
	vector<vector <double> > psi;	
	vector<vector <double> > psi_updated;
	double ERROR; 
	int maxiter;
	int iter;
};


void thomasTriDiagonalY(int i, double a[],double b[],double c[], double d[], struct Elliptic *laplace){

	int jmax = laplace->jmax-1;
	double dprime[jmax+1];
	double cprime[jmax+1];
	dprime[1] = d[1];
	cprime[1] = c[1];

	//FORWARD LOOP
	for (int j = 2; j<=jmax;j++){
		dprime[j] = d[j] - ((b[j]*a[j-1])/(dprime[j-1]));
		cprime[j] = c[j] - ((cprime[j-1]*b[j])/(dprime[j-1]));
	}

	laplace->psi_updated[i][jmax] = cprime[jmax]/dprime[jmax];
	
	//BACKWARD LOOP
	for (int j = jmax-1;j>=2;j--){
		laplace->psi_updated[i][j] = (cprime[j]-(a[j]*laplace->psi_updated[i][j+1]))/(dprime[j]);
	}
}


int main(){

	//--------------------------------------------------------------//
	Elliptic laplace;
	laplace.dx = 0.2;
	laplace.dy = 0.2;
	laplace.imax = 6;
	laplace.jmax = 6;
	laplace.ERROR_MAX = 0.01; //max allowable error
	laplace.maxiter = 500;
	laplace.iter = 0;
	double psi1 = 0.0;
	double psi2 = 0.0;
	double psi3 = 100.00;
	double psi4 = 0.0;
	double beta = (laplace.dx)/(laplace.dy);
	double TimesBetaSquare = (2.*(1.+pow(beta,2.)));
	double betaSquare = pow(beta,2.);

	// vector<vector<double> > psi;
	laplace.psi.resize(laplace.imax+1);
	laplace.psi_updated.resize(laplace.imax+1);

	// psi.resize(imax+1);
	for (int i = 1; i<=laplace.imax;i++){
		laplace.psi[i].resize(laplace.jmax+1);
		laplace.psi_updated[i].resize(laplace.jmax+1);
	}

	//bottom BC,initially ZERO everywhere

	for (int j = 1; j<=2;j++){
		laplace.psi[1][j] = 0.0;

	}

	// but 100 from j = 6 to j= jmax
		
	for (int j = 3; j<=laplace.jmax;j++){
			laplace.psi[1][j] = 100.00;
			laplace.psi_updated[1][j] = 100.00;

	}


	for (int j = 1; j<=laplace.jmax;j++){
		laplace.psi[laplace.imax][j] = 0.0;
	}


	for (int i = 2; i<=laplace.imax-1;i++){
		laplace.psi[i][1] = 0.0;
	}


	//THOMAS PARAMETERS

	//----------------------------------------------------------------//
	//Define vectors for Thomas (in Y)

	double ay[laplace.jmax+1]={0.0}; //above
	double by[laplace.jmax+1]={0.0}; //below
	double cy[laplace.jmax+1]={0.0};//rhs
	double diagonalY[laplace.jmax+1]={0.0};//diagonal

	//Fill out values
	//For ay,by,Dy from 1 to JMAX for now
	//Rewrite value below
	for (int j=1; j<=laplace.jmax;j++){
		ay[j] = 1.;
		by[j] = 1.;
		diagonalY[j] = -1.0*TimesBetaSquare;
	}



	//Special values rewrite

	//Don't use zero elements; zero out these
	//Really dont need, because loop at 2 to jmax-1
	ay[0] = 0.0;
	by[0] = 0.0;
	cy[0] = 0.0;
	diagonalY[0] = 0.0;

	//All the 1st elements = 0
	ay[1] = 0.0;
	by[1] = 0.0;
	cy[1] = psi1; //which is 0
	diagonalY[1] = 1.0;

	//All last values = 0 

	ay[laplace.jmax] = 0.0;
	by[laplace.jmax] = 0.0;
	cy[laplace.jmax] = 0.0;
	diagonalY[laplace.jmax] = 1.0;


	//Diagonals, d has full, a misses last, b misses first
	 
	ay[laplace.jmax-1] = 0.0;// a misses last
	
	by[2] = 0.0; //b misses first


	//----------------------------------------------------------------------------//

	printf("Test printing out 1st row:\n");

	//Testing to print out 1st row

	for (int j = 1; j<=laplace.jmax;j++){
		printf("%.2f\t",laplace.psi[1][j] );
	}

	printf("\n");


	//Test printing out the thomas vectors
	printf("Test printing out thomas vectors\n\n");
	printf("i         |   A   |    D      |   B         |    C   | \n");
	printf("-----------------------------------------------------------\n");
	for (int j =1; j<=laplace.jmax;j++){
		printf("%d\t  |%5.2f  |\t%5.2f  |\t%5.2f|\t%6.2f|\n",j,ay[j],diagonalY[j],by[j],cy[j]);
	}

	printf("\n");


	printf("Test printing out original psi matrix, before any changes\n\n");

	for (int i = 1; i<=laplace.imax;i++){
		for (int j = 1; j<=laplace.jmax;j++){
			printf("%.2f\t", laplace.psi[i][j]);
		}
		printf("\n");
	}

	printf("\n");


	//--------------------------------------------------------------//



	//LINE GAUSS SEIDEL LOOP

	do{
		laplace.ERROR = 0.0;		
		// Loop through iTH row
		for (int i = 2; i<=laplace.imax-1;i++){
			//Loop through jTH row
			for (int j = 2; j<=laplace.jmax-1;j++){

				//--------------------------------------------------------------

				// //BOUNDARY CONDITIONS

				//at J = 2, zero 

				if (j == 2){
					cy[j] = (
						(-betaSquare*(laplace.psi[i+1][j]))+
						(-betaSquare*laplace.psi_updated[i-1][j])-
						(psi1)
					);

				} 


				//at J = JMAX -1, will have psi[JMAX], which needs dPsi/dx = 0
				else if (j == laplace.jmax-1){
					// laplace.psi[i][laplace.jmax] = 1*(laplace.psi[i][laplace.jmax-1]);

					// cy[j] = cy[j+1];
					cy[j] = (
						(-betaSquare*(laplace.psi[i+1][j]))+
						(-betaSquare*laplace.psi_updated[i-1][j])-
						(laplace.psi_updated[i][j+1])
					);
				}

				else{										
					cy[j] = (
						-(betaSquare*laplace.psi_updated[i+1][j])-
						(betaSquare*laplace.psi_updated[i-1][j]));
					
				}

				// //Calling thomas, should ovewrite Psi_updated
				
			}
			thomasTriDiagonalY(i,ay,by,cy,diagonalY,&laplace);		

			laplace.psi_updated[i][laplace.jmax] = 
				laplace.psi_updated[i][laplace.jmax-1]; 		

			for(int j=2; j<=laplace.jmax-1; j++) {
				laplace.ERROR += 
					abs((laplace.psi_updated[i][j] - laplace.psi[i][j]));
			}

		}

		// printf("%f\n",laplace.ERROR);

		
		//Updating Pupdated with P

		for (int i = 2; i<=laplace.imax-1;i++){
			for (int j =2; j<=laplace.jmax-1;j++){
				laplace.psi[i][j] = laplace.psi_updated[i][j];
			}
		}	

		//Update the iteration counter
		laplace.iter = laplace.iter + 1;
			
	}while(laplace.ERROR > laplace.ERROR_MAX);




	printf("\n");
	printf("Print psi, after changes\n");



	//Testing to print out results 
	printf("Max iteration is %d \n", laplace.iter);
	for (int i = 1; i<=laplace.imax;i++){
		printf("%d\t",i);
		for (int j =1; j<=laplace.jmax;j++){
			printf("%6.2f\t",laplace.psi_updated[i][j]);
		}
		printf("\n");
	}	
}

