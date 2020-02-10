#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <ctime>
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
	double w ; //relaxation parameter
	double beta;
	double oneby2BetaSquare;
	double betaSquare;
	double psi1;
	double psi3;
	double OneMinusOmega;
	double OmegaBetaSquare;
	double TimesBetaSquare;
	double Omegaby2BetaSquare;
	float timeElapsed;
};


void InitializePsi(struct Elliptic *laplace, float wi){
	laplace->dx = 0.2;
	laplace->dy = 0.2;
	laplace->imax = 21;
	laplace->jmax = 31;
	laplace->ERROR_MAX = 0.01; //max allowable error
	laplace->ERROR=0.0;
	laplace->maxiter = 500;
	laplace->iter = 0;
	laplace->psi.resize(laplace->imax+1);
	laplace->psi_updated.resize(laplace->imax+1);
	laplace->w = wi;
	laplace->beta = (laplace->dx)/(laplace->dy);
	laplace->oneby2BetaSquare = (1./(2*(1+pow(laplace->beta,2.))));
	laplace->betaSquare = pow(laplace->beta,2.);
	laplace->psi1 = 0.0;
	laplace->psi3 = 100.00;
	laplace->OneMinusOmega = (1-laplace->w);
	laplace->OmegaBetaSquare = laplace->w*laplace->betaSquare;
	laplace->TimesBetaSquare = 2*(1+pow(laplace->beta,2.));
	laplace->Omegaby2BetaSquare = (laplace->w/(2*(1+pow(laplace->beta,2.))));
		
	for (int i = 1; i<=laplace->imax;i++){
		laplace->psi[i].resize(laplace->jmax+1);
		laplace->psi_updated[i].resize(laplace->jmax+1);
	}

	//SET zero to everywhere
	for (int i = 1; i<=laplace->imax;i++){
		for (int j = 1; j<=laplace->jmax;j++){
			laplace->psi[i][j] = 0.0;
			laplace->psi_updated[i][j] = 0.0;		}
	}
	//bottom BC,initially ZERO everywhere			
	for (int j = 7; j<=laplace->jmax;j++){
		laplace->psi[1][j] = laplace->psi3;
		laplace->psi_updated[1][j] = laplace->psi3;
	}
}





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


void PointGaussSeidel(struct Elliptic *laplace){

	do{
		laplace->ERROR = 0.0;			
		for (int i = 2; i<=laplace->imax-1;i++){
			for (int j = 2; j<=laplace->jmax-1;j++){	
					//Finite Difference, using psi_updated as Latest Data

					// laplace->psi_updated[i][j] = (laplace->oneby2BetaSquare)*(laplace->psi[i+1][j]+laplace->psi_updated[i-1][j]+(laplace->betaSquare)*(laplace->psi[i][j+1]+laplace->psi_updated[i][j-1]));

					laplace->psi_updated[i][j] = (laplace->oneby2BetaSquare)*(laplace->psi[i][j+1]+laplace->psi_updated[i][j-1]+(laplace->betaSquare)*(laplace->psi[i+1][j]+laplace->psi_updated[i-1][j]));						

					//Calculate error, keep doing this until satisfy ERROR MAX
					laplace->ERROR += abs((laplace->psi_updated[i][j] - laplace->psi[i][j]));					
			}
		}

		//Updating Pupdated with P

		for (int i = 2; i<=laplace->imax-1;i++){
			for (int j =2; j<=laplace->jmax-1;j++){
				laplace->psi[i][j] = laplace->psi_updated[i][j];
			}
		}


		// Make sure BC satisfied, dpsi/dx = 0
		for (int i = 2; i<=laplace->imax-1;i++){
				laplace->psi[i][laplace->jmax] = laplace->psi[i][laplace->jmax-1];
		}

		//Update the iteration counter
		laplace->iter = laplace->iter + 1;

		// printf("iter =  %d\n\n", laplace->iter);
		//Testing to print out results 
		for (int i = 1; i<=laplace->imax;i++){
			for (int j =1; j<=laplace->jmax;j++){
				// printf("%6.2f\t", laplace->psi[i][j]);
			}
			// printf("\n");
		}


		
		}while(laplace->ERROR > laplace->ERROR_MAX);
		// printf("\n");

		printf("PGS Converged! Max iter is: %d \n", laplace->iter);

		

		//PRINTING TO FILE

		FILE * outfile1; 
		outfile1 = fopen("ResultsPointGS.dat","w");


		// printf("----------------Final Results----------------\n");

		// printf("Number of iteration: %d\n\n", laplace->iter);
		//Testing to print out results 
		for (int i = 1; i<=laplace->imax;i++){
			for (int j =1; j<=laplace->jmax;j++){
				// printf("%6.2f\t", laplace->psi[i][j]);
				fprintf(outfile1,"%6.5f\t", laplace->psi[i][j]);
			}
			// printf("\n");
			fprintf(outfile1,"\n");
		}
}


void LineGaussSeidel(struct Elliptic *laplace){

	//THOMAS PARAMETERS
	//----------------------------------------------------------------//
	//Define vectors for Thomas (in Y)

	int jmax = laplace->jmax;

	double ay[jmax+1]={0.0}; //above
	double by[jmax+1]={0.0}; //below
	double cy[jmax+1]={0.0};//rhs
	double diagonalY[jmax+1]={0.0};//diagonal

	//Fill out values
	//For ay,by,Dy from 1 to JMAX for now
	//Rewrite value below
	for (int j=1; j<=laplace->jmax;j++){
		ay[j] = 1.;
		by[j] = 1.;
		diagonalY[j] = -1.0*laplace->TimesBetaSquare;
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
	cy[1] = laplace->psi1; //which is 0
	diagonalY[1] = 1.0;

	//All last values = 0 

	ay[jmax] = 0.0;
	by[jmax] = 0.0;
	cy[jmax] = 0.0;
	diagonalY[jmax] = 1.0;


	//Diagonals, d has full, a misses last, b misses first
	 
	ay[jmax-1] = 0.0;// a misses last
	
	by[2] = 0.0; //b misses first


	//LINE GAUSS SEIDEL LOOP

	do{
		laplace->ERROR = 0.0;		
		// Loop through iTH row
		for (int i = 2; i<=laplace->imax-1;i++){
			//Loop through jTH row
			for (int j = 2; j<=laplace->jmax-1;j++){

				// //BOUNDARY CONDITIONS

				//at J = 2, zero 

				if (j == 2){
					cy[j] = (
						(-laplace->betaSquare*(laplace->psi[i+1][j]))+
						(-laplace->betaSquare*laplace->psi_updated[i-1][j])-
						(laplace->psi1)
					);
				}
				//at J = JMAX -1, will have psi[JMAX], which needs dPsi/dx = 0
				else if (j == laplace->jmax-1){					
					cy[j] = (
						(-laplace->betaSquare*(laplace->psi[i+1][j]))+
						(-laplace->betaSquare*laplace->psi_updated[i-1][j])-
						(laplace->psi_updated[i][j+1])
					);
				}

				else{										
					cy[j] = (
						-(laplace->betaSquare*laplace->psi_updated[i+1][j])-
						(laplace->betaSquare*laplace->psi_updated[i-1][j]));					
				}

								
			}
			thomasTriDiagonalY(i,ay,by,cy,diagonalY,laplace);	

			laplace->psi_updated[i][laplace->jmax] = 
				laplace->psi_updated[i][laplace->jmax-1]; 		

			for(int j=2; j<=laplace->jmax-1; j++) {
				laplace->ERROR += 
					abs((laplace->psi_updated[i][j] - laplace->psi[i][j]));
			}
		}

		//Updating Pupdated with P

		for (int i = 1; i<=laplace->imax;i++){
			for (int j =1; j<=laplace->jmax;j++){
				laplace->psi[i][j] = laplace->psi_updated[i][j];
			}
		}	

		//Update the iteration counter
		laplace->iter = laplace->iter + 1;
		
	}while(laplace->ERROR > laplace->ERROR_MAX);

	printf("Line GS Converged! Max iter is: %d, w = %f\n", laplace->iter,laplace->w);
	

	//PRINTING TO FILE

	FILE * outfile2; 
	outfile2 = fopen("ResultsLineGS.dat","w");

	for (int i = 1; i<=laplace->imax;i++){
		for (int j =1; j<=laplace->jmax;j++){
			// printf("%6.2f\t", laplace->psi[i][j]);
			fprintf(outfile2,"%6.2f\t", laplace->psi[i][j]);
		}
		fprintf(outfile2,"\n");
		// printf("\n");
	}
}







void LineSOR(struct Elliptic *laplace){

	//start clock
	clock_t t;

	t = clock();

	int jmax = laplace->jmax;
	//THOMAS PARAMETERS
	//Define vectors for Thomas (in Y)

	double ay[jmax+1]={0.0}; //above
	double by[jmax+1]={0.0}; //below
	double cy[jmax+1]={0.0};//rhs
	double diagonalY[jmax+1]={0.0};//diagonal

	//Fill out values
	//For ay,by,Dy from 1 to JMAX for now
	//Rewrite value below
	for (int j=1; j<=laplace->jmax;j++){
		ay[j] = laplace->w;
		by[j] = laplace->w;
		diagonalY[j] = -1.0*laplace->TimesBetaSquare;
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
	cy[1] = laplace->psi1; //which is 0
	diagonalY[1] = 1.0;

	//All last values = 0 

	ay[jmax] = 0.0;
	by[jmax] = 0.0;
	cy[jmax] = 0.0;
	diagonalY[jmax] = 1.0;


	//Diagonals, d has full, a misses last, b misses first
	 
	ay[jmax-1] = 0.0;// a misses last
	
	by[2] = 0.0; //b misses first


	//LINE GAUSS SEIDEL LOOP

	do{
		laplace->ERROR = 0.0;		
		// Loop through iTH row
		for (int i = 2; i<=laplace->imax-1;i++){
			//Loop through jTH row
			for (int j = 2; j<=laplace->jmax-1;j++){

				//--------------------------------------------------------------

				// //BOUNDARY CONDITIONS

				//at J = 2, zero 

				if (j == 2){
					cy[j] = -(laplace->OneMinusOmega*laplace->TimesBetaSquare*laplace->psi[i][j]) 
						-((laplace->OmegaBetaSquare)*(((laplace->psi[i+1][j]))+
						(laplace->psi_updated[i-1][j])))-
						((laplace->w)*(laplace->psi1))
					;
				}

				//at J = JMAX -1, will have psi[JMAX], which needs dPsi/dx = 0
				else if (j == laplace->jmax-1){
					cy[j] = -(laplace->OneMinusOmega*laplace->TimesBetaSquare*laplace->psi[i][j]) 
						-((laplace->OmegaBetaSquare)*(((laplace->psi[i+1][j]))+
						(laplace->psi_updated[i-1][j])))-
						((laplace->w)*(laplace->psi_updated[i][j+1]));
				}

				else{										
					cy[j] = 
						-(laplace->OneMinusOmega*laplace->TimesBetaSquare*laplace->psi[i][j]) 
						-((laplace->OmegaBetaSquare)*(((laplace->psi[i+1][j]))+
						(laplace->psi_updated[i-1][j])));
					
				}				
			}
			thomasTriDiagonalY(i,ay,by,cy,diagonalY,laplace);		

			//BC for dPsi/dx = 0;
			laplace->psi_updated[i][laplace->jmax] = 
				laplace->psi_updated[i][laplace->jmax-1]; 		

			for(int j=1; j<=laplace->jmax; j++) {
				laplace->ERROR += 
					abs((laplace->psi_updated[i][j] - laplace->psi[i][j]));
			}

		}

		//Updating Pupdated with P

		for (int i = 1; i<=laplace->imax;i++){
			for (int j =1; j<=laplace->jmax;j++){
				laplace->psi[i][j] = laplace->psi_updated[i][j];
			}
		}	

		//Update the iteration counter
		laplace->iter = laplace->iter + 1;
				
	}while(laplace->ERROR > laplace->ERROR_MAX);

	t = clock()-t;

	laplace->timeElapsed = (float)t/(CLOCKS_PER_SEC);
	
	//End clock


	printf("Line SOR Converged! Max iter is: %d, w = %f at %f seconds \n", laplace->iter,laplace->w, laplace->timeElapsed);



	// FILE *outfile3; 

	// outfile3 = fopen("ResultsLSOR.dat","w");
	// for (int i = 1; i<=laplace->imax;i++){
	// 	for (int j = 1; j<=laplace->jmax;j++){
	// 		fprintf(outfile3,"%6.5f\t", laplace->psi[i][j]);
	// 	}
	// fprintf(outfile3,"\n");
	// }
	//print to Test different w 

	FILE *outfile4; 

	outfile4 = fopen("omegaLSOR.dat","w");
	for (int i = 1; i<=laplace->imax;i++){
		for (int j = 1; j<=laplace->jmax;j++){
			fprintf(outfile4,"%6.5f\t", laplace->psi[i][j]);
		}
	fprintf(outfile4,"\n");
	}
}


void PointSOR(struct Elliptic *laplace){

	//start clock
	clock_t t;
	t = clock();
	do{
		laplace->ERROR = 0.0;			
		for (int i = 2; i<=laplace->imax-1;i++){
			for (int j = 2; j<=laplace->jmax-1;j++){	
					//Finite Difference, using psi_updated as Latest Data

					laplace->psi_updated[i][j] = ((1-laplace->w)*(laplace->psi[i][j]))+((laplace->Omegaby2BetaSquare)*(laplace->psi[i][j+1]+laplace->psi_updated[i][j-1]+(laplace->betaSquare)*(laplace->psi[i+1][j]+laplace->psi_updated[i-1][j])));	

					 //Calculate error, keep doing this until satisfy ERROR MAX
					laplace->ERROR += abs((laplace->psi_updated[i][j] - laplace->psi[i][j]));					
			}
		}

		//Updating Pupdated with P

		for (int i = 2; i<=laplace->imax-1;i++){
			for (int j =2; j<=laplace->jmax-1;j++){
				laplace->psi[i][j] = laplace->psi_updated[i][j];
			}
		}


		// Make sure BC satisfied, dpsi/dx = 0
		for (int i = 2; i<=laplace->imax-1;i++){
				laplace->psi[i][laplace->jmax] = laplace->psi[i][laplace->jmax-1];
		}

		//Update the iteration counter
		laplace->iter = laplace->iter + 1;
			
		}while(laplace->ERROR >= laplace->ERROR_MAX);

		t = clock()-t;

		laplace->timeElapsed = (float)t/(CLOCKS_PER_SEC);


		// FILE *outfile5; 

		// outfile5 = fopen("ResultsPSOR.dat","w");
		// for (int i = 1; i<=laplace->imax;i++){
		// 	for (int j = 1; j<=laplace->jmax;j++){
		// 		fprintf(outfile5,"%6.5f\t", laplace->psi[i][j]);
		// 	}
		// 	fprintf(outfile5,"\n");
		// }


		printf("POINT SOR Converged! Max iter is: %d, w = %f at %f seconds \n", laplace->iter,laplace->w, laplace->timeElapsed);


		//File to test different w 

		FILE *outfile6; 

		outfile6 = fopen("omegaPSOR.dat","w");
		for (int i = 1; i<=laplace->imax;i++){
			for (int j = 1; j<=laplace->jmax;j++){
				fprintf(outfile6,"%6.5f\t", laplace->psi[i][j]);
			}
			fprintf(outfile6,"\n");
		}




}


