//BTCS IMPLICIT 

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

//DECLARE STRUCTURE

struct Hyperbolic
{
	double dx;
	double dt;
	int imax;
	double tmax; //max time
	double nmax; //max time steps
	double alpha; //speed of sound
	double c;
	double c_square;
	vector <double> u0; //initial u //USED IN LOOPS
	vector <double> u; //u at n+1 level, USED IN LOOPS
	vector <double> uhalf; //half way point, for print
	vector <double> ufull; //end solution, for print
	vector <double> uInitial; //initial sol, for print 
 
};

void InitializeU(struct Hyperbolic *wave){

	//BASIC PROPERTIES	
	wave->dx = 1.0;
	wave->imax = 71;
	wave->alpha = 200.00;
	wave->tmax = 0.15;
	wave->nmax = (wave->tmax/wave->dt)+1.;

	wave->c = (wave->alpha*wave->dt)/(wave->dx);
	wave->c_square = pow(wave->c,2.);

	//SETTING U AND UN
	wave->u0.resize(wave->imax+1);
	wave->u.resize(wave->imax+1);
	wave->uhalf.resize(wave->imax+1);
	wave->ufull.resize(wave->imax+1);
	wave->uInitial.resize(wave->imax+1);
	
	

	//SETTING UP INITAL CONDITIONS, U0
	for (int i = 1; i<=wave->imax;i++){
		
		if (i == 15){
			wave->u0[i]= 20.00;
		}

		else if (i>=5 && i<15){
			wave->u0[i] = (2*i)-10;
		}
		else if (i>15 && i<=25){
			wave->u0[i] = (-2*i) + 50;
		}
	}
}


void thomasTriDiagonal(int nmax, double a[],double b[],double c[], double d[],struct Hyperbolic *wave){


	int imax = wave->imax-1;
	double dprime[imax+1];
	double cprime[imax+1];
	dprime[1] = d[1];
	cprime[1] = c[1];

	//FORWARD LOOP
	for (int i = 2; i<=imax;i++){
		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));
	}



	wave->u0[imax] = cprime[imax]/dprime[imax];

	// u[imax] = cprime[imax]/dprime[imax];


	//BACKWARD LOOP
	for (int i = imax-1;i>=2;i--){
		wave->u0[i] = (cprime[i]-(a[i]*wave->u0[i+1]))/(dprime[i]);

		// u[i] = (cprime[i]-(a[i]*u[i+1]))/(dprime[i]);
	}
}

void BTCS_implicit(struct Hyperbolic *wave){



	//THOMAS + TRIDIAGONAL PARAMETERS

	int imax = wave->imax;
	int nmax = wave->nmax;

	double halfC = 0.5*wave->c;

	double a[imax+1] = {0.0}; //above
	double b[imax+1] = {0.0}; //below
	double c[imax+1]; //rhs
	double d[imax+1] = {0.0};//diagonal

	for (int i=2; i<=imax-1;i++){
		a[i] = -1.*halfC;
		b[i] = halfC;
		d[i] = -1.;
	}

	//Special values

	
	a[imax-1] = 0.0;
	a[1] = 0.0;

	b[1] = 0.0;
	b[imax] = 0.0;

	d[1] = 1.;
	d[imax] = 1.;


	c[1] = wave->u[1];
	c[imax] = wave->u[imax];

	for (int i = 1; i<=imax-1;i++){
		c[i] = -1.*wave->u0[i];
	}




	//Create a U initial 

	for (int i = 1; i<=imax;i++){
		wave->uInitial[i] = wave->u0[i];
	}

	// thomasTriDiagonal(nmax,a,b,c,d,wave);

	//Finite Difference Loop 

	for (int n = 0; n<=nmax;n++){
		thomasTriDiagonal(nmax,a,b,c,d,wave);
		for (int i = 1;i<=imax;i++){
			c[i] = -1.*wave->u0[i];

			
			if (n == nmax/2){
				wave->uhalf[i] = wave->u0[i];
			}

			else if (n == nmax){
				wave->ufull[i] = wave->u0[i];
			}				
		}
	}

	if (wave->dt == 0.005){
		FILE* outfile1;
		outfile1 = fopen("05implicit.dat","w");
		fprintf(outfile1,"Initial Condition       Halfway Point        Final Condition\n");
		fprintf(outfile1,"---------------------------------------------------------------\n");
		for (int i = 1; i<=wave->imax;i++){
			// printf("%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			fprintf(outfile1,"%10.7f\t\t\t\t%10.7f\t\t\t\t%10.7f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			// fprintf(outfile1,"%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);

		}
	}

	else if (wave->dt == 0.0025){
		FILE* outfile2;
		outfile2 = fopen("025implicit.dat","w");
		fprintf(outfile2,"Initial Condition       Halfway Point        Final Condition\n");
		fprintf(outfile2,"---------------------------------------------------------------\n");
		for (int i = 1; i<=wave->imax;i++){
			// printf("%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			fprintf(outfile2,"%10.7f\t\t\t\t%10.7f\t\t\t\t%10.7f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			// fprintf(outfile2,"%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);

		}
	}
	else if (wave->dt == 0.00125){
		FILE* outfile3;
		outfile3 = fopen("0125implicit.dat","w");
		fprintf(outfile3,"Initial Condition       Halfway Point        Final Condition\n");
		fprintf(outfile3,"---------------------------------------------------------------\n");	
		for (int i = 1; i<=wave->imax;i++){
			// printf("%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			fprintf(outfile3,"%10.7f\t\t\t\t%10.7f\t\t\t\t%10.7f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			// fprintf(outfile3,"%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);

		}
	}
}
		
void runBTCSCode(){
	Hyperbolic wave;
	wave.dt = 0.00125;
	InitializeU(&wave);
	BTCS_implicit(&wave);
}

int main(){
	runBTCSCode();	
}
