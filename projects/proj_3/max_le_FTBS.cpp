//FTBS EXPLICIT 

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
	wave->nmax = (wave->tmax/wave->dt);
	wave->c = (wave->alpha*wave->dt)/(wave->dx);

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


void FTBSExplicit(struct Hyperbolic *wave){	

	//Finite Difference Loop 

	//Time loop, FULL TIME
	for (int n = 0; n<=wave->nmax;n++){
		//copy loop 
		for (int i = 1; i<=wave->imax;i++){
			wave->u[i] = wave->u0[i];			
		}

		//FDE
		for (int i = 2; i<=wave->imax;i++){
			wave->u0[i] = wave->u[i] - (wave->c)*(wave->u[i]-wave->u[i-1]);			
		}

		if (n == 0){
			for (int i = 1; i<=wave->imax;i++){
				wave->uInitial[i] = wave->u[i];
			}

		}

		else if (n == (wave->nmax/2)){
			for (int i = 1; i<=wave->imax;i++){
				wave->uhalf[i] = wave->u[i]; 
			}
		}
		else if (n == wave->nmax){
			for (int i = 1; i<=wave->imax;i++){
				wave->ufull[i] = wave->u[i];
			}
		}		
	}

	if (wave->dt == 0.005){
		FILE* outfile1;
		outfile1 = fopen("005FTBS.dat","w");
		fprintf(outfile1,"Initial Condition       Halfway Point        Final Condition\n");
		fprintf(outfile1,"---------------------------------------------------------------\n");
		for (int i = 1; i<=wave->imax;i++){
			// printf("%f\t\t%f\t\t%f\n",wave.uInitial[i],wave.uhalf[i],wave.ufull[i]);
			fprintf(outfile1,"%10.7f\t\t\t\t%10.7f\t\t\t\t%10.7f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
		}
	}


	else if (wave->dt == 0.0025){
		FILE* outfile2;
		outfile2 = fopen("0025FTBS.dat","w");
		fprintf(outfile2,"Initial Condition       Halfway Point        Final Condition\n");
		fprintf(outfile2,"---------------------------------------------------------------\n");
		for (int i = 1; i<=wave->imax;i++){
			// printf("%f\t\t%f\t\t%f\n",wave.uInitial[i],wave.uhalf[i],wave.ufull[i]);
			fprintf(outfile2,"%10.7f\t\t\t\t%10.7f\t\t\t\t%10.7f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
		}
	}

	else if (wave->dt == 0.00125){
		FILE* outfile3;
		outfile3 = fopen("0125FTBS.dat","w");
		fprintf(outfile3,"Initial Condition       Halfway Point        Final Condition\n");
		fprintf(outfile3,"---------------------------------------------------------------\n");
		for (int i = 1; i<=wave->imax;i++){
			// printf("%f\t\t%f\t\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			// fprintf(outfile3,"%10.7f\t\t\t\t%10.7f\t\t\t\t%10.7f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);
			fprintf(outfile3,"%f\t%f\t%f\n",wave->uInitial[i],wave->uhalf[i],wave->ufull[i]);


		}
	}
}


void runFTBSCode(){
	Hyperbolic wave;
	wave.dt = 0.00125;
	InitializeU(&wave);
	FTBSExplicit(&wave);
}


int main(){
	runFTBSCode();

}
