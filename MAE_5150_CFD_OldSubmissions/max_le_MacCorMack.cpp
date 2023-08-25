//MAC CORMACK

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;


struct Hyperbolic
{
	double dx;
	double dt;
	int imax;
	double tmax; //max time
	int nmax; //max time steps
	double dtbydx; //dt/dx quantity
	vector <double> u_star; //intermediate u value
	vector <double> un; //u at nth level, for corrector method 
	vector <double> u; // future u, final solution (n+1)
	vector <double> e_vector; //E = u square /2
	vector <double> u_initial;
	vector <double> x_vector_MCM;
	double interval;
};


void InitializeU(struct Hyperbolic *burgers){

	//BASIC PROPERTIES	
	burgers->dx = 1.0;
	burgers->imax = 41;
	burgers->tmax = 2.4;
	burgers->nmax = (burgers->tmax/burgers->dt)+1;
	burgers->dtbydx = (burgers->dt)/(burgers->dx);
	burgers->interval = 0.4;

	//SETTING UP AND Es
	burgers->u_star.resize(burgers->imax+1);
	burgers->un.resize(burgers->imax+1);
	burgers->e_vector.resize(burgers->imax+1);
	burgers->u.resize(burgers->imax+1);
	burgers->u_initial.resize(burgers->imax+1);
	burgers->x_vector_MCM.resize(burgers->imax+1);


	burgers->x_vector_MCM[1] = 0;
	for (int i = 1; i<=burgers->imax;i++){
		burgers->x_vector_MCM[i+1] = burgers->x_vector_MCM[i]+1.0;
	}
	

	//Zeros out the vectors 

	for (int i = 1; i<=burgers->imax;i++){
		burgers->u[i] = 0.0;
		
	}





	//SETTING UP INITIAL CONDITIONS

	for (int i = 1; i<=burgers->imax;i++){
		if (i>=1 && i<=20){
			burgers->u[i] = 5.0;
		}

		else if (i>20 && i<=burgers->imax){
			burgers->u[i] = 0.0;
		}
	}

	//Zeros out Un 

	for (int i = 1; i<=burgers->imax;i++){
		burgers->un[i] = 0.0;
	}	

}


double FluxVariable(double u){
	u = 0.5*pow(u,2.);
	return u;
}

double MacCormack(struct Hyperbolic *burgers){


	InitializeU(burgers);
	int imax = burgers->imax;
	int nmax = burgers->nmax;

	double u_combined[100][100]={0.0};

	//Generating nmax_vector;

	double nmax_vector[7];

	double factor = burgers->interval/burgers->dt;

	for (int i = 2; i<=7;i++){
		nmax_vector[i] = (factor)*double(i)-factor;
	}	

	//Put into Initial

	for (int i = 1; i<=imax;i++){
		burgers->u_initial[i] = burgers->u[i];
	}

	for (int k = 1; k<=imax;k++){
		u_combined[k][1] = burgers->u_initial[k];
	}

	burgers->u_star[1] = burgers->u[1]-(burgers->dtbydx)*(FluxVariable(burgers->u[2])-FluxVariable(burgers->u[1]));	

	//Time loop
		
	for (int n = 1; n<=nmax;n++){	


		//Predictor step
		for (int i = 2; i<=imax-1;i++){

			burgers->u_star[i] = burgers->u[i]-(burgers->dtbydx)*0.5*(pow(burgers->u[i+1],2.)-pow(burgers->u[i],2.));
		}

		//Corrector step 

		for (int k = 2; k<=imax-1;k++){
			burgers->un[k] = 0.5*(burgers->u[k]+burgers->u_star[k])-(0.25*burgers->dtbydx*(pow(burgers->u_star[k],2.)-pow(burgers->u_star[k-1],2.)));

		}

		//Update U 
		for (int p = 2; p<=imax-1;p++){

			burgers->u[p] = burgers->un[p];
		}

		if (n == nmax_vector[2]){
			for (int i = 1; i<=imax;i++){
				u_combined[i][2] = burgers->u[i];
			}
		}

		if (n == nmax_vector[3]){
			for (int i = 1; i<=imax;i++){
				u_combined[i][3] = burgers->u[i];
			}
		}

		if (n == nmax_vector[4]){
			for (int i = 1; i<=imax;i++){
				u_combined[i][4] = burgers->u[i];
			}
		}

		if (n == nmax_vector[5]){
			for (int i = 1; i<=imax;i++){
				u_combined[i][5] = burgers->u[i];
			}
		}

		if (n == nmax_vector[6]){
			for (int i = 1; i<=imax;i++){
				u_combined[i][6] = burgers->u[i];
			}
		}

		if (n == nmax_vector[7]){
			for (int i = 1; i<=imax;i++){
				u_combined[i][7] = burgers->u[i];
			}
		}

	}		

	for (int i = 1; i<=imax;i++){
		u_combined[i][8] = burgers->x_vector_MCM[i];
	}
	


	

	if (burgers->dt == 0.1){
		FILE *outfile1;
		outfile1 = fopen("dt01burgers.dat","w");
		fprintf(outfile1,"Table shows solution at dt = 0.1 sec\n\n");
		fprintf(outfile1,"0.0 sec\t\t0.4 sec\t\t0.8 sec\t\t1.2 sec\t\t1.6 sec\t\t2.0 sec\t\t2.4 sec\t\tX POSITION\n");
		fprintf(outfile1,"------------------------------------------------------------------------------------------------\n");
		for (int i = 1; i<=41;i++){
			for (int j = 1;j<=8;j++){
				if ( j == 8){
					fprintf(outfile1,"%4.2f\n",u_combined[i][j]);
				}
				else{
				fprintf(outfile1,"%f\t", u_combined[i][j]);
				}
			}
		fprintf(outfile1,"\n");
		}
	}
	

	else if (burgers->dt == 0.2){
		FILE *outfile2;
		outfile2 = fopen("dt02burgers.dat","w");
		fprintf(outfile2,"Table shows solution at dt = 0.2 sec\n\n");
		fprintf(outfile2,"0.0 sec\t\t0.4 sec\t\t0.8 sec\t\t1.2 sec\t\t1.6 sec\t\t2.0 sec\t\t2.4 sec\t\tX POSITION\n");
		fprintf(outfile2,"------------------------------------------------------------------------------------------------\n");
		for (int i = 1; i<=41;i++){
			for (int j = 1;j<=8;j++){
				if ( j == 8){
					fprintf(outfile2,"%4.2f\n",u_combined[i][j]);
				}
				else{
				fprintf(outfile2,"%f\t", u_combined[i][j]);
				}
			}
		fprintf(outfile2,"\n");
		}

	}
}

void runMacCorMack(){
	Hyperbolic burgers;
	burgers.dt = 0.2;
	InitializeU(&burgers);
	MacCormack(&burgers);	
}


int main(){
	runMacCorMack();		
}
