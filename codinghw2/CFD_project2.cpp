#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "CFD_project2.h"
#include <ctime>
using namespace std;

int main(){


	int nw = 21;
	double dw = 0.1;
	double w0 = 0.1;
	// float wi;

	float wi = 1.3;

	// Elliptic PGS;
	// InitializePsi(&PGS,wi);
	// PointGaussSeidel(&PGS);

	// Elliptic LGS;
	// InitializePsi(&LGS,wi);
	// LineGaussSeidel(&LGS);


	// PSOR best at 1.78
	// Elliptic PSOR;
	// InitializePsi(&PSOR,wi);
	// PointSOR(&PSOR);

	Elliptic LSOR;
	InitializePsi(&LSOR,wi);
	LineSOR(&LSOR);


	// // LSOR best at 1.3
	// {	float wi = 1.3;
	// 	Elliptic LSOR;
	// InitializePsi(&LSOR,wi);
	// LineSOR(&LSOR);}


	


	//TABLE FOR BLOWING UP 


	// {

	// 	int nw = 20;
	// 	double dw = 0.1;
	// 	double w0 = 0.1;
	// 	float wi;

	// 	FILE *testPSOR;
	// 	testPSOR = fopen("testPSOR.dat","w");


	// 	// fprintf(testPSOR,"Different values of W for POINT SOR\n\n");
	// 	// fprintf(testPSOR,"  w     |  Iteration  | Time[sec] |\n");
	// 	// fprintf(testPSOR,"-----------------------------------\n");
	// 	for(int i=0; i<nw; i++) {
	// 		Elliptic PSOR;
	// 		wi = w0 + dw*( (float) i);
	// 		InitializePsi(&PSOR,wi);
	// 		PointSOR(&PSOR);
	// 		// fprintf(testPSOR, "%.3f\t|\t%4.0d\t  |\t%f  |\n",PSOR.w,PSOR.iter,PSOR.timeElapsed);
	// 		fprintf(testPSOR, "%.3f\t%4.0d\t%f\n",PSOR.w,PSOR.iter,PSOR.timeElapsed);
		
	// 	}
	// 	fclose(testPSOR);
	// }

	// printf("\n");	


	// {	int nw = 20;
	// 	double dw = 0.1;
	// 	double w0 = 0.1;
	// 	float wi;

	// 	FILE *testLSOR;
	// 	testLSOR = fopen("testLSOR.dat","w");

		
	// 	// fprintf(testLSOR,"Different values of W for Line SOR\n\n");
	// 	// fprintf(testLSOR,"  w     |  Iteration  | Time[sec] |\n");
	// 	// fprintf(testLSOR,"-----------------------------------\n");
	// 	for(int i=0; i<nw; i++) {
	// 		Elliptic LSOR;
	// 		wi = w0 + dw*( (float) i);
	// 		InitializePsi(&LSOR,wi);
	// 		LineSOR(&LSOR);
	// 		// fprintf(testLSOR, "%.3f\t|\t%4.0d\t  |\t%f  |\n",LSOR.w,LSOR.iter,LSOR.timeElapsed);
	// 		fprintf(testLSOR, "%.3f\t%4.0d\t%f\n",LSOR.w,LSOR.iter,LSOR.timeElapsed);
	// 	}
	// 	fclose(testLSOR);
	// }
	return 0;

}
