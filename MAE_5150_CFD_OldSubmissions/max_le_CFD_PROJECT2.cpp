//MAE 5150: Coding Project 2
//Max Le
//MAIN PROGRAM
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "CFD_PROJECT2_HEADERS.h"
#include <ctime>
using namespace std;
float wi;//dummy variable to test W later on

int main(){

	//CALL POINT GAUSS SEIDEL
	// Elliptic PGS;
	// InitializePsi(&PGS,wi);
	// PointGaussSeidel(&PGS);

	//CALL LINE GAUSS SEIDEL
	// Elliptic LGS;
	// InitializePsi(&LGS,wi);
	// LineGaussSeidel(&LGS);

	//CALL POINT SOR AT OPTIMUM VALUE OF W 
	{
		float wi = 1.801;	
		Elliptic PSOR;
		InitializePsi(&PSOR,wi);
		PointSOR(&PSOR);
	}

	//CALL LINE SOR AT OPTIMUM VALUE OF W
	{
		float wi = 1.3;
		Elliptic LSOR;
		InitializePsi(&LSOR,wi);
		LineSOR(&LSOR);
	}

	//CALL THIS FUNCTION TO SHOW TABLES OF W VS. ITERATION
	//IF CALL THIS FUNCTION, THEN COMMENT OUT THE INDIVIDUAL PSOR/LSOR AT OPTIMUM W 
	// PrintTablesAndCompare();	
	return 0;
}
//END MAIN PROGRAM