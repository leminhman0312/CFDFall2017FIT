#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
using namespace std;

struct heat {
    vector< vector<float> > u;
    vector< vector<float> > u_dummy_matrix;
	int imax;
	int jmax;
	double deltax;
	double deltay;
	double deltat;
};

void thomasTriDiagonalX(int j, double a[],double b[],double c[], double d[], struct heat *heat){
	int imax = heat->imax-1;

    double dprime[imax+1];
	double cprime[imax+1];
	dprime[1] = d[1];
	cprime[1] = c[1];

	//FORWARD LOOP
	for (int i = 2; i<=imax;i++){
		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));
	}

	heat->u_dummy_matrix[imax][j] = cprime[imax]/dprime[imax];
	
	//BACKWARD LOOP
	for (int i = imax-1;i>=2;i--){
		heat->u_dummy_matrix[i][j] = (cprime[i]-(a[i]*heat->u_dummy_matrix[i+1][j]))/(dprime[i]);
	}
}

void thomasTriDiagonalY(int i, double a[],double b[],double c[], double d[], struct heat *heat){
	int jmax = heat->jmax-1;

    double dprime[jmax+1];
	double cprime[jmax+1];
	dprime[1] = d[1];
	cprime[1] = c[1];

	//FORWARD LOOP
	for (int j = 2; j<=jmax;j++){

		dprime[j] = d[j] - ((b[j]*a[j-1])/(dprime[j-1]));
		cprime[j] = c[j] - ((cprime[j-1]*b[j])/(dprime[j-1]));
		
	}

	heat->u_dummy_matrix[i][jmax] = cprime[jmax]/dprime[jmax];
	
	//BACKWARD LOOP
	for (int j = jmax-1;j>=2;j--){
		heat->u_dummy_matrix[i][j] = (cprime[j]-(a[j]*heat->u_dummy_matrix[i][j+1]))/(dprime[j]);
	}
}







int main(){
	//Define basic parameters
	double dt = 0.01;
	double hours = 24; //hour
	double nmax_double = hours/dt;
	int nmax = int(nmax_double);
	double deltax = 0.1; 
	double deltay = 0.1;
	double alpha = 0.645; //thermal diffusivity
	double maxlength_X = 3.5; //3.5 by 3.5
	double minlength_X = 0.0; //assume start at 0
	double maxlength_Y = 3.5; //
	double minlength_Y = 0.0; //assume start at 0
	double numberofpointsX = ((maxlength_X-minlength_X)/(deltax))+1.;
	double numberofpointsY = ((maxlength_Y-minlength_Y)/(deltay))+1.;
	int numberofpointsX_int = ceil(numberofpointsX);
	int numberofpointsY_int = ceil(numberofpointsY);
	int imax = numberofpointsX_int;
	int jmax = numberofpointsY_int;
	double t0 = 0.0; //initial temperature
	double t1 = 200.00; //boundary temperature
	double t2 = 200.00; //boundary temperature
	double t3 = 0.0;//initial temperature
	double t4 = 0.0;//initial temperature

  heat heat;

  heat.imax = imax;
  heat.jmax = jmax;
  heat.deltax = deltax;
  heat.deltay = deltay;
  heat.deltat = dt;

	heat.u.resize(imax+1);
	heat.u_dummy_matrix.resize(imax+1);
	heat.u[0].resize(jmax+1);
	for(int i=1; i<=imax; i++) {
		heat.u[i].resize(jmax+1);
		heat.u_dummy_matrix[i].resize(jmax+1);
		heat.u[i][1] = t2;
		for(int j = 1; j<=jmax; j++){
			heat.u[1][j] = t1;
		}
	}

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
	for (int i=1; i<=imax;i++){
		ax[i] = -d1;
		bx[i] = -d1;
		diagonalX[i] = plus2d1;
	}
    ax[1] = 0.0;
    ax[imax] = 0.0;
    bx[1] = 0.0;
    bx[imax] = 0.0;
    diagonalX[1] = 1.;
    diagonalX[imax] = 1.;

	//Put values into Udummy = U
	for (int i = 1; i<=imax;i++){
		for (int j = 1; j<=jmax;j++){
			heat.u_dummy_matrix[i][j] = heat.u[i][j];
		}
	}

	//Define vectors for Thomas in Y

	double ay[imax+1]={0.0}; //above
	double by[imax+1]={0.0}; //below
	double cy[imax+1];//rhs
	double diagonalY[imax+1]={0.0};//diagonal

	//Fill out values
	//For ax,bx
	for (int i=1; i<=jmax;i++){
		ay[i] = -d2;
		by[i] = -d2;
		diagonalY[i] = plus2d2;
	}
    ay[1] = 0.0;
    ay[jmax] = 0.0;
    by[1] = 0.0;
    by[imax] = 0.0;
    diagonalY[1] = 1.;
    diagonalY[jmax] = 1.;

    cx[1] = t2;
    cx[imax] = t4;
    cy[1] = t1;
    cy[imax] = t3;
    

    FILE *outfile;
    outfile = fopen("01resultsTables.txt","w");
    FILE *outfile3;
    outfile3 = fopen("05resultsTables.txt","w");   

	for (int n = 1; n<=nmax;n++){
        // Updating cx, embedded BC, X sweep
        for (int j = 2;j<=jmax-1;j++){
            for (int i = 2;i<=imax-1;i++){
                if (i == 2){
                    cx[i] = d2*heat.u[i][j+1] + (minus2d2*heat.u[i][j])+(d2*heat.u[i][j-1])+(d1*t2); //BC
                } 
                else if (i==imax-1){
                    cx[i] = d2*heat.u[i][j+1] + (minus2d2*heat.u[i][j])+(d2*heat.u[i][j-1])+(d1*t4);
                }
                else {
                    cx[i] = d2*heat.u[i][j+1] + (minus2d2*heat.u[i][j])+(d2*heat.u[i][j-1]);
                }
                cx[i] = d2*heat.u[i][j+1] + (minus2d2*heat.u[i][j])+(d2*heat.u[i][j-1]);
                
            }
            thomasTriDiagonalX(j,ax,bx,cx,diagonalX,&heat);
	    }	   

        //update cy,embedded BC, Y sweep
        for (int i = 2;i<=imax-1;i++){
            for (int j = 2;j<=jmax-1;j++){
                if (j == 2){
                    cy[j] = d1*heat.u[i+1][j] + (minus2d1*heat.u[i][j])+(d1*heat.u[i-1][j])+(d2*t1); //BC
                } 
                else if (j == jmax-1){
                    cy[j] = d1*heat.u[i+1][j] + (minus2d1*heat.u[i][j])+(d1*heat.u[i-1][j])+(d2*t3); //BC				
                }
                else {
                    cy[j] = d1*heat.u[i+1][j] + (minus2d1*heat.u[i][j])+(d1*heat.u[i-1][j]);
                }
                cy[j] = d1*heat.u_dummy_matrix[i+1][j] + (minus2d1*heat.u_dummy_matrix[i][j])+(d1*heat.u_dummy_matrix[i-1][j]);
            }
            thomasTriDiagonalY(i,ay,by,cy,diagonalY,&heat);
        }
        for(int i=1; i<=imax; i++) {
					for(int j=1; j<=jmax; j++) {
		                heat.u[i][j] = heat.u_dummy_matrix[i][j];
					}
				}      

       //Tables to print, change "n" to desire time: 10= 0.1hr, 20 = 0.2hr,... etc

				if (n == 10){
			
	        fprintf(outfile,"Table show temperature distribution at t = %f hours\n",n*dt);
	        fprintf(outfile,"\n");
	        fprintf(outfile,"\t");
	        for(int j=1; j<=jmax+1; j++) {
	            if(((j-1)%5)==0) {
	                fprintf(outfile,"\t | %8.3f",(j-1)*deltay);
	            }
	        }
	        fprintf(outfile, "  |");
	        fprintf(outfile,"\n");
	        fprintf(outfile, "---------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n" );
	        fprintf(outfile, "         |           |           |           |           |           |           |           |           |");
	        fprintf(outfile,"\n");
	        for(int i=1; i<=imax; i++) {
	          fprintf(outfile,"%.3f    |\t",(((i)*deltax)-deltax));	              
            for(int j=1; j<=jmax+1; j++) {
                if(((j-1)%5)==0) {
                    fprintf(outfile,"%7.3f  |\t",heat.u[i][j]);
                }
            }
            fprintf(outfile,"\n");
        	}
        }	

        if (n == 40){
			
	        fprintf(outfile3,"Table show temperature distribution at t = %f hours\n",n*dt);
	        fprintf(outfile3,"\n");
	        fprintf(outfile3,"\t");
	        for(int j=1; j<=jmax+1; j++) {
	            if(((j-1)%5)==0) {
	                fprintf(outfile3,"\t | %8.3f",(j-1)*deltay);
	            }
	        }
	        fprintf(outfile3, "  |");
	        fprintf(outfile3,"\n");
	        fprintf(outfile3, "---------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|-----------|\n" );
	        fprintf(outfile3, "         |           |           |           |           |           |           |           |           |");
	        fprintf(outfile3,"\n");
	        for(int i=1; i<=imax; i++) {
	          fprintf(outfile3,"%.3f    |\t",(((i)*deltax)-deltax));	              
            for(int j=1; j<=jmax+1; j++) {
                if(((j-1)%5)==0) {
                    fprintf(outfile3,"%7.3f  |\t",heat.u[i][j]);
                }
            }
            fprintf(outfile3,"\n");
        	}
        }
     }

     		//This is to print entire tables for plotting.
     		//Change hours for time step
     		//For ease of import to Matlab, only temperature profiles are printed
     		//X and Y are generated in script
        FILE *outfile2;
        outfile2 = fopen("CFD_hw1q2_ResultsPlottingMatlab.txt","w"); 
				fprintf(outfile2,"\n");
				for(int i=1; i<=imax; i++) {	          
				for(int j=1; j<=jmax; j++) {
				fprintf(outfile2,"%f\t",heat.u[i][j]);
				}
				fprintf(outfile2,"\n");
				}

}






	
	







