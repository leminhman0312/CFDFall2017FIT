#include <iostream>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;


struct NavierStoke{
	double dx;
	double dy;
	double Re; //Reynolds number
	double p0; //constant pressure for top and bottom surface
	double u0; //top velocity of top surface
	int tmax; //max time
	int imax; 
	int jmax;
  int convergenceIter;
  double nmax; //max time step
	double dt; 
	double w_optimal; //omega for LineSOR
	double errormax;//error max for SOR
	double error;
	int lsormax;
	vector< vector<double> > u; //U velocity component 
	vector< vector<double> > v; //V velocity component
	vector< vector<double> > u_updated; //Holder for U velocity component 
	vector< vector<double> > v_updated; //Holder for V velocity component
	vector< vector<double> > p; //Pressure field
	vector< vector<double> > p_updated; //Pressure field holder
	vector< vector<double> > dilation_N; //ensure continuity at N

};


void Initialize(struct NavierStoke *cavity){
	
	//DECLARE BASIC PARAMETERS
	cavity->dx = 0.00625;
	cavity->dy = 0.00625;
	cavity->Re = 1000; 
  cavity->nmax = 15000;
	cavity->p0 = 3350; 
	cavity->u0 = 1.; 	
	cavity->imax = 11; 
	cavity->jmax = 11;
	cavity->dt = 0.003; 
  cavity->tmax = 10;
	// cavity->tmax = cavity->nmax*cavity->dt;
	cavity->w_optimal = 1.; 
	cavity->errormax = 0.001;
	cavity->lsormax = 1000;

	//DECLARE FIELDS

	cavity->u.resize(cavity->jmax+1);
	cavity->u_updated.resize(cavity->jmax+1);
	cavity->v.resize(cavity->jmax+1);
	cavity->v_updated.resize(cavity->jmax+1);
	cavity->p.resize(cavity->jmax+1);
	cavity->p_updated.resize(cavity->jmax+1);
	cavity->dilation_N.resize(cavity->jmax+1);
	for (int j = 1; j<=cavity->jmax;j++){
		cavity->u[j].resize(cavity->imax+1);
		cavity->u_updated[j].resize(cavity->imax+1);
		cavity->v[j].resize(cavity->imax+1);
		cavity->v_updated[j].resize(cavity->imax+1);
		cavity->p[j].resize(cavity->imax+1);
		cavity->p_updated[j].resize(cavity->imax+1);
		cavity->dilation_N[j].resize(cavity->imax+1);
	}


	//zeros out all fields
	for (int j = 1; j<=cavity->jmax;j++){
		for (int i = 1; i<=cavity->imax;i++){
			cavity->u[j][i] = 0.0;
			cavity->u_updated[j][i] = 0.0;
			cavity->v[j][i] = 0.0;
			cavity->v_updated[j][i] = 0.0;
			cavity->p[j][i] = 0.0;
			cavity->p_updated[j][i] = 0.0;
		}
	}

	//Top lid moving at u = u0;
	//Bottom and Top pressures = p0
	
	for (int i = 1; i<=cavity->imax;i++){
		cavity->u_updated[1][i] = cavity->u0; //top layer
		cavity->u[1][i] = cavity->u0; //top layer 

		cavity->p_updated[1][i] = cavity->p0;//top
		cavity->p_updated[cavity->jmax][i] = cavity->p0; //bottom
		cavity->p[1][i] = cavity->p0; //top 
		cavity->p[cavity->jmax][i] = cavity->p0; //bottom

	}
}

//Thomas Algorithm
void thomasTriDiagonal(int j, double a[],double b[],double c[], double d[], struct NavierStoke *cavity){

	int imax = cavity->imax;
	double dprime[imax+1];
	double cprime[imax+1];
	dprime[2] = d[2];
	cprime[2] = c[2];

	//FORWARD LOOP
	for (int i = 3; i<=imax-1;i++){
		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));
	}

	cavity->p_updated[j][imax-1] = cprime[imax-1]/dprime[imax-1];
	
	//BACKWARD LOOP
	for (int i = imax-2;i>=3;i--){
		cavity->p_updated[j][i] = (cprime[j]-(a[j]*cavity->p_updated[j][i+1]))/(dprime[j]);
	}  
}




double forwardCentralSecond(double plus2, double plus1, double same, double delta){
	return (-plus2+4.*plus1-3.*same)/(2.*delta);
}

double backwardCentralSecond(double minus2, double minus1, double same, double delta){
	return (minus2-4.*minus1+3.*same)/(2.*delta);
}

double square_func(double a){
	return pow(a,2.);
}

double SecondOrderCentral(double plus, double minus, double delta){
	return (plus-minus)/(2.*delta);
}

double SecondDerivativeOrderCentral(double plus, double minus, double same, double delta){
	return (plus-2.*same+minus)/(pow(delta,2.));
}


double FowardDiff(double plus, double same, double delta){
	return (plus-same)/(delta);
}

double BackwardDiff(double same, double minus, double delta){
	return (same-minus)/delta;
}

//CALCULATE DILATATION

void Compute_Dilatation(struct NavierStoke *cavity){
	int imax = 	cavity->imax;
	int jmax = cavity->jmax;

	

	double dx = cavity->dx;
	double dy = cavity->dy;
	//Interior nodes
	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){
			cavity->dilation_N[j][i] = SecondOrderCentral(cavity->u[j][i+1],cavity->u[j][i-1],dx)+SecondOrderCentral(cavity->v[j+1][i],cavity->v[j-1][i],dy);
			}
	}


	//Left wall and Right wall

	for (int j = 2; j<=jmax-1;j++){
		//LEFT WALL

		cavity->dilation_N[j][1] = forwardCentralSecond(cavity->u[j][3],cavity->u[j][2],cavity->u[j][1],dx)+SecondOrderCentral(cavity->v[j+1][1],cavity->v[j-1][1],dy);

		//RIGHT WALL

		cavity->dilation_N[j][imax] = backwardCentralSecond(cavity->u[j][imax-2],cavity->u[j][imax-1],cavity->u[j][imax],dx)+SecondOrderCentral(cavity->v[j+1][imax],cavity->v[j-1][imax],dy);
	}


	// Top and Bottom wall

	for (int i = 2; i<=imax-1;i++){
		//TOP WALL

		cavity->dilation_N[1][i] = SecondOrderCentral(cavity->u[1][i+1],cavity->u[1][i-1],dx)+forwardCentralSecond(cavity->v[3][i],cavity->v[2][i],cavity->v[1][i],dy);

		//BOTTOM WALL

		cavity->dilation_N[jmax][i] = SecondOrderCentral(cavity->u[jmax][i+1],cavity->u[jmax][i-1],dx)+backwardCentralSecond(cavity->v[jmax-2][i],cavity->v[jmax-1][i],cavity->v[jmax][i],dy);
	}


	//CORNERS

	//TOP LEFT

	cavity->dilation_N[1][1] = forwardCentralSecond(cavity->u[1][3],cavity->u[1][2],cavity->u[1][1],dx)+forwardCentralSecond(cavity->v[3][1],cavity->v[2][1],cavity->v[1][1],dx);

	//TOP RIGHT

	cavity->dilation_N[1][imax] = backwardCentralSecond(cavity->u[1][imax-2],cavity->u[1][imax-1],cavity->u[1][imax],dx)+forwardCentralSecond(cavity->v[3][imax],cavity->v[2][imax],cavity->v[1][imax],dx);

	//BOTTOM LEFT

	cavity->dilation_N[jmax][1] = forwardCentralSecond(cavity->u[jmax][3],cavity->u[jmax][2],cavity->u[jmax][1],dx)+backwardCentralSecond(cavity->v[jmax-2][1],cavity->v[jmax-1][1],cavity->v[jmax][1],dy);

	//BOTTOM RIGHT

	cavity->dilation_N[jmax][imax] = backwardCentralSecond(cavity->u[jmax][imax-2],cavity->u[jmax][imax-1],cavity->u[jmax][imax],dx)+backwardCentralSecond(cavity->v[jmax-2][imax],cavity->v[jmax-1][imax],cavity->v[jmax][imax],dy);

}




//LineSOR- SOLVING THE PRESSURE POISSON EQUATION
void LineSOR(struct NavierStoke *cavity){
	double beta = cavity->dx/cavity->dy;
	double BetaSquare = pow(beta,2.);
	int imax = cavity->imax;
	int jmax = cavity->jmax;





	//THOMAS PARAMETERS
	double ax[imax+1]; //above
	double bx[imax+1]; //below
	double cx[imax+1];//rhs
	double diagonalX[imax+1];//diagonal
	double dx = cavity->dx;
	double dy = cavity->dy;

	
	//Fill out values
	for (int i=1; i<=imax;i++){
		ax[i] = cavity->w_optimal;
		bx[i] = cavity->w_optimal;
		diagonalX[i] = -2.*(1.+pow(beta,2.));

	}

  //Special values rewrite 

	// For A
  ax[imax] = 1.;
	ax[imax-1] = 0.0;
  ax[1]= 0.;
  ax[0] = 0.;


	//For B
	bx[0] = 0.;
  bx[1] = -1.;
  bx[imax] = 0.;
	bx[2] = 0.0;
  
	//For C 

	cx[1] = 0.;
  
	//For D
	diagonalX[0] = 0.;
  diagonalX[1] = 1.;
	diagonalX[imax] = 1.;




	//Terms for the RHS 

	double first[imax+1];
	double second[imax+1];
	double third[imax+1];
	double forth[imax+1];
	double fifth[imax+1];
	double sixth[imax+1];
	double seventh[imax+1];
	

  


	for(int lsoriter=0; lsoriter <= cavity->lsormax; lsoriter++) {
		cavity->error = 0.0;

		//loop through jTH
		for (int  j = 2; j<=jmax-1;j++){
			//loop through iTH
			for (int i = 2; i<=imax-1;i++){			


				//EMBEDDED BOUNDARY CONDITIONS
					first[i] = -(1.-cavity->w_optimal)*(2.*(1+pow(beta,2.)))*cavity->p[j][i];

					second[i] = (cavity->w_optimal*BetaSquare)*(cavity->p[j+1][i]+cavity->p_updated[j-1][i]);

					third[i] = (cavity->dilation_N[j][i])/cavity->dt;

					forth[i] = SecondDerivativeOrderCentral(square_func(cavity->u[j][i+1]),square_func(cavity->u[j][i-1]),square_func(cavity->u[j][i]),dx);

					fifth[i] = SecondDerivativeOrderCentral(square_func(cavity->v[j+1][i]),square_func(cavity->v[j-1][i]),square_func(cavity->v[j][i]),dy);

					sixth[i] = ((1./cavity->Re)*(SecondDerivativeOrderCentral(cavity->dilation_N[j][i+1],cavity->dilation_N[j][i-1],cavity->dilation_N[j][i],dx)+ 	SecondDerivativeOrderCentral(cavity->dilation_N[j+1][i],cavity->dilation_N[j-1][i],cavity->dilation_N[j][i],dy)));

					seventh[i] = (((cavity->u[j+1][i+1]*cavity->v[j+1][i+1])+(cavity->u[j-1][i-1]*cavity->v[j-1][i-1])-(cavity->u[j+1][i-1]*cavity->v[j+1][i-1])-(cavity->u[j-1][i+1]*cavity->v[j-1][i+1])))/(2.*dx*dy);


					cx[i] = (first[i]-second[i]+third[i]-forth[i]-fifth[i]+sixth[i]-seventh[i]);								
			}
      
      //dPsi/dX = 0 on left and right wall 
      diagonalX[2] = diagonalX[2]+bx[2];
      diagonalX[imax-1] = diagonalX[imax-1]+ax[imax-1];
      
			//CALL THOMAS
			thomasTriDiagonal(j,ax,bx,cx,diagonalX,cavity);
			//CALCULATE ERROR

			cavity->p_updated[j][1] = cavity->p_updated[j][2];
			cavity->p_updated[j][imax] = cavity->p_updated[j][imax-1];

			for (int i = 1; i<=imax;i++){
				cavity->error += abs(cavity->p_updated[j][i]-cavity->p[j][i]);
			}

		}

		// for(int j=1; j<jmax; j++) {
		// 	cavity->p_updated[1][j] = cavity->p_updated[2][j];
		// 	cavity->p_updated[imax][j] = cavity->p_updated[imax-1][j];
		// }

	
		cavity->p = cavity->p_updated;
		
		//UPDATE ITERATION COUNTER

		lsoriter = lsoriter +1;
    // printf("%d\n",lsoriter);
    if(cavity->error < cavity->errormax){
      
			break;
		}
    cavity->convergenceIter = lsoriter;
	}
}



//MAIN MOMENTUM SOLVER-> GET U AND V 
void MomentumSolver(struct NavierStoke *cavity){
	int imax = cavity->imax;
	int jmax = cavity->jmax;
	double dt = cavity->dt;
	double dx = cavity->dx;
	double dy = cavity->dy;
	double reynolds = cavity->Re;
	double epsilonX; //ex for U equation
	double epsilonY; //ey for V equation

//CALCULATE U
	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){

			if (cavity->u[j][i]>=0){
				epsilonX == 1.;
			}
			else if (cavity->u[j][i]<0){
				epsilonX == -1.;
			}

			else if (cavity->p_updated[j][i]>=0){
				epsilonX == 1.;
			}
			else if (cavity->p_updated[j][i]<0){
				epsilonX == -1.;
			}			
			
			cavity->u_updated[j][i] = cavity->u[j][i]+dt*			
			
			
			(
				((1./reynolds)*(
					SecondDerivativeOrderCentral(cavity->u[j][i+1],cavity->u[j][i-1],cavity->u[j][i],dx)+
					SecondDerivativeOrderCentral(cavity->u[j+1][i],cavity->u[j-1][i],cavity->u[j][i],dy)))-
				
        (0.5*(1.+epsilonX)*BackwardDiff(pow(cavity->u[j][i],2.),pow(cavity->u[j][i-1],2.),dx))-
				
        (0.5*(1.-epsilonX)*FowardDiff(pow(cavity->u[j][i+1],2.),pow(cavity->u[j][i],2.),dx))-
				
        (0.5*(1.+epsilonX)*BackwardDiff(cavity->p_updated[j][i],cavity->p_updated[j][i-1],dx))-
				
        (0.5*(1.-epsilonX)*FowardDiff(cavity->p_updated[j][i+1],cavity->p_updated[j][i],dx))-
				
        (SecondOrderCentral(cavity->u[j+1][i]*cavity->v[j+1][i],cavity->u[j-1][i]*cavity->v[j-1][i],dy))
        
			);
		}
	}

	//CALCULATE V
	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){		

			
			if (cavity->v[j][i]>0){
				epsilonY == 1.;
			}

			else if (cavity->v[j][i]<0){
				epsilonY == -1.;
			}

			else if (cavity->p_updated[j][i]>=0){
				epsilonY == 1.;
			}
			else if (cavity->p_updated[j][i]<0){
				epsilonY == -1.;
			}			


			cavity->v_updated[j][i] = cavity->v[j][i] + dt*			
			
			(
				((1./reynolds)*(SecondDerivativeOrderCentral(cavity->v[j][i+1],cavity->v[j][i-1],cavity->v[j][i],dx)+SecondDerivativeOrderCentral(cavity->v[j+1][i],cavity->v[j-1][i],cavity->v[j][i],dy)))

      -(SecondOrderCentral(cavity->u[j][i+1]*cavity->v[j][i+1],cavity->u[j][i-1]*cavity->v[j][i-1],dx))      
      
      -(0.5*(1+epsilonY)*BackwardDiff(pow(cavity->v[j][i],2.),pow(cavity->v[j-1][i],2.),dy))
      
      -(0.5*(1-epsilonY)*FowardDiff(pow(cavity->v[j+1][i],2.),pow(cavity->v[j][i],2.),dy))
      
      -(0.5*(1+epsilonY)*BackwardDiff(cavity->p_updated[j][i],cavity->p_updated[j-1][i],dy))
      
      -(0.5*(1-epsilonY)*FowardDiff(cavity->p_updated[j+1][i],cavity->p_updated[j][i],dy))
			
			);

			

		}
	}
	
	// Update U and V
	cavity->u = cavity->u_updated;
	cavity->v = cavity->v_updated;
}





//Main code
int main(){
	NavierStoke cavity;
	Initialize(&cavity);


  for (int time = 0; time<=0;time++){
    printf("t = %d LSOR iteration = %d\n",time,cavity.convergenceIter);
		Compute_Dilatation(&cavity); //get D
		LineSOR(&cavity);// get P
		MomentumSolver(&cavity); //get U and V
	}


  FILE *outfile1;
  outfile1 = fopen("Yew_results.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile1,"%.2f\t",cavity.u_updated[j][i]);
		}
		fprintf(outfile1,"\n");
	}

  FILE *outfile2;
  outfile2 = fopen("Vee_results.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile2,"%.2f\t",cavity.v_updated[j][i]);
		}
		fprintf(outfile2,"\n");
	}

  FILE *outfile3;
  outfile3 = fopen("Pressure_results.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile3,"%.2f\t",cavity.p_updated[j][i]);
		}
		fprintf(outfile3,"\n");
	}






	return 0;
}