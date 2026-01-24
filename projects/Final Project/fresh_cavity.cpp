#include <cmath> 
#include <stdio.h>

//Thomas Algorithm
void thomasTriDiagonal(int j, int imax, double a[],double b[],double c[], double d[], double u[][100]){

	double dprime[imax];
	double cprime[imax];
	dprime[1] = d[1];
	cprime[1] = c[1];

	//FORWARD LOOP
	for (int i = 2; i<=imax;i++){
		dprime[i] = d[i] - ((b[i]*a[i-1])/(dprime[i-1]));
		cprime[i] = c[i] - ((cprime[i-1]*b[i])/(dprime[i-1]));
	}

	u[j][imax] = cprime[imax]/dprime[imax];
	
	//BACKWARD LOOP
	for (int i = imax-1;i>=2;i--){		
		u[j][i] = (cprime[i]-(a[i]*u[j][i+1]))/(dprime[j]);
	}  
}



void Compute_Dilatation(int imax, int jmax, double dilation_N[][100], double u[][100],double v[][100]){
		

	double dx = dx;
	double dy = dy;

	// INTERIOR NODES 
	// ALL 2ND CENTRAL DIFF 

	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){
			dilation_N[j][i] = ((u[j][i+1]-u[j][i-1])/(2.*dx))+
			((v[j+1][i]-v[j-1][i])/(2.*dy));
			}
	}


	

	for (int j = 2; j<=jmax-1;j++){
		//LEFT WALL
		//FORWARD IN U, CENTRAL IN V 

		dilation_N[j][1] = ((-u[j][3]+4.*u[j][2]-3.*u[j][1])/(2.*dx))+((v[j+1][1]-v[j-1][1])/(2.*dy));

		//RIGHT WALL
		//BACKWARD IN U, CENTRAL IN V 
		dilation_N[j][imax] = ((u[j][imax-2]-4.*u[j][imax-1]+3.*u[j][imax])/(2.*dx))+
		((v[j+1][imax]-v[j-1][imax])/(2.*dy));
	}


	

	for (int i = 2; i<=imax-1;i++){
		//TOP WALL
		//CENTRAL IN U, FORWARD IN V

		dilation_N[1][i] = ((u[1][i+1]-u[1][i-1])/(2.*dx))+
		((-v[3][i]+4.*v[2][i]-3.*v[1][i])/(2.*dy));

		//BOTTOM WALL
		//CENTRAL IN U, BACKWARD IN V
		dilation_N[jmax][i] = (u[jmax][i+1]-u[jmax][i-1])/(2.*dx)+
		((v[jmax-2][i]-4.*v[jmax-1][i]+3.*v[jmax][i])/(2.*dy));
	}


	//CORNERS

	//TOP LEFT
	//FORWARD IN U, FORWARD IN V 
	dilation_N[1][1] = ((-u[1][3]+4.*u[1][2]-3.*u[1][1])/(2.*dx))+
	((-v[3][1]+4.*v[2][1]-3.*v[1][1])/(2.*dx));

	//TOP RIGHT
	//BACKWARD IN U, FORWARD IN V

	dilation_N[1][imax] = ((u[1][imax-2]-4.*u[1][imax-1]+3.*u[1][imax])/(2.*dx))+
	((-v[3][imax]+4.*v[2][imax]-3.*v[1][imax])/(2.*dx));

	//BOTTOM LEFT
	//FORWARD IN U, BACKWARD IN V
	dilation_N[jmax][1] = ((-u[jmax][3]+4.*u[jmax][2]-3.*u[jmax][1])/(2.*dx))+
	((v[jmax-2][1]-4.*v[jmax-1][1]+3.*v[jmax][1])/(2.*dy));

	//BOTTOM RIGHT
	//BACKWARD IN U, BACKWARD IN V
	dilation_N[jmax][imax] = ((u[jmax][imax-2]-4.*u[jmax][imax-1]+3.*u[jmax][imax])/(2.*dx))+((v[jmax-2][imax]-4.*v[jmax-1][imax]+3.*v[jmax][imax])/(2.*dy));

}




void LineSOR(double Re, int imax, int jmax, double dx, double dy, double dt, double w_optimal, double lsormax, double p_updated[][100],double dilation_N[][100],double u[][100],double v[][100], double p[][100]){
	double beta = dx/dy;
	double BetaSquare = pow(beta,2.);
	





	//THOMAS PARAMETERS
	double ax[imax+1]; //above
	double bx[imax+1]; //below
	double cx[imax+1];//rhs
	double diagonalX[imax+1];//diagonal
	

	
	

  //Special values rewrite 

	// // For A
  // ax[imax] = 1.;
	// ax[imax-1] = 0.0;
  // ax[1]= 0.;
  // ax[0] = 0.;


	// //For B
	// bx[0] = 0.;
  // bx[1] = -1.;
  // bx[imax] = 0.;
	// bx[2] = 0.0;
  
	// //For C 

	// cx[1] = 0.;
  
	// //For D
	// diagonalX[0] = 0.;
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
	

  


	for(int lsoriter=0; lsoriter <= lsormax; lsoriter++) {
		double error = 0.0;

		//loop through jTH
		for (int  j = 2; j<=jmax-1;j++){
			//loop through iTH
			for (int i = 2; i<=imax-1;i++){			


				//SET UP THOMAS VECTORS

					//A- thomas
					ax[i] = w_optimal;
					
					//B-thomas
					bx[i] = w_optimal;

					//D-thomas
					diagonalX[i] = -2.*(1.+pow(beta,2.));

					//C - VECTORS
					first[i] = -(1.-w_optimal)*(2.*(1+pow(beta,2.)))*p[j][i];

					second[i] = (w_optimal*BetaSquare)*(p[j+1][i]+p_updated[j-1][i]);

					third[i] = (dilation_N[j][i])/dt;

					forth[i] = (pow((u[j][i+1]),2.)-2.*pow((u[j][i]),2.)+pow((u[j][i-1]),2.))/(pow(dx,2.));

					fifth[i] = (pow((v[j+1][i]),2.)-2.*pow((v[j][i]),2.)+pow((v[j-1][i]),2.))/(pow(dy,2.));

					sixth[i] = (1./Re)*((dilation_N[j][i+1]-2.*dilation_N[j][i]+dilation_N[j][i-1])/(pow(dx,2.)))+((dilation_N[j+1][i]-2.*dilation_N[j][i]+dilation_N[j-1][i])/(pow(dy,2.)));

					seventh[i] = (((u[j+1][i+1]*v[j+1][i+1])+(u[j-1][i-1]*v[j-1][i-1])-(u[j+1][i-1]*v[j+1][i-1])-(u[j-1][i+1]*v[j-1][i+1])))/(2.*dx*dy);


					cx[i] = (first[i]-second[i]+third[i]-forth[i]-fifth[i]+sixth[i]-seventh[i]);		

										
			}
      
			//BOUNDARY CONDITIONS
      //dPsi/dX = 0 on left and right wall 
      diagonalX[1] = diagonalX[1]+bx[1];
			// bx[2] = 0.0;
      diagonalX[imax] = diagonalX[imax]+ax[imax];
			// ax[imax-1] = 0.0;
      
			//CALL THOMAS
			thomasTriDiagonal(j,imax, ax,bx,cx,diagonalX,p_updated);
			

		
			//CALCULATE ERROR
			for (int i = 1; i<=imax;i++){
				error += abs(p_updated[j][i]-p[j][i]);
			}

		}

		p_updated = p;
		
		//UPDATE ITERATION COUNTER

		lsoriter = lsoriter +1;
    // printf("%d\n",lsoriter);
    if(error <= errormax){
      
			break;
		}
    convergenceIter = lsoriter;
	}
}











int main(){

  int imax = 11;
  int jmax = 11; 
  int w = 1.315;
  double dt = 0.003;
  double dx = 0.00625;
  double dy = 0.00625;
  double u0 = 1.;
  double p0 = 3350;
  int tmax = 10;
  double errormax = 0.001;
  double lsormax = 1000;


  double u[jmax][imax];
  double u_updated[jmax][imax];
  double v[jmax][imax];
  double v_updated[jmax][imax];
  double p[jmax][imax];
  double p_updated[jmax][imax];
  double dilatation[jmax][imax];

  for (int i = 1; i<=imax;i++){
    u_updated[1][i] = u0;
    u[1][i] = u0; //top layer 

		p_updated[1][i] = p0;//top
		p_updated[jmax][i] = p0; //bottom
		p[1][i] = p0; //top 
		p[jmax][i] = p0; //bottom
  }

}
