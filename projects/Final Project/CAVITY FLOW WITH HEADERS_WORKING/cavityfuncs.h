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
	cavity->p0 = 3350.; 
	cavity->u0 = 1.; 	
	cavity->imax = 81; 
	cavity->jmax = 81;
	cavity->dt = 0.003; 
  cavity->tmax = 10;
	// cavity->tmax = cavity->nmax*cavity->dt;
	cavity->w_optimal = 1.315; 
	cavity->errormax = 0.0001;
	cavity->lsormax = 100000000;

	//DECLARE FIELDS

	cavity->u.resize(cavity->jmax+1);
	cavity->u_updated.resize(cavity->jmax+1);
	cavity->v.resize(cavity->jmax+1);
	cavity->v_updated.resize(cavity->jmax+1);
	cavity->p.resize(cavity->jmax+1);
	cavity->p_updated.resize(cavity->jmax+1);
	cavity->dilation_N.resize(cavity->jmax+1);

	//INITIALIZE VECTORS
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
			cavity->dilation_N[j][i] = 0.0;
		}

	}

	int jmax = cavity->jmax;
	//INITIAL CONDITIONS 
	for (int i = 1; i<=cavity->imax;i++){
		cavity->u_updated[jmax][i] = cavity->u0; //top layer
		cavity->u[jmax][i] = cavity->u0; //top layer 

		cavity->p_updated[1][i] = cavity->p0;//top
		cavity->p_updated[cavity->jmax][i] = cavity->p0; //bottom
		cavity->p[1][i] = cavity->p0; //top 
		cavity->p[cavity->jmax][i] = cavity->p0; //bottom

	}

	cavity->p_updated = cavity->p;
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
	for (int i = imax-2;i>=2;i--){
		cavity->p_updated[j][i] = (cprime[i]-(a[i]*cavity->p_updated[j][i+1]))/(dprime[i]);
	}  
}

//CALCULATE DILATATION

void Compute_Dilatation(struct NavierStoke *cavity){
	int imax = 	cavity->imax;
	int jmax = cavity->jmax;	

	double dx = cavity->dx;
	double dy = cavity->dy;

	// INTERIOR NODES 
	// ALL 2ND CENTRAL DIFF 

	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){
			cavity->dilation_N[j][i] = (
                    (cavity->u[j][i+1]-
                    cavity->u[j][i-1])/(2.*dx)
                )
                +
			    (
                    (cavity->v[j+1][i]-
                    cavity->v[j-1][i])/(2.*dy)
                );
			}
	}	

	for (int j = 2; j<=jmax-1;j++){
		//LEFT WALL
		//FORWARD IN U, CENTRAL IN V 

		cavity->dilation_N[j][1] = (
                (
                    -cavity->u[j][3]+
                    4.*cavity->u[j][2]-
                    3.*cavity->u[j][1]
                )/(2.*dx)
            )+
            (
                (
                    cavity->v[j+1][1]-
                    cavity->v[j-1][1]
                )/(2.*dy)
            );

		//RIGHT WALL
		//BACKWARD IN U, CENTRAL IN V 
		cavity->dilation_N[j][imax] = (
                (
                    cavity->u[j][imax-2]-
                    4.*cavity->u[j][imax-1]+
                    3.*cavity->u[j][imax]
                )/(2.*dx)
            )+
		    (
                (
                    cavity->v[j+1][imax]-
                    cavity->v[j-1][imax]
                )/(2.*dy)
            );
	}	

	for (int i = 2; i<=imax-1;i++){
		//TOP WALL
		//CENTRAL IN U, FORWARD IN V

		cavity->dilation_N[1][i] = (
            (
                cavity->u[1][i+1]-
                cavity->u[1][i-1]
            )/(2.*dx)
        )+
		(
            (
                -cavity->v[3][i]+
                4.*cavity->v[2][i]-
                3.*cavity->v[1][i]
            )/(2.*dy)
        );

		//BOTTOM WALL
		//CENTRAL IN U, BACKWARD IN V
		cavity->dilation_N[jmax][i] = (
            (
                cavity->u[jmax][i+1]-
                cavity->u[jmax][i-1]
            )/(2.*dx)
        )+
		(
            (
                cavity->v[jmax-2][i]-
                4.*cavity->v[jmax-1][i]+
                3.*cavity->v[jmax][i]
            )/(2.*dy)
        );
	}
	//CORNERS

	//TOP LEFT
	//FORWARD IN U, FORWARD IN V 
	cavity->dilation_N[1][1] = ((-cavity->u[1][3]+4.*cavity->u[1][2]-3.*cavity->u[1][1])/(2.*dx))+
	((-cavity->v[3][1]+4.*cavity->v[2][1]-3.*cavity->v[1][1])/(2.*dx));

	//TOP RIGHT
	//BACKWARD IN U, FORWARD IN V

	cavity->dilation_N[1][imax] = ((cavity->u[1][imax-2]-4.*cavity->u[1][imax-1]+3.*cavity->u[1][imax])/(2.*dx))+
	((-cavity->v[3][imax]+4.*cavity->v[2][imax]-3.*cavity->v[1][imax])/(2.*dx));

	//BOTTOM LEFT
	//FORWARD IN U, BACKWARD IN V
	cavity->dilation_N[jmax][1] = ((-cavity->u[jmax][3]+4.*cavity->u[jmax][2]-3.*cavity->u[jmax][1])/(2.*dx))+
	((cavity->v[jmax-2][1]-4.*cavity->v[jmax-1][1]+3.*cavity->v[jmax][1])/(2.*dy));

	//BOTTOM RIGHT
	//BACKWARD IN U, BACKWARD IN V
	cavity->dilation_N[jmax][imax] = ((cavity->u[jmax][imax-2]-4.*cavity->u[jmax][imax-1]+3.*cavity->u[jmax][imax])/(2.*dx))+((cavity->v[jmax-2][imax]-4.*cavity->v[jmax-1][imax]+3.*cavity->v[jmax][imax])/(2.*dy));

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
	double dt = cavity->dt;

	//Terms for the RHS 

	double first;
	double second;
	double third;
	double fourth;
	double fifth;
	double sixth;
	double seventh; 


	for(int lsoriter=0; lsoriter <= cavity->lsormax; lsoriter++) {
		cavity->error = 0.0;

		//loop through jTH
		for (int  j = 2; j<=jmax-1;j++){
			//loop through iTH
			for (int i = 2; i<=imax-1;i++){	

				//SET UP THOMAS VECTORS

					//A- thomas
					ax[i] = cavity->w_optimal;
					
					//B-thomas
					bx[i] = cavity->w_optimal;

					//D-thomas
					diagonalX[i] = -2.*(1.+pow(beta,2.));

					//C - VECTORS
					first = -(1.-cavity->w_optimal)*(
                            2.*(1+pow(beta,2.))
                        )*cavity->p[j][i];

					second = (cavity->w_optimal*BetaSquare)*(cavity->p[j+1][i]+cavity->p_updated[j-1][i]);

					third = (cavity->dilation_N[j][i])/dt;

					fourth = (
                        pow((cavity->u[j][i+1]),2.)-
                        2.*pow((cavity->u[j][i]),2.)+
                        pow((cavity->u[j][i-1]),2.)
                    )/(pow(dx,2.));

					fifth =  (
                        pow((cavity->v[j+1][i]),2.)-
                        2.*pow((cavity->v[j][i]),2.)+
                        pow((cavity->v[j-1][i]),2.)
                    )/(pow(dy,2.));

					sixth = (1./cavity->Re)*(
                        (
                            (cavity->dilation_N[j][i+1]-
                            2.*cavity->dilation_N[j][i]+
                            cavity->dilation_N[j][i-1]
                            )/(pow(dx,2.))
                        )+(
                            (cavity->dilation_N[j+1][i]-
                            2.*cavity->dilation_N[j][i]+
                            cavity->dilation_N[j-1][i]
                            )/(pow(dy,2.))
                        )
                    );

					seventh = (
                        (
                            (cavity->u[j+1][i+1]*cavity->v[j+1][i+1])+
                            (cavity->u[j-1][i-1]*cavity->v[j-1][i-1])-
                            (cavity->u[j+1][i-1]*cavity->v[j+1][i-1])-
                            (cavity->u[j-1][i+1]*cavity->v[j-1][i+1])
                        )
                    )/(2.*dx*dy);

                    // printf("%.3f %.3f %.3f %.3f %.3f %.3f %.3f \n", first, second, third, fourth, fifth, sixth, seventh);

					cx[i] = (first-second)+(cavity->w_optimal*pow(dx,2.))*(third-fourth-fifth+sixth-seventh);						
			}
			ax[1] = 0.;
			bx[1] = 0.;
			ax[imax] = 0.;
			bx[imax] = 0.;
			diagonalX[1] = 0.;
			diagonalX[imax] = 0.;
			cx[1] = 0.;
			cx[imax] = 0.;

			//boundary conditions
			ax[2] = cavity->w_optimal;
			diagonalX[2] = -2.*(1.+pow(beta,2.)) + cavity->w_optimal;
			bx[2] = cavity->w_optimal;
			bx[imax-1] = cavity->w_optimal;
			diagonalX[imax-1] = -2.*(1.+pow(beta,2.)) + cavity->w_optimal;
			ax[imax-1] = cavity->w_optimal;

      
			
      
			//CALL THOMAS
			thomasTriDiagonal(j,ax,bx,cx,diagonalX,cavity);

			cavity->p_updated[j][1] = cavity->p_updated[j][2];
			cavity->p_updated[j][imax] = cavity->p_updated[j][imax-1];  

			// CALCULATE ERROR
			for (int i = 1; i<=imax;i++){
				cavity->error += pow(cavity->p_updated[j][i]-cavity->p[j][i],2.);
              
			}

		}
		cavity->error = sqrt(cavity->error);
		cavity->p = cavity->p_updated;
		
		//UPDATE ITERATION COUNTER

		lsoriter = lsoriter +1;
    // printf("%d\n",lsoriter);
    if(cavity->error <= cavity->errormax){
      
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

	double dsqudx, dsqudy, dusqdxp, dusqdxm, dpdxp, dpdxm;
	double duvdy, delsqu, dudt;

	double dsqvdx, dsqvdy, dvsqdyp, dvsqdym, dpdyp, dpdym;
	double duvdx, delsqv, dvdt;

// CALCULATE U
	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){

			if (cavity->u[j][i]>=0){
				epsilonX = 1.;					
			}else {
				epsilonX = -1.;				
			}

			dsqudx = (cavity->u[j][i+1]-
							2.*cavity->u[j][i]
							+cavity->u[j][i-1]
						)/(pow(dx,2.));
			dsqudy = (cavity->u[j+1][i]-
							2.*cavity->u[j][i]
							+cavity->u[j-1][i]
						)/(pow(dy,2.));
			delsqu = dsqudx + dsqudy;

			dusqdxm = (
				pow(cavity->u[j][i],2.)-
					pow(cavity->u[j][i-1],2.)
				)/dx;
			dusqdxp = (
				pow(cavity->u[j][i+1],2.)-
					pow(cavity->u[j][i],2.)
				)/dx;

			dpdxm = (
				cavity->p[j][i]-
				cavity->p[j][i-1]
			)/dx;
			dpdxp = (
				cavity->p[j][i+1]-
				cavity->p[j][i]
			)/dx;
			duvdy = (
				(cavity->u[j+1][i]*cavity->v[j+1][i])-
				(cavity->u[j-1][i]*cavity->v[j-1][i])
			)/(2.*dy);

			dudt = (1./reynolds)*delsqu-
				0.5*(1.+epsilonX)*dusqdxm-
				0.5*(1.-epsilonX)*dusqdxp-
				0.5*(1.+epsilonX)*dpdxm-
				0.5*(1.-epsilonX)*dpdxp-
				duvdy;
			cavity->u_updated[j][i] = cavity->u[j][i]+dt*dudt;
		}
	}



	//CALCULATE V 

	for (int j = 2; j<=jmax-1;j++){
		for (int i = 2; i<=imax-1;i++){
			if (cavity->v[j][i]>=0){
				epsilonY = 1.;
			}

			else{
			
				epsilonY = -1.;
			}

			dsqvdx = (cavity->v[j][i+1]-
							2.*cavity->v[j][i]
							+cavity->v[j][i-1]
						)/(pow(dx,2.));
			dsqvdy = (cavity->v[j+1][i]-
							2.*cavity->v[j][i]
							+cavity->v[j-1][i]
						)/(pow(dy,2.));
			delsqv = dsqvdx + dsqvdy;

			dvsqdym = (
				pow(cavity->v[j][i],2.)-
					pow(cavity->v[j-1][i],2.)
				)/dx;
			dvsqdyp = (
				pow(cavity->v[j+1][i],2.)-
					pow(cavity->v[j][i],2.)
				)/dy;

			dpdym = (
				cavity->p[j][i]-
				cavity->p[j-1][i]
			)/dy;
			dpdyp = (
				cavity->p[j+1][i]-
				cavity->p[j][i]
			)/dy;
			duvdx = (
				(cavity->u[j][i+1]*cavity->v[j][i+1])-
				(cavity->u[j][i-1]*cavity->v[j][i-1])
			)/(2.*dx);

            dvdt = (1./reynolds)*delsqv
				-0.5*(1.+epsilonY)*dvsqdym
				-0.5*(1.-epsilonY)*dvsqdyp
				-0.5*(1.+epsilonY)*dpdym
				-0.5*(1.-epsilonY)*dpdyp
				-duvdx;

			cavity->v_updated[j][i] = cavity->v[j][i]+dt*dvdt;
		}
	}

	// FORCING
	for (int j = 1; j<=jmax;j++){
		cavity->u_updated[j][1] = 0.0; //left wall
		cavity->u_updated[j][imax] = 0.0; //right wall
		cavity->v_updated[j][1] = 0.0; //left wall
		cavity->v_updated[j][imax] = 0.0; //right wall

	}

	for (int i = 1; i<=imax;i++){
		cavity->u_updated[jmax][i] = cavity->u0; //bottom wall
		cavity->v_updated[jmax][i] = 0.0; //bottom wall

		cavity->u_updated[1][i] = 0.0;
		cavity->v_updated[1][i] = 0.0;

	}

	// Update U and V
	cavity->u = cavity->u_updated;
	cavity->v = cavity->v_updated;
}