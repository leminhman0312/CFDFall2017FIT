#include "cavityfuncs.h"
using namespace std;

NavierStoke cavity;	
//Main code
int main(){
	
	Initialize(&cavity);
	int tmax=20000;

  for (int time = 0; time<=tmax ;time++){
		Compute_Dilatation(&cavity); //get D
		LineSOR(&cavity);// get P
		MomentumSolver(&cavity); //get U and V
		printf("t = %d LSOR %.6f \t iteration = %d\n",time,cavity.error, cavity.convergenceIter);
	}


  FILE *outfile1;
  outfile1 = fopen("Yew_results_15000.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile1,"%.2f\t",cavity.u_updated[j][i]);
		}
		fprintf(outfile1,"\n");
	}

  FILE *outfile2;
  outfile2 = fopen("Vee_results_15000.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile2,"%.2f\t",cavity.v_updated[j][i]);
		}
		fprintf(outfile2,"\n");
	}

  FILE *outfile3;
  outfile3 = fopen("Pressure_results_15000.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile3,"%.2f\t",cavity.p_updated[j][i]);
		}
		fprintf(outfile3,"\n");
	}

  FILE *outfile4;
  outfile4 = fopen("Dilation.dat","w");
	for (int j = 1;j<=cavity.jmax;j++){
		for (int i = 1; i<=cavity.imax;i++){
			fprintf(outfile4,"%.2f\t",cavity.dilation_N[j][i]);
		}
		fprintf(outfile4,"\n");
	}

	fcloseall;
	return 0;
}