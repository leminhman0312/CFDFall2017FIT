#include <stdio.h>
#include <math.h>


double myfunction(double x){

	
	double a  = (M_PI*x)/(2.0);
	double f = sin(a);
	return f;

}

//EXACT DERiVATIVES

double exactprime(double x){
	double fprime = (M_PI/2.0)*cos(M_PI*x/2.0);
}
double exactdoubleprime(double x){
	double fdoubleprime = (-1.0)*(pow(M_PI,2.0)/4.0)*sin(M_PI*x/2.0);
}

double exacttripleprime(double x){
	double ftripleprime = (-1.0)*(pow(M_PI,3.0)/8.0)*cos(M_PI*x/2.0);
}
double exactforthprime(double x){
	double fforthprime = (pow(M_PI,4.0)/16.0)*sin(M_PI*x/2.0);
}





//FINITE DIFFERENCES
double prime(double dx, double x){
	double dfdx = (myfunction(x+dx) - myfunction(x-dx))/(2.0*dx);
	
}

double doubleprime(double dx, double x){
	double d2fdx2 = (myfunction(x+dx)-2.0*myfunction(x)+myfunction(x-dx))/(pow(dx,2.0));
}

double tripleprime(double dx, double x){
	double d3fdx3 = (myfunction(x+(2.0*dx))-(2.0*myfunction(x+dx))+(2.0*myfunction(x-dx))-myfunction(x-(2.0*dx)))/(2.0*(pow(dx,3.0)));
}

double quadprime(double dx, double x){
	double d4fdx4 = (myfunction(x+(2.0*dx))-(4.0*myfunction(x+dx))+(6.0*myfunction(x))-(4.0*myfunction(x-dx))+myfunction(x-(2.0*dx)))/(pow(dx,4.0));
}

int main(){

	FILE *outfile;

	outfile = fopen("max_le_cdf_hwd1_question7.txt", "w");

	double x0 = 1.5;
	double analyticalprime = exactprime(1.5);
	double analyticaldoubleprime = exactdoubleprime(1.5);
	double analyticaltripleprime = exacttripleprime(1.5);
	double analyticalforthprime = exactforthprime(1.5);

	fprintf(outfile, "TABLE OF EXACT DERIVATES AT X = 1.5\n\n" );
	fprintf(outfile, "ORDER\t| ANALYTICAL|\n");
	fprintf(outfile, "--------------------|\n" );
	fprintf(outfile, "First\t|%9.6f  |\n",analyticalprime);
	fprintf(outfile, "Second\t|%9.6f  |\n",analyticaldoubleprime);
	fprintf(outfile, "Third\t|%9.6f  |\n",analyticaltripleprime);
	fprintf(outfile, "Fourth\t|%9.6f  |\n",analyticalforthprime);

	double dx_array[7] = {0.0005,0.001,0.01,0.1,0.2,0.3,0.4};

	int length = sizeof(dx_array)/sizeof(dx_array[0]);

	//array of derivatives
	double first[length];
	double second[length];
	double third[length];
	double forth[length];

	//array of errors

	double errorprime[length];
	double errordoubleprime[length];
	double errortripleprime[length];
	double errorforthprime[length];

	
	fprintf(outfile,"\n");
	fprintf(outfile,"TABLE OF NUMERICAL DERIVATIVES AT DIFFERENT dx AND %%ERROR\n\n");

//Printing to File

	fprintf(outfile, "\t\t|\tFIRST DERIVATIVE\t |\tSECOND DERIVATIVE\t  |\tTHIRD DERIVATIVE     |\tFOURTH DERIVATIVE\t |\n");
	fprintf(outfile, "----------------------------------------------------------------------------------------------------------\n");
	fprintf(outfile,"dx\t\t|\tNumerical\t %%Error\t |\tNumerical\t  %%Error  |\tNumerical\t %%Error\t |\tNumerical\t %%Error\t |\n");
	fprintf(outfile, "----------------------------------------------------------------------------------------------------------\n");


	for (int i = 0; i<length; i++){

		first[i] = prime(dx_array[i],x0);
		second[i] = doubleprime(dx_array[i],x0);
		third[i] = tripleprime(dx_array[i],x0);
		forth[i] = quadprime(dx_array[i],x0);

		errorprime[i] = 100.00*((first[i]-analyticalprime)/(analyticalprime));
		errordoubleprime[i] = 100.00*((second[i]-analyticaldoubleprime)/(analyticaldoubleprime));
		errortripleprime[i] = 100.00*((third[i]-analyticaltripleprime)/(analyticaltripleprime));
		errorforthprime[i] = 100.00*((forth[i]-analyticalforthprime)/(analyticalforthprime));

		

		fprintf(outfile,"%2.4f\t|  %2.6f\t%2.3f   |\t%2.6f\t %2.3f   |\t%2.6f\t%2.3f   |\t%2.6f\t%6.3f   |\n", dx_array[i], first[i], errorprime[i], second[i], errordoubleprime[i], third[i], errortripleprime[i],forth[i],errorforthprime[i]);		
	}

	fclose(outfile);

}