#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

struct ChiSqLnRes {
    double m,n,CorrelCoeff,DegOfFreedom,StDm,StDn,StDFit,ChiSquare;
};
ChiSqLnRes LinearFit(double xval[], double yval[], double yerr[], const int NValues);

/*******************************************
***************** MAIN *********************
*******************************************/

int main(void)
{
	#define NELEMENTS 100
	double xval[NELEMENTS],yval[NELEMENTS],yerr[NELEMENTS];
    int nelmt=0;

	/******* Read from File *******/
	ifstream fa("Pa234.dat");
    if (fa.fail()) {
       cout << "Couldn't open file\n";
       return 1;
    }

    do {
        //Comments (beginning with '#') shall be ignored
		if ((char)fa.peek()=='#') {
			fa.ignore(80,'\n');
			continue;
		}

        fa >> xval[nelmt] >> yval[nelmt];
        if (fa.fail()) {
            cout << "Couldn't read a line from file\n";     //It's the white line at the end of the file which triggers this messages
            continue;
        }

        /* Calculate Error */
        yerr[nelmt] = sqrt(yval[nelmt]);                    //Because the yval[i]s are Poisson distributed and yval[nelmt] the only value in this distribution and therefore also the mean and therefore also the variance of the distribution

        /* MAKE LINEAR (the error must be linearized too) */
        yerr[nelmt] = abs(yerr[nelmt]/yval[nelmt]);         //NEW ERROR CALCULATED WITH ERROR PROPAGATION !!! f(x)=ln(x+-sx) -> sf = sqrt( (1/x *sx)^2 ) = |sx/x|
        yval[nelmt] = log(yval[nelmt]);

        nelmt++;
    } while ((!fa.eof()) and (nelmt<NELEMENTS));            //Stop reading if eof reached or the end of array is reached
    /**** nelmts contains number of elements in array ****/
    /******* Read from File *******/
    #undef NELEMENTS
	ChiSqLnRes LinFunc = LinearFit( xval,yval,yerr,nelmt );
	cout << "--- WEIGHTED ---"
         << "\n\n(Mean +/- 1*Sigma)"
         << "\nm:" << LinFunc.m << " +/- " << LinFunc.StDm
         << "\nn:" << LinFunc.n << " +/- " << LinFunc.StDn
         << "\n\nMinimum ChiSquare: " << LinFunc.ChiSquare
         << "\nNumber of Degrees of Freedom: " << LinFunc.DegOfFreedom
         << "\nCorrelation Coefficient: " << LinFunc.CorrelCoeff
         << "\nStandard Deviation of Fit: " << LinFunc.StDFit
         << "\nHalf Life (T12): " << -1/LinFunc.m*log(2) << "   N0: " << exp(LinFunc.n)               //log( N0*exp(-t/tau) ) = log(N0)-t/tau = log(N0)-t*ln(2)/T12
         << "\n\n => N(t)=" << exp(LinFunc.n) << "*2^(-t/" << -1/LinFunc.m*log(2) << ")";

	cin.get();
	return 0;
}


ChiSqLnRes LinearFit(double x[], double y[], double yerr[], const int N) {
    ChiSqLnRes data;
	double a11=0, a12=0, a22=0, b1=0, b2=0;
	for (int i=0; i<N; i++) {
		double w = 1/(yerr[i]*yerr[i]);
		a11 += w;
		a12 += w*x[i];
		a22 += w*x[i]*x[i];
		b1 += w*y[i];
		b2 += w*x[i]*y[i];
	}

	const double D = a11*a22 - a12*a12;
    data.n = (b1*a22-b2*a12)/D;
    data.m = (b2*a11-b1*a12)/D;
    data.CorrelCoeff = -a12/sqrt(a11*a22);
    data.DegOfFreedom = N-2;
    data.StDn = sqrt(a22/D);
	data.StDm = sqrt(a11/D);

    double sum=0, temp;
    for (int i=0; i<N; i++) {
        double w = 1/(yerr[i]);
        temp = ( y[i] - (data.m*x[i]+data.n) );
        sum += w*w*temp*temp;
    }
    data.ChiSquare = sum;
    data.StDFit = sqrt(data.ChiSquare/data.DegOfFreedom);
	return data;
}
