#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "Matrices.h"
#include "FindZero.h"
#include "ChiSqNum.h"

#define DEBUG 1

using namespace std;

double FitFunc(Matrix &a, double x) { return a[0]*exp( -x*log(2)/a[1] ); }
double LinearFunc(Matrix &a, double x) { return a[0]*x+a[1]; }

int main(void)
{
	#define NELEMENTS 100
	double xval[NELEMENTS],yval[NELEMENTS],yerr[NELEMENTS];
    int nelmt=0;
	/******* Read from File *******/
	ifstream fa("Pa234.dat");
    if (fa.fail()) {
       cout << "Couldn't open file\n";                      //Triggered by white lines like the last line
       return 1;                                            //Exit program with error-code 1
    }
    do {                                                    //Read values from file into array
        //Comments (beginning with '#') shall be ignored and if so
		if ((char)fa.peek()=='#') {
			fa.ignore(80,'\n');
			continue;
		}

        fa >> xval[nelmt] >> yval[nelmt];
        if (fa.fail()) {
			cout << "Couldn't read a line from file\n";
            continue;
        }

        yerr[nelmt] = sqrt(yval[nelmt]);                    //Because the yval[i]s are Possion distributed and the single value in this distribution and therefore also the mean and therefore also the variance

        nelmt++;
    } while ((!fa.eof()) and (nelmt<NELEMENTS));            //Stop reading if eof reached or the end of array is reached
    /**** nelmts contains number of elements in array ****/
    /******* Read from File *******/
    #undef NELEMENTS

	ChiSquareNumeric NonLinFit(FitFunc,2, xval,yval,yerr,nelmt );
    //TODO: BIGGEST PROBLEM IS GUESSING START VALUES !!!
    (*NonLinFit.Parameter)[0]=1000;
    (*NonLinFit.Parameter)[1]=80;
	NonLinFit.Calc();
	cout << "--- WEIGHTED NONLINEARIZED (CORRECT FUNCTION (MULTIDIMENSIONAL NEWTON) - GRAPH 1 IN PLOT.GP) ---"
         << "\nN0:" << (*NonLinFit.Parameter)[0] << " T12:" << (*NonLinFit.Parameter)[1]
         << "\nChiSquare: " << NonLinFit.ChiSquare(*NonLinFit.Parameter)
         << "\nNumber of Degrees of Freedom: " << NonLinFit.DegOfFreedom
         /*<< "\nStandard Deviation of m: " << NonLinFit.StDm
         << "\nStandard Deviation of n: " << NonLinFit.StDn
         << "\nCorrelation Coefficient: " << NonLinFit.CorrelCoeff
         << "\nStandard Deviation of Fit: " << NonLinFit.StDFit*/
         << "\n => N(t)=" << (*NonLinFit.Parameter)[0] << "*exp(t/" << (*NonLinFit.Parameter)[1] << ")\n\n";

	ChiSquareNumeric fit(LinearFunc,2, xval,yval,yerr,nelmt );
    fit.Parameter->SetAll(0);
	fit.Calc();
	cout << "--- WEIGHTED NONLINEARIZED (WRONG FUNCTION - GRAPH 1 IN PLOT.GP) ---"
         << "\nm:" << (*fit.Parameter)[0] << " n:" << (*fit.Parameter)[1]
         << "\nChiSquare: " << fit.ChiSquare(*fit.Parameter)
         << "\nNumber of Degrees of Freedom: " << fit.DegOfFreedom
         /*<< "\nStandard Deviation of m: " << fit.StDm
         << "\nStandard Deviation of n: " << fit.StDn
         << "\nCorrelation Coefficient: " << fit.CorrelCoeff
         << "\nStandard Deviation of Fit: " << fit.StDFit*/
         << "\n => Log(N(t))=" << (*fit.Parameter)[0] << "*x+" << (*fit.Parameter)[1] << "\n\n";

	for (int i=0; i<nelmt; i++) {
        yerr[i] = 1;
    }
    fit.Parameter->SetAll(0);
	fit.Calc();
	cout << "--- UNWEIGHTED NONLINEARIZED (WRONG FUNCTION - GRAPH 1 IN PLOT.GP) ---"
         << "\nm:" << (*fit.Parameter)[0] << " n:" << (*fit.Parameter)[1]
         << "\nChiSquare: " << fit.ChiSquare(*fit.Parameter)
         << "\nNumber of Degrees of Freedom: " << fit.DegOfFreedom
        /* << "\nStandard Deviation of m: " << fit.StDm
         << "\nStandard Deviation of n: " << fit.StDn
         << "\nCorrelation Coefficient: " << fit.CorrelCoeff
         << "\nStandard Deviation of Fit: " << fit.StDFit*/
         << "\n => Log(N(t))=" << (*fit.Parameter)[0] << "*x+" << (*fit.Parameter)[1] << "\n\n";


	for (int i=0; i<nelmt; i++) {
        yerr[i] = sqrt(yval[i]);                    //Because the yval[i]s are Possion distributed and the single value in this distribution and therefore also the mean and therefore also the variance
        yerr[i] = yerr[i]/yval[i];              //NEW ERROR CALCULATED WITH ERROR PROPAGATION f(x)=ln(x+-sx) -> sf = sqrt( (1/x *sx)^2 ) = sx/x
        yval[i] = log(yval[i]);
    }
    fit.Parameter->SetAll(0);
	fit.Calc();
	cout << "--- WEIGHTED LINEARIZED (GRAPH 2 IN TESTPLOT.GP AND ACTUAL TASK)---"
         << "\nm:" << (*fit.Parameter)[0] << " n:" << (*fit.Parameter)[1]
         << "  N0=exp(n):" << exp((*fit.Parameter)[1]) << "  T12=1/m*ln(2):" << -1/(*fit.Parameter)[0]*log(2)
         << "\nChiSquare: " << fit.ChiSquare(*fit.Parameter)
         << "\nNumber of Degrees of Freedom: " << fit.DegOfFreedom
         /*<< "\nStandard Deviation of m: " << fit.StDm
         << "\nStandard Deviation of n: " << fit.StDn
         << "\nCorrelation Coefficient: " << fit.CorrelCoeff
         << "\nStandard Deviation of Fit: " << fit.StDFit*/
         << "\n => N(t)=" << (*fit.Parameter)[0] << "*x+" << (*fit.Parameter)[1] << "\n\n"
         << "\n => N(t)=" << exp((*fit.Parameter)[1]) << "*exp(-t/" << -1/(*fit.Parameter)[0]*log(2) << ")\n\n";

	cin.get();
	return 0;
}


#include "FindZero.cpp"
#include "Matrices.cpp"
#include "ChiSqNum.cpp"
