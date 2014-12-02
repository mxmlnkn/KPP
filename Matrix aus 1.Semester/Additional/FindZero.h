#ifndef MATRICES_H
    #include "Matrices.h"
#endif

typedef Matrix(*FuncVectorNParam)(Matrix &a);
typedef double(*FuncScalarNParam)(Matrix &a);

Matrix Hessian( FuncScalarNParam f, Matrix &x0);                        //Return Hesse-Matrix evaluated at point x0 (column vector)
Matrix Jacobian( FuncVectorNParam f, Matrix &x0);                       //Return Jacobi-Matrix evaluated at point x0 (column vector)

double Newton( double(*f)(double), double x0);
double Bisection( double(*f)(double), double a, double b);
Matrix Newton( FuncVectorNParam f, Matrix &x0, double PRECISION = 1e-5, double NCALC = 1e5);
double Differentiate(FuncScalarNParam f, Matrix a, int n);
Matrix Differentiate(FuncVectorNParam f, Matrix a, int n);
double Differentiate(FuncScalarNParam f, Matrix x0, const int *Vars, int n);   //Differentiate n-fold with Vars in that sequence
Matrix Nabla(FuncScalarNParam f, Matrix &x0);

void SortDesc (double list[], int nelmt);                               //Sorts an array of >nelemt< values in descending order
double Mean (double list[], int nelmt);                                 //Returns the arithmetic mean of an array with >nelmt< elements
double StdDev (double list[], int nelmt);                               //Returns standard-deviation of an array
void CntDiff (double x[], double y[], double dy[], int counter);

double IntCenter (double(*f)(double), double a, double b, int N=100 );
double IntTrapeze (double(*f)(double), double a, double b, int N=100 );
double IntSimpson (double(*f)(double), double a, double b, int N=100 );
