#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "Matrices.h"
#include "FindZero.h"

#define DEBUG 1

using namespace std;

double cf1(Matrix &a) {                                                         //reference circumvents perfomance heavy copying of whole Matrix on every call
    return pow(a[0],3) + pow(a[1],3) - 2*a[0]*a[0] + 5*a[0]*a[1]*a[1];
}

Matrix cf2(Matrix &a) {
    return Nabla(cf1,a);
}

//f[x_, y_] := x^2*y^2*Log[Abs[x + y] + 0.5] - 2 - x^2 - y^2;
double cf21(Matrix &a) { return pow(a[0],2)*pow(a[1],2) *log(abs(a[0] + a[1]) + 0.5) - 2 - pow(a[0],2) - pow(a[1],2); }
Matrix cf22(Matrix &a) { return Nabla(cf21,a);}

//Log[Abs[x + y + 3] + 0.5] + 1/(Abs[x*y^2] + 0.3) + (x^2 + y^2)/100;
double cf31(Matrix &a) { return log(abs(a[0] + a[1]+3) + 0.5) + 1/(abs(a[0]*pow(a[1],2))+0.3) + (pow(a[0],2)+pow(a[1],2))/100 ;}
Matrix cf32(Matrix &a) { return Nabla(cf31,a);}

//((x - a)^2 + y^2) Exp[1 - x^2 - y^2];
double cf41(Matrix &a) { return (pow(a[0]-a[2],2)+pow(a[1],2) )*exp( 1- pow(a[0],2) - pow(a[1],2) ); }
Matrix cf42(Matrix &a) { return Nabla(cf41,a);}


void AnalyzeCriticalPoint(FuncScalarNParam f, Matrix &xc) {
    if (Hessian(f, xc).Det() > 0) {
        cout << "Stationary Point: ";
        if (Hessian(f, xc).Trace() > 0)
            cout << "Minimum\n";
        else if (Hessian(f, xc).Trace() < 0)
            cout << "Maximum\n";
        else
            cout << "Dunno :S\n";
    } else if (Hessian(f, xc).Det() < 0) {
        cout << "Inflection Point\n";
    } else {
        cout << "No statement about the critical point possible (complex critical point like Minimum and Maximum at the same time)\n";
    }
    return;
}

int main(void) {
    Matrix x0(2,1);
    x0(0,0)=100; x0(1,0)=10;
    cout << "f(x,y) = x^3+x^3-2*x^2+5*x*y^2\n";
    cout << "x0^-1: "; x0.Print();
    cout << "\ncf1 differentiated: ";
    cf2(x0).Transpose().Print();
    Nabla(&cf1,x0).Transpose().Print();
    cout << "d/dx: ";
    Differentiate(cf2, x0, 0).Transpose().Print();
    cout << "d/dy: ";
    Differentiate(cf2, x0, 1).Transpose().Print();
    cout << "Critical point at xc: ";
    Matrix xc = Newton(cf2, x0);
    xc.Transpose().Print();
    cout << "Jacobian Matrix at xc:\n"; Jacobian(cf2, xc).Print();
    cout << "Hessian Matrix at xc:\n";
    int var[] = {0,0};
    cout << "f(x,y)/dx^2 at xc = " << Differentiate(cf1,xc,var,2) << "\n";
    var[0]=0; var[1]=1;
    cout << "f(x,y)/(dx*dy) at xc = " << Differentiate(cf1,xc,var,2) << "\n";
    var[0]=1; var[1]=0;
    cout << "f(x,y)/(dy*dx) at xc = " << Differentiate(cf1,xc,var,2) << "\n";
    var[0]=1; var[1]=1;
    cout << "f(x,y)/dy^2 at xc = " << Differentiate(cf1,xc,var,2) << "\n";
    var[0]=0;
    cout << "f(x,y)/dx at xc = " << Differentiate(cf1,xc,var,1) << " = " << Differentiate(cf1,xc,0) << "\n";
    var[0]=1;
    cout << "f(x,y)/dy at xc = " << Differentiate(cf1,xc,var,1) << " = " << Differentiate(cf1,xc,1) << "\n";
    Hessian(cf1, xc).Print();
    AnalyzeCriticalPoint(cf1,xc);

    x0(0,0)=-10; x0(1,0)=-10;
    xc = Newton(cf22, x0);
    cout << "\nNewton found xc: "; xc.Transpose().Print();
    cout << "Nabla at xc: "; Nabla(cf21,xc).Transpose().Print();
    cout << "f2(xc)=" << cf21(xc) << "\n";
    cout << "Hessian at xc: \n"; Hessian(cf21, xc).Print();
    AnalyzeCriticalPoint(cf21,xc);

    x0(0,0)=1; x0(1,0)=-1;
    xc = Newton(cf32, x0);
    cout << "\nNewton found xc: "; xc.Transpose().Print();
    cout << "Nabla at xc: "; Nabla(cf31,xc).Transpose().Print();
    cout << "f3(xc)=" << cf31(xc) << "\n";
    cout << "Hessian at xc: \n"; Hessian(cf31, xc).Print();
    AnalyzeCriticalPoint(cf31,xc);

    Matrix x04(3,1);
    x04(0,0)=1; x04(1,0)=-0.1; x04(2,0)=1;
    xc = Newton(cf42, x04);
    cout << "\nNewton found xc: "; xc.Transpose().Print();
    cout << "Nabla at xc: "; Nabla(cf41,xc).Transpose().Print();
    cout << "f4(xc)=" << cf41(xc) << "\n";
    cout << "Hessian at xc: \n"; Hessian(cf41, xc).Print();
    AnalyzeCriticalPoint(cf41,xc);

    cin.get();
    return 0;
}


#include "FindZero.cpp"
#include "Matrices.cpp"
