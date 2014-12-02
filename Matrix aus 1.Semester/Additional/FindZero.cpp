#ifndef PRECISION
    #define PRECISION 1e-5
#endif
#define DIVERG 1e3          //If there is a problem the interval [x-DIVERG,x+DIVERG] may have to be analyzed
#define NCALC 1e4           //Amount of values to be calculated

Matrix Nabla(FuncScalarNParam f, Matrix &x0) {
    int dim = x0.GetVectorDim();
    if (dim < 1) {
        Matrix mat = (double)0;
        return mat;
    }
    Matrix res(x0.GetDim().m, x0.GetDim().n);                                   //Column Vector
    for (int i=0; i<dim; i++) {
        res[i] = Differentiate(f, x0, i);
    }
    return res;
}

Matrix Hessian( FuncScalarNParam f, Matrix &x0) {                   //Return Hesse-Matrix evaluated at point x0 (column vector)
    // Hesse(i,j) = d(f_j)/d(x_i) so Jacobi(i-th row) = Nabla(f_i) 
    if (!x0.IsVector()) {
        #if DEBUG >=1
            cout << "Can't calculate Hesse Matrix with x0 not being a column or row vector\n";
        #endif
        Matrix Zero(1,2);                                           //Is clearly an errorcode, becose Jacobi Matrix is a square Matrix
        return Zero;
    }
    int NParam = x0.GetVectorDim();
    Matrix Hesse(NParam,NParam);
    for (int i=0; i<NParam; i++) {
        for (int j=0; j<NParam; j++) {
            int vars[] = {i,j};
            Hesse(i,j)= Differentiate( f, x0, vars, 2);
        }
    }

    return Hesse;
}

Matrix Jacobian( FuncVectorNParam f, Matrix &x0) {                  //Return Jacobi-Matrix evaluated at point x0 (column vector)
    /* Jacobi(i,j) = d(f_j)/d(x_i) so Jacobi(i-th row) = Nabla(f_i) */
    if (!x0.IsVector()) {
        #if DEBUG >=1
            cout << "Can't calculate Jacobi Matrix with x0 not being a column or row vector\n";
        #endif
        Matrix Zero(1,2);                                           //Is clearly an errorcode, becose Jacobi Matrix is a square Matrix
        return Zero;
    }
    int NParam = x0.GetVectorDim();
    Matrix Jac(NParam,NParam);
    for (int i=0; i<NParam; i++) {
        Matrix dxi = Differentiate(f, x0, i);
        //copy dxi into i-th column
        for (int j=0; j<NParam; j++) {
            Jac(j,i)=dxi[j];
        }
    }
    return Jac;
}

/* This function heavily needs a mathemtical function which has zeroes, so unless you want the program to
never end or an stack overflow to happen (because this function recursively calls to itself), test if the function
has in fact zeroes */
double Newton( double(*f)(double), double x0) {
	double xn=x0, xn1=x0, diff, temp;
	const double dx = PRECISION;
	int N=0;
	do {
        N++;
		temp = xn1;
		diff = f(xn) / ( (f(xn+dx)-f(xn)) /dx);
		xn1 = xn - diff;
		xn = temp;
        if ( (N>=1e7) )
        {
            cout << "N: " << N << " - xn1: " << xn1 << "xn: " << xn << endl;
            N=0;
            if ( f(xn1+DIVERG) * f(xn1) <= 0)
                xn1+=DIVERG/2;
            else if ( f(xn1-DIVERG) * f(xn1) <= 0)                        //thest wether the zero is to the left of x0
                xn1-=DIVERG/2;
            else
                return numeric_limits<double>::quiet_NaN();               //if there are no zeroes in [x0-DIVERG,x0+DIVERG] then return NAN
        }
        if ( (isnan(xn1)) or (isinf(xn1)) )
            break;
	} while (abs( f(xn1) ) > PRECISION);
	
    if (isnan(xn1) or isinf(xn1))
        {
            cout << "N: " << N << " - xn1: " << xn1 << endl;
    }
	return xn1;
}

/* In contrast Newton's method the bisection-method always converges. Though it may not recognize Zeroes
for example if the function is x^2-4 and a=-3 and b=3, because there is an even number of zeroes in [a,b]
It also may have problem if there are more than two zeroes in [a,b], so be carefull to test this, before
calling this function (In this case there should always be only one zero in the Kepler-function) */
double Bisection( double(*f)(double), double a, double b )
{
	double c=(a+b)/2;
    if ( (f(a) <= 0) and (f(b) >= 0) ) //test whether where are zeroes in [a,b]
    {
    	while (abs(b-a) > PRECISION)
    	{
    		c=(a+b)/2;
    		(f(c) > 0) ? (b=c) : (a=c);
    	}
	   return (a+b)/2;
    }
    else
    {
        return numeric_limits<double>::quiet_NaN();       //return NAN, to indicate that no zero could be found
    }
}

#define dx 1e-3                                                                 //TODO: If too small then the result will be 0 and if too large the result will be inprecise :(
Matrix Differentiate(FuncVectorNParam f, Matrix a, int n) {                     //n counts from 0
    MDIM dim=a.GetDim();
    if (n<0 || (dim.n==1 && n>dim.m) || (dim.m==1 && n>dim.n) || (dim.n!=1 && dim.m!=1) ) {
        #if DEBUG >=1
            cout << "Can't differentiate " << dim.m << "x" << dim.n << " " << n << "th Variable\n";
        #endif
        
        Matrix Zero(1,1);                                                       //Todo: =double Constructor for Matrix
        return Zero;
    }

    a[n] += dx; 
    Matrix res = f(a);
    a[n] -= 2*dx;
    res = (res-f(a))/(2*dx);

    return res;
}

double Differentiate(FuncScalarNParam f, Matrix a, int n) {                     //n counts from 0
    MDIM dim=a.GetDim();
    if (n<0 || (dim.n==1 && n>dim.m) || (dim.m==1 && n>dim.n) || (!a.IsVector()) ) {
        #if DEBUG >=1
            cout << "Can't differentiate " << n << "th Variable\n";
        #endif
        return 0;
    }
    a[n] += dx; 
    double res = f(a);
    a[n] -= 2*dx;
    res = (res-f(a))/(2*dx);
    return res;
}

double Differentiate(FuncScalarNParam f, Matrix x0, const int *Vars, int n) {
    if (n==1) {
        return Differentiate(f,x0,Vars[0]);
    }
    else {                                                                      //Recursive call to this function because df^2/(dx*dy) = d(df/dx)/dy
        int nth = Vars[0];
        x0[nth] += dx; 
        double res = Differentiate(f,x0,&Vars[1],n-1);
        x0[nth] -= 2*dx;
        res = (res-Differentiate(f,x0,&Vars[1],n-1))/(2*dx);
        return res;
    }
}

Matrix Newton( FuncVectorNParam f, Matrix &x0, double NewtonPrec, double NewtonCalc) {
    MDIM dim = x0.GetDim();
    int NParam = dim.m;
    Matrix xn=x0, Jac(NParam,NParam);                                           //Column Vectors

	int N=0;
	do {
        N++;
        for (int i=0; i<NParam; i++) {
            Matrix dxi = Differentiate(f, x0, i);
            //copy dxi into i-th column
            for (int j=0; j<NParam; j++) {
                Jac(j,i)=dxi[j];
            }
        }
        xn = xn - Jac.Invert() * f( xn );

        if ( (N>=NCALC) ) {
            #if DEBUG >= 1
                cout << "Multidimensional Newton didn't converge(Norm = " << f(xn).Norm() << ") after " << NewtonCalc << " steps!!!\n";
            #endif
            break;
        }

        for (int i=0; i<NParam; i++)
            if ( (isnan(xn[i])) or (isinf(xn[i])) )
                break;
	} while (f(xn).Norm() > NewtonPrec);
	return xn;
}













void SortDesc (double list[], int nelmt)                   //Sorts an array of >nelemt< values in descending order
{
    int swapped;                                            //Counter containing amount of executed swappings
    do
    {
        swapped=0;                                          //Set back counter at beginning of every cycle
        for (int i=nelmt-1; i>0; i--)
        {
            if (list[i-1] < list[i])                        //Swap elements if prior element is smaller
            {
                double swapbuf = list[i];
                list[i] = list[i-1];
                list[i-1] = swapbuf;
                swapped++;
            }
        }
    } while (swapped!=0);                                   //Exit loop if there are no "wrongly" ordered pair of elements to be found
    return;
}

double Mean (double list[], int nelmt)                      //Returns the arithmetic mean of an array with >nelmt< elements
{
    double mean=0;
    for (int i=0; i<nelmt; i++)
    {
        mean += list[i];                                    //Add all elements together and save it into "mean"
    }
    return mean/nelmt;                                      //Return the arithmetic mean (sum of all elements divided by number of elements)
}

double StdDev (double list[], int nelmt)                   //Returns standard-deviation of an array
{
    double var = 0;
    double mn = Mean(list, nelmt);
    
    for (int i=0; i<nelmt; i++)                             //Add up the the differences of each element and the arithmetic mean raised to the power 2
    {
        var += pow( list[i]-mn ,2);
    }
    
    return sqrt(var/(nelmt-1));                             //Return standard deviation, which is the squareroot of the variance


/****************************************************************
 ***  Alternative method to calculate the standard deviation  ***
 ****************************************************************
	double lspw2[nelmt];
    for (int i=0; i<nelmt; i++)
    {
        lspw2[i] = pow(list[i],2);
    }
    return sqrt((mean(lspw2,nelmt)-pow(mean(list,nelmt),2))*nelmt/(nelmt-1));
*/
}




void CntDiff (double x[], double y[], double dy[], int counter)
//Centered Difference
//because of the nature of the algorithm the FIRST and LAST element of dy will be set to 0
//counter contains the amount of element, not the maximum index of the arrays
{
	dy[0]=0;
	dy[counter-1]=0;

	for (int i=1; i<counter-1; i++)
	{
		dy[i]=( (y[i+1]-y[i-1]) / (x[i+1]-x[i-1]) );
	}
	return;
}

#undef DIVERG
#undef NCALC
