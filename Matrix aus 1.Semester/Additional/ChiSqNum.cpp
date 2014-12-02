void *ChiSquareNumeric::ActiveObject;

ChiSquareNumeric::ChiSquareNumeric (FuncScalarNParam1Var f, const int NParam,             //f with NParam Parameters
double xval[], double yval[], double yerr[], const int N) {
    this->Parameter = new Matrix(NParam,1);
    this->FitFunction = f;
    this->xval = xval;
    this->yval = yval;
    this->yerr = yerr;
    this->NValues = N;
    SetActive();
    return;
}

ChiSquareNumeric::~ChiSquareNumeric (void) {
    delete Parameter;
    return;
}

inline void ChiSquareNumeric::SetActive(void) {
    ActiveObject = this;
    return;
}

double ChiSquareNumeric::ChiSquare(Matrix &a) {
    ChiSquareNumeric *THIS = static_cast<ChiSquareNumeric*>(ChiSquareNumeric::ActiveObject);
    double sum=0,temp;
    for (int i=0; i<THIS->NValues; i++) {
        double w = 1/(THIS->yerr[i]);
        temp = (THIS->yval[i] - THIS->FitFunction(a, THIS->xval[i]));
        sum += w*w*temp*temp;
    }
    return sum;
}

Matrix ChiSquareNumeric::NecessaryConditions(Matrix &a) {                //linear system \vec{f} == 0 is sought
    return Nabla( ChiSquare,a );
}

void ChiSquareNumeric::Calc(void) { //Calculate ChiSq-Minimum
    /* Minimize Chi-Square... I know, there is a Gauss-Newton Method which would
       be better, but this is a proof of concept, that it would work with fundamental
       functions. The optimization can be done later. */
    SetActive();
    *Parameter = Newton( (FuncVectorNParam) &NecessaryConditions, *Parameter);

    DegOfFreedom = NValues-Parameter->GetVectorDim();
    //CorrelCoeff = -a12/sqrt(a11*a22);                           //C(0,1)/(D*sqrt(C(0,0)/D)*sqrt(C(1,1)/D));
    //StDn = sqrt(a22/D);
	//StDm = sqrt(a11/D);
    //StDFit = sqrt(ChiSquare(x,y,yerr,N)/DegOfFreedom);
	return;
}
