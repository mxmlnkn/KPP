typedef Matrix(*FuncVectorNParam)(Matrix &a);
typedef double(*FuncScalarNParam)(Matrix &a);
typedef double(*FuncScalarNParam1Var)(Matrix &a, double x);

class ChiSquareNumeric {
    private:
        static void *ActiveObject;                                              //Stores this-Pointer of current Object needed for certain functions

        int NValues;
        FuncScalarNParam1Var FitFunction;
        double *xval,*yval,*yerr;                                               //Only Pointer to the actual array (array will not be copied, therefore take care not to delete the referenced array before calling to "Calc")
        
	public:
		double CorrelCoeff,DegOfFreedom,StDm,StDn,StDFit;                       //These public variables store the results of Calc
        Matrix *Parameter;                                                      //This array instead gets created by new

        ChiSquareNumeric (FuncScalarNParam1Var f, const int NParam,             //f with NParam Parameters
                double xval[], double yval[], double yerr[], const int N);
        ~ChiSquareNumeric (void);
        inline void SetActive(void);

        //static objects need to be static because I have to give a function parameter of them to Differentiate() and Newton()
        static double ChiSquare(Matrix &a);                                     //Even though the argument is void it makes use of nParam and Parameter in this class
        static Matrix NecessaryConditions(Matrix &a);
        void Calc(void);                                                        //Calculate ChiSq-Minimum with N values
};
