#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include "Matrices.h"

#define DEBUG 1                     //I wouldn't set it to 2, because then a flood of error messages will erupt into the console xD

using namespace std;

int main (void) {
        Matrix mat1(3,3), mat2(3,3);
        mat1.Load("Matrix1.dat");
        mat1.Print();
        cout << endl;
        mat2.Load("Matrix2.dat");
        mat2.Print();
        (mat1+mat2).Print();        //Matrix gets deleted after being printed :)
        cout << "m1 + m2:\n";
        Matrix mat3=(mat1+mat2);
        mat3.Print();
        cout << "m1 - m2:\n";
        mat3=mat2-mat1;
        mat3.Print();
        cout << "Product:\n";
        (mat1*mat2).Print(); 
        cout << "Determinante: m1: " << mat1.Det() << " m2: " << mat2.Det() << " m1-m2: " << mat3.Det() << "\n";
    
        Matrix id = mat1;
        id.SetDiagonal(1);
        cout << "Temporary Matrix to calculate Inverse of m1:\n";
        mat1.Augment( id ).Print();
        cout << "m1: Norm: " << mat1.Norm() << "  Rank: " << mat1.Rank() << "  Trace: " << mat1.Trace() << "\n";
        cout << "m1 formed into Trapeze Matrix:\n";
        mat1.RowEchelon().Print();

        cout << "Inverse of m1:\n";
        mat1.Invert().Print();
        cout << "Inverse of m2:\n";
        mat2.Invert().Print();
        cout << "Inverse of m1-m2:\n";
        mat3.Invert().Print();
        cout << "Test Inverse Matrix by Multiplicating it with the original Matrix:\n";
        cout << "m1 * m1^-1:\n";
        (mat1*mat1.Invert()).Print();
        cout << "m2^-1 * m2:\n";
        (mat2.Invert()*mat2).Print();
        cout << "(m1-m2) * (m1-m2)^-1:\n";
        (mat3*mat3.Invert()).Print();

        cout << "Transposed of m1:\n";
        mat1.Transpose().Print();
    
        cout << "Inverse of M4:\n";
        Matrix m4(2,2);
        m4.Load("m4.dat");
        cout << "Loaded\n";
        (m4.Invert()).Print();
        m4.Print();
        m4.Save("m4save.dat");
    
        cout << "Transposed of m5:\n";
        Matrix m5("m5.dat");                        //Automatically Load from file
        m5.Transpose().Print();
        cout << "see above: Rank: " << m5.Rank() << "  Trace: " << mat1.Trace() << "\n";
        m5.DelRow(2);
        cout << "Norm of "; m5.Print(); cout << " = " << m5.Norm() << " = " << m5.Transpose().Norm() << "\n";

        cout << "Identity:\n";
        Matrix m6(3,5);
        m6.SetDiagonal(3);
        m6.Print();
        cout << "Identity:\n";
        Matrix m7(7,7);
        m7.SetDiagonal();
        m7.Print();
        m7.SetAll(7);
        cout << "SetAll:\n";
        m7.Print();

        Matrix m8("m3.dat");
        cout << "Rank of: \n"; m8.Print(); cout << endl; m8.RowEchelon().Print(); cout << "is: " << m8.Rank();
        cout << "\n 0 == " << m8.RowEchelon()(1,0) << " <=> " << (m8.RowEchelon()(1,0) == 0) << "  but 0==0 <=> " << (0==0) << ". That's why in Rank I made the condition (0!=abs(element))";

        cin.get();

        return 0;
}

#include "Matrices.cpp"
