
#include "Matrix.hpp"
#include "ILU.hpp"

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

using namespace std;
using namespace NBody;

int main()
{
    // Create the matrix

    Matrix A(3, 3);

    A.SetValue(0, 0,  4.0);
    A.SetValue(0, 1, -1.0);

    A.SetValue(1, 0, -1.0);
    A.SetValue(1, 1,  4.0);
    A.SetValue(1, 2, -1.0);

    A.SetValue(2, 1, -1.0);
    A.SetValue(2, 2,  4.0);

    A.Assemble();

    // Calculate ILU factorization

    ILU ilu(A);

    // Display L

    cout << "\nL Matrix:\n";

    const Matrix& L = ilu.GetL();

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            cout << L.GetValue(i, j) << "\t";
        }

        cout << endl;
    }

    // Display U

    cout << "\nU Matrix:\n";

    const Matrix& U = ilu.GetU();

    for (size_t i = 0; i < 3; i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            cout << U.GetValue(i, j) << "\t";
        }

        cout << endl;
    }

    // Apply ILU preconditioner

    vector<double> b = {1.0, 2.0, 3.0};

    vector<double> x = ilu.Apply(b);

    cout << "\nPreconditioned solution:\n";

    for (double value : x)
    {
        cout << value << "\t";
    }

    cout << endl;

    // Check solution

    // vector<double> expected = {0.5, 1.0, 1.0};

    vector<double> expected = {
    13.0 / 28.0,
    6.0 / 7.0,
    27.0 / 28.0
    };

    double tolerance = 1e-10;

    for (size_t i = 0; i < 3; i++)
    {
        assert(abs(x[i] - expected[i]) < tolerance);
    }

    cout << "\nILU test passed!" << endl;

    return 0;
}
