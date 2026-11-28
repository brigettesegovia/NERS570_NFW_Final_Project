#ifndef __NBody_ILU__
#define __NBody_ILU__

#include "Matrix.hpp"
#include <vector>

using namespace std;

namespace NBody
{
    class ILU
    {
        private:
            // Incomplete LU factors
            Matrix _L;
            Matrix _U;

            // Matrix dimensions
            size_t _n;

        public:
            // Constructor
            ILU(const Matrix& A);

            // Destructor
            ~ILU() = default;

            // Apply preconditioner
            vector<double> Apply(const vector<double>& rhs) const;

            // Forward substitution
            vector<double> ForwardSolve(const vector<double>& rhs) const;

            // Backward substitution
            vector<double> BackwardSolve(const vector<double>& rhs) const;

            // Access the lower and upper triangular matrices
            const Matrix& GetL() const noexcept { return _L; }
            const Matrix& GetU() const noexcept { return _U; }

    };
}

#endif
