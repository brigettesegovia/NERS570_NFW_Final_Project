#include "ILU.hpp"

#include <vector>
#include <map>
#include <cassert>
#include <stdexcept>
#include <cmath>

using namespace std;

namespace NBody
{

    /**
     * @brief Constructs an ILU(0) preconditioner
     *
     * Calculates the incomplete LU factorization of A,
     * preserving the original sparsity pattern.
     *
     * @param A: input square matrix in CSR format
     */
    ILU::ILU(const Matrix& A)
        : _L(A.GetNrows(), A.GetNcols()),
          _U(A.GetNrows(), A.GetNcols()),
          _n(A.GetNrows())
    {
        // Make sure the matrix is square
        if (A.GetNrows() != A.GetNcols())
        {
            throw invalid_argument("ILU requires a square matrix.");
        }

        // Temporary sparse row storage for factorization
        // Each map stores: column index -> matrix value
        vector<map<size_t, double>> rows(_n);

        // Copy the original matrix into temporary storage
        for (size_t i = 0; i < _n; i++)
        {
            for (const auto& [j, val] : A.GetRowEntries(i))
            {
                rows[i][j] = val;
            }
        }

        // ILU(0) FACTORIZATION

        for (size_t i = 0; i < _n; i++)
        {
            // Process entries below the diagonal
            for (auto it = rows[i].begin();
                 it != rows[i].end() && it->first < i;
                 ++it)
            {
                size_t j = it->first;

                // Get the diagonal pivot U(j,j)
                auto pivot_it = rows[j].find(j);

                if (pivot_it == rows[j].end() ||
                    pivot_it->second == 0.0)
                {
                    throw runtime_error("Zero or missing pivot in ILU.");
                }

                double pivot = pivot_it->second;

                // Calculate L(i,j)
                double lij = it->second / pivot;

                // Store the lower triangular factor
                it->second = lij;

                // Update existing entries only
                // This preserves the original sparsity pattern
                for (auto upper = rows[j].upper_bound(j);
                     upper != rows[j].end();
                     ++upper)
                {
                    size_t k = upper->first;

                    auto entry = rows[i].find(k);

                    // Only update if the entry already exists
                    if (entry != rows[i].end())
                    {
                        entry->second -= lij * upper->second;
                    }
                }
            }
        }

        // EXTRACT L AND U INTO CSR MATRICES

        for (size_t i = 0; i < _n; i++)
        {
            // L has a unit diagonal
            _L.SetValue(i, i, 1.0);

            for (const auto& [j, val] : rows[i])
            {
                if (j < i)
                {
                    // Lower triangular entries
                    _L.SetValue(i, j, val);
                }
                else
                {
                    // Diagonal and upper triangular entries
                    _U.SetValue(i, j, val);
                }
            }
        }

        // Assemble both matrices into CSR format
        _L.Assemble();
        _U.Assemble();
    }


    /**
     * @brief Forward substitution
     *
     * Solves Lz = rhs
     *
     * L has a unit diagonal.
     *
     * @param rhs: right-hand-side vector
     * @return solution vector z
     */
    vector<double> ILU::ForwardSolve(
        const vector<double>& rhs) const
    {
        assert(rhs.size() == _n);

        vector<double> z(_n, 0.0);

        for (size_t i = 0; i < _n; i++)
        {
            double sum = rhs[i];

            // Only use entries below the diagonal
            for (const auto& [j, val] : _L.GetRowEntries(i))
            {
                if (j < i)
                {
                    sum -= val * z[j];
                }
            }

            // Unit diagonal, so no division required
            z[i] = sum;
        }

        return z;
    }


    /**
     * @brief Backward substitution
     *
     * Solves Uy = z
     *
     * @param rhs: right-hand-side vector z
     * @return solution vector y
     */
    vector<double> ILU::BackwardSolve(
        const vector<double>& rhs) const
    {
        assert(rhs.size() == _n);

        vector<double> y(_n, 0.0);

        // Iterate backwards
        for (int i = static_cast<int>(_n) - 1; i >= 0; i--)
        {
            double sum = rhs[i];
            double diagonal = 0.0;

            // Iterate through the upper triangular row
            for (const auto& [j, val] :
                 _U.GetRowEntries(static_cast<size_t>(i)))
            {
                if (j > static_cast<size_t>(i))
                {
                    sum -= val * y[j];
                }
                else if (j == static_cast<size_t>(i))
                {
                    diagonal = val;
                }
            }

            if (diagonal == 0.0)
            {
                throw runtime_error("Zero diagonal in ILU upper factor.");
            }

            y[i] = sum / diagonal;
        }

        return y;
    }


    /**
     * @brief Applies the ILU preconditioner
     *
     * Solves LUy = rhs using forward and backward substitution.
     *
     * @param rhs: input residual vector
     * @return preconditioned vector
     */
    vector<double> ILU::Apply(
        const vector<double>& rhs) const
    {
        // First solve Lz = rhs
        vector<double> z = ForwardSolve(rhs);

        // Then solve Uy = z
        vector<double> y = BackwardSolve(z);

        return y;
    }

}
