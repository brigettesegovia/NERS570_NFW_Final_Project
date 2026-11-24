#ifndef __NBodySerial_Matrix__
#define __NBodySerial_Matrix__

#include <stddef.h>
#include <vector>
#include <map>          
#include <utility> 

using namespace std;

namespace NBody
{
    // stores sparse matrices in CSR format for quicker & more efficient GMRES
    // adapted from the Matrix class in the NERS570 SpMV repo (HW3),
    // my group's JDS class (also from HW3),
    // and my own CSR code from HW2
    class Matrix
    {
        // keeps track of the state of the matrix
        enum MatrixState {undefined, initialized, building, assembled};

        // keeps track of nonzero values as we add them in the building state
        using key_type = pair<size_t, size_t>; // defines pairs of size_t as a valid key type (like a Python tuple?)
        using CoeffMap = map<key_type, double>; // behaves like a Python dictionary


        private:
            MatrixState _state = MatrixState::undefined; // current state of the matrix
            size_t _nrows = 0;                           // number of rows
            size_t _ncols = 0;                           // number of columns
            size_t _nnz = 0;                             // total number of nonzero elements
            CoeffMap _buildCoeff;                        // stores nonzero values pre-assembly
            vector<double> _values;                       // stores nonzero values post-assembly
            vector<size_t> _row_ptr;                      // indices of the first nonzero in each row in values
            vector<size_t> _col_idx;                      // column indices of nonzeros
            

        public:
            // constructor
            Matrix(const size_t nrows, const size_t ncols);
            // destructor
            ~Matrix() = default;

            // getters
            [[nodiscard]] MatrixState     GetState()    const noexcept { return _state; }
            [[nodiscard]] size_t          GetNrows()    const noexcept { return _nrows; }
            [[nodiscard]] size_t          GetNcols()    const noexcept { return _ncols; }
            [[nodiscard]] size_t          GetNnz()      const noexcept { return _nnz; }
            [[nodiscard]] const CoeffMap& GetCoeffMap() const noexcept { return _buildCoeff; }
            [[nodiscard]] double          GetValue(const size_t i, const size_t j) const noexcept;

            // setter
            void SetValue(const size_t i, const size_t j, double val);

            // assembly
            void Assemble();
            void Disassemble();
            
            // linear algebra functions
            vector<double> MatVec(const vector<double> vecin) const;
            //Matrix Transpose() const;
    };
}

#endif