#include "Matrix.hpp"

#include <vector>
#include <cassert>

using namespace std;

namespace NBody
{
    
    /**
     * @brief: Constructs & initializes a 3D CSR `Matrix`
     * 
     * Requires the number of rows and columns to be set
     * 
     * @param nrows: number of rows (size_t)
     * @param ncols: number of columns (size_t) 
     */
    Matrix::Matrix(const size_t nrows, const size_t ncols) :
    _nrows(nrows), _ncols(ncols)
    {
        assert(this->_state == MatrixState::undefined);
        this->_state = MatrixState::initialized;
        assert(this->_state == MatrixState::initialized);
    }

    /**
     * @brief: Sets a nonzero value in the CSR `Matrix`
     * 
     * @param i: row index (size_t)
     * @param j: column index (size_t)
     * @param value: nonzero value (double)
     */
    void Matrix::SetValue(const size_t i, const size_t j, double val)
    {
        assert(this->_state >= MatrixState::initialized);
        assert(i < this->_nrows);
        assert(j < this->_ncols);
        
        //disassemble matrix if needed
        if (this->_state == MatrixState::assembled) { this->Disassemble(); }
        // otherwise set matrix state to building!
        else { this->_state = MatrixState::building; }

        // store value
        this->_buildCoeff[make_pair(i,j)] = val;
        this->_nnz = this->_buildCoeff.size();

        assert(this->_state == MatrixState::building);
    }

    /**
     * @brief Returns a value at the requested i and j 
     * 
     * @param i: row index (size_t)
     * @param j: column index (size_t)
     * @return value (double)
     */
    double Matrix::GetValue(const size_t i, const size_t j) const noexcept
    {
        assert(this->_state > MatrixState::initialized);
        assert(i < this->_nrows);
        assert(j < this->_ncols);

        // get value from buildCoeff if still building
        if (this->_state == MatrixState::building) 
        {
            // search for index pair in _buildCoeff
            auto it = this->_buildCoeff.find(std::make_pair(i, j));
            // If the entry is found, return its value; otherwise, return 0
            if (it != this->_buildCoeff.end())
            {
                return it->second;
            }
            return 0.0;
        }
        
        // otherwise, extract it from arrays
        else
        {
            // first, check if value is a nonzero (i.e. in the value array)
            // loop through nonzeros in row
            for (size_t k=0; k < _row_ptr[i+1]-_row_ptr[i]; k++)
            {
                // if the column index = j, return the corresponding value
                if (_col_idx[_row_ptr[i]+k]==j) { return _values[_row_ptr[i]+k]; }
            }
            // otherwise return zero
            return double{0.0};
        }
    }

    /**
     * @brief Assembles the CSR `Matrix`
     * 
     * Converts the `_buildCoeff` into the storage arrays
     */
    void Matrix::Assemble()
    {
        assert(this->_state <= MatrixState::building);

        // temporary working vectors
        vector<size_t> temp_row_idx(this->_nnz); // stores row indices of elements
        vector<size_t> temp_col_idx(this->_nnz); // stores column indices of elements
        vector<double> temp_val(this->_nnz);     // stores values of elements
        //vector<size_t> idxs(this->_nnz);         // stores indices to sort other arrays by later

        // allocate memory for row_ptr - one more element than nrows!
        this->_row_ptr.resize(this->_nrows+1);

        // first convert into COO and count row sizes
        size_t entry_idx = 0;
        for (auto &[key, val] : this->_buildCoeff)
        {
            this->_row_ptr[key.first+1]++;
            temp_row_idx[entry_idx] = key.first;
            temp_col_idx[entry_idx] = key.second;
            temp_val[entry_idx] = val;
            entry_idx++;
        }

        // no longer need _buildCoeffs! deallocate their memory
        this->_buildCoeff.clear();

        // fill values and col_idxs 
        this->_col_idx.resize(this->_nnz);
        this->_values.resize(this->_nnz);
        for (size_t k=0; k<this->_nnz; k++)
        {
            this->_col_idx[k] = temp_col_idx[k];
            this->_values[k] = temp_val[k];
        }
        // convert row_ptr to actual row_ptr array by doing a prefix sum
        for (size_t k=1; k<=this->_nrows; k++)
        {
            this->_row_ptr[k] = this->_row_ptr[k] + this->_row_ptr[k-1];
        }
        
        this->_state = MatrixState::assembled;
    }

    /**
     * @brief Disassembles the CSR `Matrix`
     * 
     * Converts the storage arrays back into `_buildCoeff`s
     */    
    void Matrix::Disassemble()
    {
        assert(this->_state==MatrixState::assembled);

        // rebuild _buildCoeffs
        // loop through rows
        for (size_t i=0; i<this->_nrows; i++)
        {
            // figure out how many elements are in this row
            size_t rowsize = _row_ptr[i+1] - _row_ptr[i];

            // loop through all nonzeros in row
            for (size_t k=0; k<rowsize; k++)
            {
                // row index is i
                // column index:
                size_t j = _col_idx[_row_ptr[i]+k];
                this->_buildCoeff[make_pair(i,j)] = _values[_row_ptr[i]+k];
            }
        }
        // clear storage arrays
        _values.clear(); _col_idx.clear(); _row_ptr.clear();

        this->_state = MatrixState::building;
    }

    /**
     * @brief multiplies a CSR `Matrix` with an input vector
     * 
     * @param vecin: input vector (vector<double>)
     * @returns: product of the matrix and `vecin` (vector<double>)
     */
    vector<double> Matrix::MatVec(const vector<double> vecin)
    {
        // make sure matrix is ready to go
        assert(this->_state == MatrixState::assembled);
        // make sure our input vector is the correct size
        // i.e. its length = ncols
        assert(vecin.size() == this->_ncols);

        // output vector - initialize to be all zeros
        // in case a row of the matrix turns out to be all zeros
        // and thus gets skipped over during the matvec routine
        vector<double> vecout(this->_nrows, 0);

        // loop through rows
        for (size_t i=0; i<this->_nrows; i++)
        {
            // check how many in this row
            size_t rowsize = _row_ptr[i+1]-_row_ptr[i];

            // loop through all nonzeros in row
            for (size_t k=0; k<rowsize; k++)
            {
                // row index is i
                // column index:
                size_t j = _col_idx[_row_ptr[i]+k];
                // do the multiplication
                vecout[i] = vecout[i] + _values[_row_ptr[i]+k]*vecin[j];
            }
        }

        return vecout;
    }

    
}