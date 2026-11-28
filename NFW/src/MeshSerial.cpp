#include "MeshSerial.hpp"

#include <stddef.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <tuple>
#include <algorithm>

// NOTE FROM HANGA: I MAY HAVE TO REWRITE MY STUFF WITH PROPER ARRAYS
// INSTEAD OF VECTORS OF VECTORS BECAUSE APPARENTLY THOSE ARE SLOW AF

using namespace std;

namespace NBody
{
    using namespace std;

    // constructors

    /**
     * @brief Default constructor for MeshSerial
     * 
     * Initializes a serial mesh.
     * 
     * Defaulted in MeshSerial.hpp
     * 
     */

    /**
     * @brief Constructor for MeshSerial with mesh sizes specifed
     * 
     * Also constructs the 7-point stencil `Matrix` for computing the Laplacian
     */
    MeshSerial::MeshSerial(size_t nx, size_t ny, size_t nz, double dx, double dy, double dz):
    _nx(nx), _ny(ny), _nz(nz), _dx(dx), _dy(dy), _dz(dz)
    {
        // construct 7-point stencil
        this->_stencil = _ConstructStencil();
        // initialize density as "point mass" (>0 in center cell, 0 everywhere else)
        // THIS IS FOR TESTING ONLY! REMOVE LATER!
        _density.resize(this->_nx*this->_ny*this->_nz, 0.0);
        size_t even_x = this->_nx%2; size_t even_y = this->_ny%2; size_t even_z = this->_nz%2;
        vector<size_t> cx(even_x+1,0); vector<size_t> cy(even_y+1,0); vector<size_t> cz(even_z+1,0);
        cx[0] = this->_nx / 2; if (not even_x) cx[1] = cx[0] + 1;
        cy[0] = this->_ny / 2; if (not even_y) cy[1] = cy[0] + 1;
        cz[0] = this->_nz / 2; if (not even_z) cz[1] = cz[0] + 1;
        double mtot = 1.0;
        double dens_div = double(cx.size() * cy.size() * cz.size());
        for (size_t i=0; i<cx.size(); i++)
        {
            for (size_t j=0; j<cy.size(); j++)
            {
                for (size_t k=0; k<cz.size(); k++) _density[Calc1DIndex(cx[i],cy[j],cz[k])] = mtot/dens_div;
            }
        }


    }

    /**
     * @brief Copy constructor for MeshSerial
     * 
     * Creates a copy of an existing serial mesh.
     * 
     * @param other: existing MeshSerial object to copy
     * 
     * Defaulted in MeshSerial.hpp
     */
    
    /**
     * @brief Copy assignment operator for MeshSerial
     *
     * Copies the contents of an existing serial mesh into this mesh.
     *
     * @param other: existing MeshSerial object to copy
     * @return MeshSerial&: reference to the updated mesh
     *
     * Defaulted in MeshSerial.hpp
     */

    /**
     * @brief Move constructor for MeshSerial
     *
     * Transfers the contents of an existing serial mesh into a new mesh.
     *
     * @param other: existing MeshSerial object to move
     *
     * Defaulted in MeshSerial.hpp
     */

    /**
     * @brief Move assignment operator for MeshSerial
     *
     * Transfers the contents of an existing serial mesh into this mesh.
     *
     * @param other: existing MeshSerial object to move
     * @return MeshSerial&: reference to the updated mesh
     *
     * Defaulted in MeshSerial.hpp
     */

    // destructor

    /**
     * @brief Destructor for MeshSerial
     *
     * Destroys the serial mesh.
     *
     * Defaulted in MeshSerial.hpp
     */

    // getters

    /**
     * @brief Get the gravitational potential matrix
     *
     * @return const Matrix&: reference to the potential matrix
     */

    /**
     * @brief Get the density matrix
     *
     * @return const Matrix&: reference to the density matrix
     */

    // methods

    /**
     * @brief Constructs a 7-point stencil for computing the Laplacian 
     * of the potential. This is a private method meant to be used internally
     * when the serial `Mesh` is initialized.
     */
    Matrix MeshSerial::_ConstructStencil()
    {
        // initialize stencil
        size_t stencil_size = this->_nx * this->_ny * this->_nz;
        Matrix stencil_temp (stencil_size, stencil_size);
        // calculate 1/(dx*dy*dz) (i.e. 1/volume of one cell) factor
        // needed to calculate Laplacian
        double inv_cell_vol = 1.0 / (this->_dx*this->_dy*this->_dz);

        // loop through every i,j,k position in mesh
        for (size_t i=0; i<this->_nx; i++)
        {
            for (size_t j=0; j<this->_ny; j++)
            {
                for (size_t k=0; k<this->_nz; k++)
                {
                    // convert from 3D to 1D index
                    size_t m = MeshSerial::Calc1DIndex(i, j, k);
                    // set diagonal 
                    // (subtract -6 times the value of potential at i,j,k)
                    stencil_temp.SetValue(m, m, -6.0*inv_cell_vol);
                    // set various off-diagonals 
                    // (add the values of the potential at +/-1 of i,j,k in all 3 dimensions)
                    // have to account for phi=0 boundary conditions
                    if (i<this->_nx-1) stencil_temp.SetValue(m, MeshSerial::Calc1DIndex(i+1,j,k), inv_cell_vol);
                    if (i>0)           stencil_temp.SetValue(m, MeshSerial::Calc1DIndex(i-1,j,k), inv_cell_vol);
                    if (j<this->_ny-1) stencil_temp.SetValue(m, MeshSerial::Calc1DIndex(i,j+1,k), inv_cell_vol);
                    if (j>0)           stencil_temp.SetValue(m, MeshSerial::Calc1DIndex(i,j-1,k), inv_cell_vol);
                    if (k<this->_nz-1) stencil_temp.SetValue(m, MeshSerial::Calc1DIndex(i,j,k+1), inv_cell_vol);
                    if (k>0)           stencil_temp.SetValue(m, MeshSerial::Calc1DIndex(i,j,k-1), inv_cell_vol);
                }
            }
        }
        // finally assemble and return stencil
        stencil_temp.Assemble();
        return stencil_temp;
    }

    /**
     * @brief Constructs the transpose of the Givens rotation matrix.
     * Intended for use in the least squares problem.
     * 
     * @param col: (vector<double>) the vector representing the column being modified
     * @param i: (size_t) the (row) index of the column element being used to zero out 
     *           another column element
     * @param j: (size_t) the (row) index of the column element being zeroed out
     * 
     * @returns (Matrix) the transpose of the Givens rotation matrix
     */
    NBody::Matrix MeshSerial::ConstructGivens(vector<double> col, size_t i, size_t k)
    {
        // the (transpose of the) Givens rotation matrix is EXTREMELY sparse so this warrants a CSR Matrix
        // at least it will be sparse for any reasonably large n
        NBody::Matrix GT(col.size(), col.size());

        // calculate c & s coefficicents used to zero out j-th element of column (v vector) 
        // with i-th element of column 
        double c =  col[i] / sqrt(col[i]*col[i] + col[k]*col[k]);
        double s = -col[k] / sqrt(col[i]*col[i] + col[k]*col[k]);

        for (size_t t=0; t<col.size(); t++)
        {
            if (t==i) 
            {
                GT.SetValue(i, i, c);
                GT.SetValue(i, k, -s);
            }
            else if (t==k)
            {
                GT.SetValue(k, k, c);
                GT.SetValue(k, i, s);
            }
            else GT.SetValue(t, t, 1.0);
        }

        GT.Assemble();
        return GT;
    }

    /**
     * @brief computes the least squares solution to argmin||beta*e1 - Hj*ym||_2
     * 
     * i.e. finds the ym that minimizes the residual of beta*e1 - Hj*ym using
     * QR decomposition. Q and R are constructed using Givens rotations.
     * 
     * @param H: the (m+1)x(m) Hessenberg matrix generated by GMRES
     * @param beta: the 2-norm of the initial guess potential x0's residual
     * 
     * @returns (tuple<vector<double>, double>) a tuple containing the least squares 
     *          solution vector ym and the 2-norm of its residual
     */
    tuple<vector<double>, double> MeshSerial::LeastSquares(vector<vector<double>> H, const double beta)
    {
        cout << "doing least squares now...\n" << endl;
        size_t m = H[0].size();
        vector<double> ym(m, 0.0);
        // solving H*ym = beta*e1 with least squares 
        // QR decomposition of H using Givens rotations
        // prepare Q and R matrices
        // although we actually need transpose of Q
        // H = (m+1)x(m) -> G = (m+1)x(m+1) -> R = (m+1)x(m), Q = (m+1)x(m+1)
        // most efficient way to do this is probably just to do first iteration of loop
        // otherwise would have to initialize QT as an identity matrix which probably takes way too much time
        
        // since H is an upper Hessenberg matrix, we only need
        // to zero out the one non-zero subdiagonal

        // get transpose of Givens matrix
        // v = H[:,j]
        // cout << "j=0\n" << endl;
        vector<double> v(m+1);
        for (size_t t=0; t<v.size(); t++) v[t] = H[t][0];
        Matrix GT = ConstructGivens(v, 0, 1);

        // R = product of all the GTs with H (R = GTn*GTn-1*...GT1*H)
        vector<vector<double>> R = GT.MatMul(H);
        // Q = G1*G2*...*Gn -> QT = GTn*GTn-1*...*GT1
        Matrix QT = GT;
        // cout << QT.GetState() << "\n";

        // if (m==2)
        // {
        //     cout << "R=\n";
        //     for (size_t i=0; i<m+1; i++)
        //     {
        //         for (size_t j=0; j<m; j++) cout << R[i][j] << ' ';
        //         cout << "\n";
        //     }
        //     cout << "Q=\n";
        //     for (size_t i=0; i<m+1; i++)
        //     {
        //         for (size_t j=0; j<m+1; j++) cout << QT.GetValue(j,i) << ' ';
        //         cout << "\n";
        //     }
        // }

        // do remainder of loop
        for (size_t j=1; j<m; j++)
        {
            // cout << "j=" << j << "\n";
            // get transpose of Givens matrix
            for (size_t t=0; t<v.size(); t++) v[t] = R[t][j];
            Matrix GT = ConstructGivens(v, j, j+1);
            // cout << "givens constructed\n" << endl;
            // if (m==2)
            // {
            //     cout << "GT=\n";
            //     for (size_t i=0; i<m+1; i++)
            //     {
            //         for (size_t j=0; j<m+1; j++) cout << GT.GetValue(i,j) << ' ';
            //         cout << "\n";
            //     }
            // }

            // R = GTn*GTn-1*...GT1*H
            R = GT.MatMul(R);
            //cout << "R made\n" << endl;
            // Q = G1*G2*...*Gn -> QT = GTn*GTn-1*...*GT1
            // note: QT is a CSR Matrix because it's the product
            // of a bunch of CSR Matrices
            // cout << "fart\n"; 
            QT = GT.MatMul(QT);
            //cout << "QT made\n" << endl;
        }
        // if (m==2)
        // {
        //     cout << "R=\n";
        //     for (size_t i=0; i<m+1; i++)
        //     {
        //         for (size_t j=0; j<m; j++) cout << R[i][j] << ' ';
        //         cout << "\n";
        //     }
        //     cout << "Q=\n";
        //     for (size_t i=0; i<m+1; i++)
        //     {
        //         for (size_t j=0; j<m+1; j++) cout << QT.GetValue(j,i) << ' ';
        //         cout << "\n";
        //     }
        // }
        
        // TESTING ONLY: verify that QR = H
        // REMOVE LATER!!!!
        // cout << "testing H = QR\n" << endl;
        // for (size_t i=0; i<m+1; i++)
        // {
        //     for (size_t j=0; j<m; j++)
        //     {
        //         // cout << "(i,j)=(" << i << "," << j << ") ";
        //         double hij_calc = 0.0;
        //         for (size_t k=0; k<m+1; k++) hij_calc += QT.GetValue(k,i)*R[k][j];
        //         // cout << "hij_calc=" << hij_calc << " ";
        //         // cout << "H[i][j]=" << H[i][j] << "\n";
        //         assert(abs(hij_calc-H[i][j])<1e-3);
        //     }
        // }

        //now have R * ym = Q^-1 * beta*e1 = Q^T * beta*e1
        // solve Q * z = beta*e1 -> compute z = Q^T * beta * e1
        // where e1 is the first basis vector of the Krylov subspace
        // i.e. e1 = (1, 0, ..., 0)
        // idk why i'm calling it z it just needed a name
        vector<double> z (m+1, 0.0);
        for (size_t t=0; t<m+1; t++) z[t] = beta*QT.GetValue(t,0);

        // now solve R * ym = z with backwards substitution
        // ym[i] = (z[i] - sum(j=i+1 -> m) of R[i,j] * ym[j]) / R[i,i]
        // start at the highest i and go down to zero
        // because need to know all ym[j]'s beneath a certain
        // ym[i] to calculate it
        // USE INT HERE SO IT DOESN'T CRAP ITSELF
        // WHEN IT TRIES SUBTRACTING 1 FROM i=0
        for (int i=m-1; i>=0; i--)
        {
            // sum over all components "beneath" current ym[i]
            // cout << "i=" << i << "\n";
            double back_sum = 0.0;
            for (int j=i+1; j<m; j++) back_sum += R[i][j]*ym[j];
            // calculate ym[i]
            ym[i] = (z[i] - back_sum)/R[i][i];
        }
        
        // calculate residual while we're at it
        double res = 0.0;
        vector<double> r(m+1, 0.0);
        r[0] = beta; 
        // get residual vector r
        for (size_t i=0; i<m+1; i++)
        {
            // don't need to sum over whole row because 
            // we know H is 0 below the first subdiagonal
            // vvv this is because H isn't square
            size_t jstart = max(i-1, size_t(0));
            for (size_t j=jstart; j<m; j++) r[i] -= H[i][j]*ym[j];
            // add r[i]^2 to residual (squared)
            res += r[i]*r[i];
        }
        // get final residual
        res = sqrt(res);

        return {ym, res};
    }

    /**
     * @brief Does one iteration of GMRES to solve Ax = b for x
     * 
     * i.e. one set of iterations for a Krylov space of size `m`
     * 
     * @param A: (Matrix) the matrix A representing out linear system
     * @param phi0: (vector<double>) an initial guess at the potential vector we're solving for
     * @param b: (vector<double>) the RHS vector
     * @param m: the size of the Krylov space for the current set of iterations
     * 
     * @returns (tuple<vector<double>, double>) the final potential for this iteration and the 
     *          residual of the least squares solution `ym` used to compute the potential
     */
    tuple<vector<double>, double> MeshSerial::GMRES(const vector<double> phi0, const size_t m, const double tol)
    {
        // NOTE: ADD IN 4*PI*G???
        // calculate initial residual
        double mesh_size = this->_nx * this->_ny * this->_nz;
        vector<double> r0(mesh_size, 0.0);
        vector<double> Ax = this->_stencil.MatVec(phi0);
        for (size_t t=0; t<mesh_size; t++) r0[t] = this->_density[t] - Ax[t];
        // & 2-norm of residual (beta)
        double beta = 0.0;
        for (size_t t=0; t<mesh_size; t++) beta += r0[t]*r0[t];
        beta = sqrt(beta);
        // calculate first basis vector for Krylov subspace from first residual
        // store basis vectors in vector of vectors
        // NOTE: these vectors make up the COLUMNS of the Krylov subspace's basis matrix
        vector<vector<double>> V(m+1, vector<double>(mesh_size, 0.0));
        for (size_t t=0; t<V[0].size(); t++) V[0][t] = r0[t]/beta;

        // stores Hessenberg matrix elements
        // note: push_back saves on memory but is slower - hopefully this cuts down on time
        vector<vector<double>> H(m+1, vector<double>(m, 0.0));

        // Arnoldi iteration until m
        for (size_t j=0; j<m; j++)
        {
            // cout << "j=" << j << "\n" << endl;
            vector<double> wj = this->_stencil.MatVec(V[j]);
            
            // construct orthonormal basis with modified Gram-Schmidt
            for (size_t i=0; i<=j; i++)
            {
                // cout << "i=" << i << "\n" << endl;
                // dot product -> Hessenberg matrix element (i,j)
                for (size_t t=0; t<wj.size(); t++) H[i][j] += wj[t]*V[i][t];
                // update wj by subtracting hij*vi
                for (size_t t=0; t<wj.size(); t++) wj[t] -= H[i][j]*V[i][t];
            }
            
            // 2-norm of wj -> Hessenberg matrix element (j+1,j)
            for (size_t t=0; t<mesh_size; t++) H[j+1][j] += wj[t]*wj[t];
            // if the 2-norm of wj is within tol of zero, stop iterating here and get xj & residual
            if (H[j+1][j]<tol)
            {
                // Resize H to only have a size of (j+2)x(j+1)
                H.resize(j+2);
                for (vector<double> h: H) h.resize(j+1);
                // solve least squares for ym that minimizes the residual
                // get residual too while we're at it
                tuple<vector<double>, double> restuple = LeastSquares(H, beta);
                vector<double> ym = get<vector<double>>(restuple);
                double res = get<double>(restuple);

                // get the resulting approximate solution phim
                vector<double> phim = phi0;
                for (size_t t=0; t<phim.size(); t++)
                {
                    for (size_t n=0; n<j+1; n++) phim[t] += V[n][t]*ym[n];
                }
                // return phim and residual
                return {phim, res};

            }
            // otherwise set the next basis vector of the Krylov subspace and continue iterating
            for (size_t t=0; t<V[j+1].size(); t++) V[j+1][t] = wj[t]/H[j+1][j];
        }

        // // check what H looks like:
        // if (m<6)
        // {
        //     cout << "H=\n";
        //     for (size_t i=0; i<m+1; i++)
        //     {
        //         for (size_t j=0; j<m; j++) cout << H[i][j] << ' ';
        //         cout << '\n';
        //     }
        // }

        // if we didn't end early:
        // solve least squares for ym that minimizes the residual
        // get residual too while we're at it
        tuple<vector<double>, double> restuple = LeastSquares(H, beta);
        vector<double> ym = get<vector<double>>(restuple);
        double res = get<double>(restuple);

        // get the resulting approximate solution phim
        vector<double> phim = phi0;
        for (size_t t=0; t<phim.size(); t++)
        {
            for (size_t n=0; n<m+1; n++) phim[t] += V[n][t]*ym[n];
            if (t%1000==0) cout << phim[t] << "\n" << endl; 
        }
        // return phim & residual
        return {phim, res};
    }


    /**
     * @brief Calculate the gravitational potential across the mesh
     * 
     * Calculate the gravitational potential at each mesh cell using
     * (restarted) GMRES iteration.
     * 
     * @param restart: (size_t) how many iterations of GMRES to allow before 
     *                 resetting the Krylov subspace
     * @param maxiter: (size_t) maximum number of total iterations of GMRES
     *                 to complete before aborting
     * @param tol: (double) tolerance for determining convergence (residual 
     *              must be less than `tol`)
     */
    void MeshSerial::CalculatePotential(const size_t restart, const size_t maxiter, const double tol)
    {
        cout << "Calculating gravitational potential...\n" << endl;

        // temporary vector to store potential
        // TO DO:
        // calc from something OTHER than size of potential vector
        // in case it's not been calculated yet
        // NOTE: would it be faster to use the previous timestep's potential as a starting guess?
        // try!!!!
        vector<double> potential_temp(this->_nx*this->_ny*this->_nx, 0.0);
        // whether to stop loop
        bool stop = false; 
        // iterate at most until we hit maxiter
        // so keep track of global iterations
        size_t iter = 1;
        // this was probably not the most efficient implementation ever...
        while ((iter <= maxiter) && (!stop))
        {
            // do main GMRES loop until hit restart
            // then last set of iterations' result 
            // becomes initial guess for new one
            size_t m = 2;
            while ((m<=restart) && (iter <= maxiter) && (!stop))
            {
                cout << "iter=" << iter << "\n" << endl;
                // get approximation of potential & residual from GMRES iteration
                tuple<vector<double>, double> temptuple = GMRES(potential_temp, m, tol);
                // extract them from output tuple
                potential_temp = get<vector<double>>(temptuple);
                double res = get<double>(temptuple);

                cout << "phi[centerish] = " << potential_temp[potential_temp.size()/2] << "\n" << endl;
                cout << "residual = " << res << "\n" << endl;

                // check if residual is within tolerance - stop looping if so
                if (res < tol) stop = true;
                m++;
                iter++;
            }

        }
        // set potential to potential_temp no matter what for now?
        this->_potential = potential_temp;
        // what do if failure to converge?
        if (iter >= maxiter)
        {
            cout << "GMRES failed to converge within maximimum number of iterations (" << maxiter << ")!\n" << endl;
        }
        else
        {
            cout << "GMRES converged in " << iter-1 << " iterations!\n" << endl;
        }
    }

    /**
     * @brief Calculate the density across the mesh
     * 
     * Calculates the density at each mesh cell based on particles contained in the mesh.
     * 
     * Currently a placeholder for the numerical implementation.
     */
    void MeshSerial::CalculateDensity()
    {
        cout << "Calculating mesh density..." << endl;
    }

    /**
     * @brief Update the positions of particles contained in the mesh
     * 
     * Updates particle positions based on their velocities and the gravitational forces acting on them
     * 
     * Currently a placeholder for the numerical implementation
     */
    void MeshSerial::UpdateParticlePositions()
    {
        cout << "Updating particle positions..." << endl;
    }
    
    /** 
     * @brief Update the velocities of particles contained in the mesh
     * 
     * Updates particle velocities based on the gravitational acceleration calculated from the mesh
     * 
     * Currently a placeholder for the numerical implementation
     */
    void MeshSerial::UpdateParticleVelocities()
    {
        cout << "Updating particle velocities..." << endl;
    }

    /**
     * @brief Initialize particles within the mesh
     * 
     * Initializes the particles and associates them with their corresponding mesh cells.
     * 
     * Currently a placeholder for the numerical implementation
     */
    void MeshSerial::InitializeParticles()
    {
        cout << "Initializing particles in serial mesh..." << endl;
    }

    size_t MeshSerial::Calc1DIndex(const size_t i, const size_t j, const size_t k)
    {
        return (this->_nx)*(this->_nx)*i + (this->_ny)*j + k;
    }
    
}