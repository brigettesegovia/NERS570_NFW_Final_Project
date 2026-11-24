#include "MeshSerial.hpp"

#include <stddef.h>
#include <vector>
#include <cassert>
#include <iostream>
#include <cmath>
#include <tuple>

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
     * @brief computes the least squares solution to argmin||beta*e1 - Hj*yj||2
     * 
     */
    tuple<vector<double>, double> MeshSerial::LeastSquares(vector<vector<double>> H, const double beta, const size_t j)
    {
        vector<double> yj(j+1, 0.0);
        double res = 0.0;
        // do least squares
        return {yj, res};
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
     */
    tuple<vector<double>, double> MeshSerial::GMRES(const vector<double> phi0, const size_t m, const double tol)
    {
        // calculate initial residual
        vector<double> r0(this->_density.size(), 0.0);
        vector<double> Ax = this->_stencil.MatVec(phi0);
        for (size_t t=0; t<r0.size(); t++) r0[t] = this->_density[t] - Ax[t];
        // & 2-norm of residual (beta)
        double beta = 0.0;
        for (size_t t=0; t<r0.size(); t++) beta += r0[t]*r0[t];
        beta = sqrt(beta);
        // calculate first basis vector for Krylov subspace from first residual
        // store basis vectors in vector of vectors
        // NOTE: these vectors make up the COLUMNS of the Krylov subspace's basis matrix
        vector<vector<double>> V(m+1, vector<double>(r0.size(), 0.0));
        for (size_t t=0; t<V[0].size(); t++) V[0][t] = r0[t]/beta;

        // stores Hessenberg matrix elements
        // note: push_back saves on memory but is slower - hopefully this cuts down on time
        vector<vector<double>> H(m+1, vector<double>(m, 0.0));

        // Arnoldi iteration until m
        for (size_t j=0; j<m; j++)
        {
            vector<double> wj = this->_stencil.MatVec(V[j]);
            
            // construct orthonormal basis with modified Gram-Schmidt
            for (size_t i=0; i<j; i++)
            {
                // dot product -> Hessenberg matrix element (i,j)
                for (size_t t=0; t<wj.size(); t++) H[i][j] += wj[t]*V[i][t];
                // update wj by subtracting hij*vi
                for (size_t t=0; t<wj.size(); t++) wj[t] -= H[i][j]*V[i][t];
            }
            
            // 2-norm of wj -> Hessenberg matrix element (j+1,j)
            for (size_t t=0; t<r0.size(); t++) H[j+1][j] += wj[t]*wj[t];
            // if the 2-norm of wj is within tol of zero, stop iterating here and get xj & residual
            if (H[j+1][j]<tol)
            {
                // solve least squares for ym that minimizes the residual
                // get residual too while we're at it
                tuple<vector<double>, double> restuple = LeastSquares(H, beta, j);
                vector<double> ym = get<vector<double>>(restuple);
                double res = get<double>(restuple);

                // get the resulting approximate solution phim
                vector<double> phim(phi0.size(), 0.0);
                for (size_t t=0; t<phim.size(); t++)
                {
                    for (size_t n=0; n<j+1; n++) phim[t] += V[n][t]*ym[n];
                }
                // return phim - will calculate residual outside function
                return {phim, res};

            }
            // otherwise set the next basis vector of the Krylov subspace and continue iterating
            for (size_t t=0; t<V[j+1].size(); t++) V[j+1][t] = wj[t]/H[j+1][j];
        }

        // if we didn't end early:
        // solve least squares for ym that minimizes the residual
        // get residual too while we're at it
        tuple<vector<double>, double> restuple = LeastSquares(H, beta, m);
        vector<double> ym = get<vector<double>>(restuple);
        double res = get<double>(restuple);

        // get the resulting approximate solution phim
        vector<double> phim(phi0.size(), 0.0);
        for (size_t t=0; t<phim.size(); t++)
        {
            for (size_t n=0; n<m+1; n++) phim[t] += V[n][t]*ym[n];
        }
        // return phim - will calculate residual outside function
        return {phim, res};
    }


    /**
     * @brief Calculate the gravitational potential across the mesh
     * 
     * Calculate the gravitational potential at each mesh cell using
     * (restarted) GMRES iteration.
     * 
     * @param restart: (size_t) how many iterations of GMRES to allow before 
     * resetting the Krylov subspace
     * @param maxiter: (size_t) maximum number of total iterations of GMRES
     * to complete before aborting
     * @param tol: (double) tolerance for determining convergence
     */
    void MeshSerial::CalculatePotential(const size_t restart, const size_t maxiter, const double tol)
    {
        cout << "Calculating gravitational potential..." << endl;

        // temporary vector to store potential
        // TO DO:
        // calc from something OTHER than size of potential vector
        // in case it's not been calculated yet
        vector<double> potential_temp(this->_potential.size(), 0.0);

        // iterate at most until we hit maxiter
        // so keep track of global iterations
        size_t iter = 0;
        while (iter < maxiter)
        {
            // do main GMRES loop until hit restart
            // then last set of iterations' result 
            // becomes initial guess for new one
            size_t m = 0;
            while ((m<=restart) && (iter < maxiter))
            {
                // get approximation of potential & residual from GMRES iteration
                tuple<vector<double>, double> temptuple = GMRES(potential_temp, m, tol);
                // extract them from output tuple
                potential_temp = get<vector<double>>(temptuple);
                double res = get<double>(temptuple);

                // check if residual is within tolerance
                // if (res < tol) 
                // TO DO: figure this loop out

            }

        }
        // what do if failure to converge?
        cout << "Failed to converge within maphimimum number of iterations (" << maxiter << ")!\n" << endl;
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

    size_t Calc1DIndex(const size_t i, const size_t j, const size_t k)
    {
        cout << "Calculating 1D index..." << endl;
        return 0;
    }
    
}