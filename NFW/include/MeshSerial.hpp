#ifndef __NBody_MeshSerial__
#define __NBody_MeshSerial__

#include <stddef.h>

#include "Mesh.hpp"
#include "Matrix.hpp"

using namespace std;

namespace NBody
{
    class MeshSerial : public Mesh
    {
        private:
            // mesh attributes
            vector<double> _potential;
            vector<double> _density; // NOTE: fold 4piG into density inherently?
            // size of mesh along x, y, and z directions
            // in number of mesh cells
            size_t _nx;
            size_t _ny;
            size_t _nz;
            // size of mesh cells along each dimension
            double _dx;
            double _dy;
            double _dz;
            // value of G set to 1 for now, will update with actual units later once everything is working.
            double _G = 1;
            // for GMRES - 7 point stencil matrix 
            // with 0 boundary conditions
            // NOTE TO SELF: ADD X,Y,Z SIZE 
            // OF GRID TO MESH CLASS
            Matrix _stencil{1,1}; // 
            // function to create stencil
            Matrix _ConstructStencil();
            // NOTE: WILL HAVE TO REWRITE A BUNCH OF STUFF
            // TO USE ARRAYS INSTEAD OF VECTORS OF VECTORS
            // FUCK

        public:
            // Constructors
            MeshSerial() = default;
            // ADD CONSTRUCTOR WHERE GIVE IT SIM DIMENSIONS?
            // MODIFY LATER
            MeshSerial(size_t nx, size_t ny, size_t nz, double dx, double dy, double dz);

            // Copy and move constructors
            MeshSerial(const MeshSerial& other) = default;
            MeshSerial& operator=(const MeshSerial& other) = default;
            MeshSerial(MeshSerial&& other) noexcept = default;
            MeshSerial& operator=(MeshSerial&& other) noexcept = default;

            // Destructor
            ~MeshSerial() override = default;

            // Getters
            [[nodiscard]] const vector<double>& getPotential() const noexcept { return _potential; }
            [[nodiscard]] const vector<double>& getDensity() const noexcept { return _density; }

            // Methods

            void CalculatePotential(const size_t restart, const size_t maxiter, const double tol); //override;
            void CalculateDensity() override;
            void UpdateParticlePositions() override;
            void UpdateParticleVelocities() override;
            void InitializeParticles() override;
            size_t Calc1DIndex(const size_t i, const size_t j, const size_t k);
            // needed for GMRES 
            tuple<vector<double>, double> GMRES(const vector<double> x0, const size_t m, const double tol); // return only result vector for now
            tuple<vector<double>, double> LeastSquares(vector<vector<double>> H, const double beta);
            NBody::Matrix ConstructGivens(vector<double> v, size_t i, size_t j);
    };
}

#endif