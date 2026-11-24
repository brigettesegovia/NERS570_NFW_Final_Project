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
            vector<double> _density;
            // for GMRES - 7 point stencil matrix 
            // with 0 boundary conditions
            // NOTE TO SELF: ADD X,Y,Z SIZE 
            // OF GRID TO MESH CLASS
            const NBody::Matrix _stencil;

        public:
            // Constructors
            MeshSerial() = default;
            // ADD CONSTRUCTOR WHERE GIVE IT SIM DIMENSIONS

            // Copy and move constructors
            MeshSerial(const MeshSerial& other) = default;
            MeshSerial& operator=(const MeshSerial& other) = default;
            MeshSerial(MeshSerial&& other) noexcept = default;
            MeshSerial& operator=(MeshSerial&& other) noexcept = default;

            // Destructor
            ~MeshSerial() override = default;

            // Getters
            [[nodiscard]] const Matrix& getPotential() const noexcept; //{ return _potential; }
            [[nodiscard]] const Matrix& getDensity() const noexcept; //{ return _density; }

            // Methods

            void CalculatePotential(const size_t restart, const size_t maxiter, const double tol); //override;
            void CalculateDensity() override;
            void UpdateParticlePositions() override;
            void UpdateParticleVelocities() override;
            void InitializeParticles() override;
            size_t Calc1DIndex(const size_t i, const size_t j, const size_t k);
            tuple<vector<double>, double> GMRES(const vector<double> x0, const size_t m, const double tol); // return only result vector for now
            tuple<vector<double>, double> LeastSquares(vector<vector<double>> H, const double beta, const size_t j);
    };
}

#endif