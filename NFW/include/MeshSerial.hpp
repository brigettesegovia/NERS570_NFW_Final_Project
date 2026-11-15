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
            Matrix _potential;
            Matrix _density;

        public:
            // Constructors
            MeshSerial() = default;

            // Copy and move constructors
            MeshSerial(const MeshSerial& other) = default;
            MeshSerial& operator=(const MeshSerial& other) = default;
            MeshSerial(MeshSerial&& other) noexcept = default;
            MeshSerial& operator=(MeshSerial&& other) noexcept = default;

            // Destructor
            ~MeshSerial() override = default;

            // Getters
            [[nodiscard]] const Matrix& getPotential() const noexcept { return _potential; }
            [[nodiscard]] const Matrix& getDensity() const noexcept { return _density; }

            // Methods
            void CalculatePotential() override;
            void CalculateDensity() override;
            void UpdateParticlePositions() override;
            void UpdateParticleVelocities() override;
            void InitializeParticles() override;
    };
}

#endif