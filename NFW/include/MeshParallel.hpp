#ifndef __NBody_MeshParallel__
#define __NBody_MeshParallel__

#include <stddef.h>

#include "Mesh.hpp"

using namespace std;

namespace NBody
{
    class MeshParallel : public Mesh
    {
        private:
            // mesh attributes
            double*** _potential;
            double*** _density;

        public:
            // Constructors
            MeshParallel() = default;

            // Copy and move constructors
            MeshParallel(const MeshParallel& other) = default;
            MeshParallel& operator=(const MeshParallel& other) = default;
            MeshParallel(MeshParallel&& other) noexcept = default;
            MeshParallel& operator=(MeshParallel&& other) noexcept = default;

            // Destructor
            ~MeshParallel() override = default;

            // Getters
            [[nodiscard]] double*** getPotential() const noexcept { return _potential; }
            [[nodiscard]] double*** getDensity() const noexcept { return _density; }

            // Methods
            void CalculatePotential() override;
            void CalculateDensity() override;
            void UpdateParticlePositions() override;
            void UpdateParticleVelocities() override;
            void InitializeParticles() override;
    };
}

#endif