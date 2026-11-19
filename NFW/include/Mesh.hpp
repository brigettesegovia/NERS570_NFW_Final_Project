#ifndef __NBody_Mesh__
#define __NBody_Mesh__

#include <stddef.h>
#include <vector>
#include <array>

#include "DMparticle.hpp"

using namespace std;

namespace NBody
{
    class Mesh
    { 
        protected: 
            // mesh attributes
            // position
            size_t*** _position;

            // particle association
            size_t*** _particleAssociation;

            // particles contained in the mesh
            vector<DMparticle> _particles;

        public:
            // constructor
            Mesh() = default;

            // copy and move constructors
            Mesh(const Mesh& other) = default;
            Mesh& operator=(const Mesh& other) = default;
            Mesh(Mesh&& other) noexcept = default;
            Mesh& operator=(Mesh&& other) noexcept = default;

            // destructor
            virtual ~Mesh() = default;

            // getters
            [[nodiscard]] size_t*** getPosition() const noexcept { return _position; }
            
            // abstract methods
            virtual void CalculatePotential() = 0;
            virtual void CalculateDensity() = 0;
            virtual void UpdateParticlePositions() = 0;
            virtual void UpdateParticleVelocities() = 0;
            virtual void InitializeParticles() = 0;

    };
}

#endif