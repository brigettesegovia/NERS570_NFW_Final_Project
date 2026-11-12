#ifndef __NBody_DMparticle__
#define __NBody_DMparticle__

#include <stddef.h>
#include <vector>

using namespace std;

namespace NBody
{
    class DMparticle
    {
        private:
            // everything zeroed out by default!
            double _mass = 0;                // particle mass
            vector<double> _pos = {0, 0, 0}; // 3D particle position
            vector<double> _vel = {0, 0, 0}; // 3D particle velocity
            // other attributes?
        
        public:
            // constructors
            DMparticle() = default; // default
            DMparticle(const double mass); // give it mass only
            DMparticle(const double mass, const vector<double> pos, const vector<double> vel); // give it mass, position, & velocity
            // copy and move constructors
            DMparticle(const DMparticle& other) = default;                 // copy
            DMparticle& operator=(const DMparticle& other) = default;      // copy assignment
            DMparticle(DMparticle&& other) noexcept = default;             // move
            DMparticle& operator=(DMparticle&& other) noexcept = default;  // move assginment

            // destructor
            ~DMparticle() = default;

            // getters - defined fully in header because they're simple
            [[nodiscard]] double         getMass() const noexcept { return _mass; }
            [[nodiscard]] vector<double> getPos()  const noexcept { return _pos; }
            [[nodiscard]] vector<double> getVel()  const noexcept { return _vel; }

            // setters
            void setMass(const double mass);
            // can set position/velocity either with a vector or three separate values
            void setPos(const vector<double> pos);
            void setPos(const double pos_x, const double pos_y, const double pos_z);
            void setVel(const vector<double> vel);
            void setVel(const double vel_x, const double vel_y, const double vel_z);

            // other methods?
    };
}

#endif