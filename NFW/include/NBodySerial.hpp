#ifndef __NBodySerial__
#define __NBodySerial__

#include "MeshSerial.hpp"
#include "Matrix.hpp"
#include "DMparticle.hpp"

#include <stddef.h>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>

using namespace std;

// =========================================================================
// Serial leapfrog integration for the particle-mesh (PM) N-body code
// =========================================================================
//
// Every timestep follows the flowchart in the proposal:
//
//   assign particles to the mesh -> density -> potential (del^2 phi = 4 pi G rho)
//   -> interpolate force back to each particle (a = -grad phi) -> leapfrog update
//
// The leapfrog integrator (LeapfrogSerial) only needs the acceleration of each
// particle, so it talks to an abstract AccelerationField. Two are provided here:
//
//   AnalyticNFWField : fixed analytic NFW potential. No self-gravity; used to
//                      test the integrator (circular orbits, energy conservation).
//   PMGravitySerial  : simple serial particle mesh self-gravity (cloud-in-cell
//                      deposit + red-black SOR Poisson solve + finite-difference
//                      gradient). A stand-in until MeshSerial is filled in; once
//                      it is, MeshSerial can inherit from AccelerationField and
//                      plug straight into LeapfrogSerial.
//
// Nothing here assumes the particles start uniformly distributed in 6D phase
// space: the integrator evolves whatever positions/velocities it is handed.
//
// Units used throughout:
//   length   : kpc
//   velocity : km/s
//   mass     : Msun
//   time     : kpc / (km/s) ~= 0.978 Gyr
// =========================================================================

namespace NBody
{
    constexpr double G_GRAV = 4.30091e-6;             // gravitational constant [kpc (km/s)^2 / Msun]
    constexpr double PI     = 3.14159265358979323846;

    // =====================================================================
    // Class declarations
    // =====================================================================

    /**
     * @brief Abstract interface for anything that can give the gravitational
     *        acceleration felt by each particle
     */
    class AccelerationField
    {
        public:
            virtual ~AccelerationField() = default;

            // fill acc[p] = {ax, ay, az} for every particle p at its current position
            virtual void CalculateAccelerations(const vector<DMparticle>& particles,
                                                vector<array<double, 3>>& acc) = 0;

            // total gravitational potential energy of the particles
            virtual double PotentialEnergy(const vector<DMparticle>& particles) const = 0;

            // optional one-line status for the run log (e.g. solver iterations)
            [[nodiscard]] virtual string Status() const { return ""; }
    };


    /**
     * @brief Fixed (external) analytic NFW potential
     *
     * Particles feel the NFW halo but not each other, so a test particle's
     * orbit and energy are known exactly. Used to check the integrator.
     */
    class AnalyticNFWField : public AccelerationField
    {
        private:
            double _rho_s;  // characteristic density [Msun / kpc^3]
            double _r_s;    // scale radius           [kpc]
            double _G;      // gravitational constant

        public:
            // constructor
            AnalyticNFWField(const double rho_s, const double r_s, const double G = G_GRAV);

            // getters
            [[nodiscard]] double getRhoS() const noexcept { return _rho_s; }
            [[nodiscard]] double getRs()   const noexcept { return _r_s; }

            // NFW profile quantities as a function of radius
            [[nodiscard]] double Density(const double r) const;
            [[nodiscard]] double EnclosedMass(const double r) const;
            [[nodiscard]] double Potential(const double r) const;
            [[nodiscard]] double CircularVelocity(const double r) const;

            // AccelerationField interface
            void CalculateAccelerations(const vector<DMparticle>& particles,
                                        vector<array<double, 3>>& acc) override;
            double PotentialEnergy(const vector<DMparticle>& particles) const override;
    };


    /**
     * @brief Simple serial particle-mesh self-gravity
     *
     * Cubic box of side L centered on the origin with n grid points per side
     * (spacing h = L/(n-1)). Each call to CalculateAccelerations does:
     *
     *   1. cloud-in-cell (CIC) mass assignment       -> _density
     *   2. red-black SOR solve of del^2 phi = 4 pi G rho -> _potential
     *      (isolated halo: boundary potential = point mass at the grid's
     *       center of mass, interior warm-started from the previous step)
     *   3. a = -grad phi by central differences, CIC-interpolated to particles
     *
     * Particles in the outermost layer of cells (or outside the box) do not
     * deposit mass; they feel the monopole (point-mass) force of the grid mass.
     */
    class PMGravitySerial : public AccelerationField
    {
        private:
            size_t _n;        // grid points along each axis (n^3 total)
            double _L;        // side length of the box [kpc]
            double _h;        // grid spacing = L/(n-1) [kpc]
            double _xmin;     // coordinate of the first grid point = -L/2
            double _G;        // gravitational constant
            double _tol;      // relative residual tolerance for the Poisson solve
            size_t _maxIter;  // max SOR iterations per solve
            double _omega;    // SOR over-relaxation factor

            vector<double> _density;    // mass density at each grid point [Msun/kpc^3]
            vector<double> _potential;  // potential at each grid point    [(km/s)^2]

            // bookkeeping from the most recent call
            double            _massOnGrid = 0;             // mass deposited on the grid
            array<double, 3>  _comOnGrid  = {0, 0, 0};     // center of mass of that mass
            size_t            _nOutside   = 0;             // particles that fell off the grid
            size_t            _lastIterations = 0;         // SOR iterations used
            double            _lastResidual   = 0;         // relative residual reached
            bool              _potentialInitialized = false;

            // flatten (i, j, k) -> index into _density / _potential
            [[nodiscard]] size_t idx(const size_t i, const size_t j, const size_t k) const noexcept
            { return (i*_n + j)*_n + k; }

            // helpers
            bool   Locate(const array<double, 3>& x, array<size_t, 3>& cell, array<double, 3>& frac) const;
            static double CICWeight(const array<double, 3>& frac, const size_t di, const size_t dj, const size_t dk);
            double InterpolatePotential(const array<size_t, 3>& cell, const array<double, 3>& frac) const;
            double MonopolePotential(const array<double, 3>& x) const;
            array<double, 3> MonopoleAcceleration(const array<double, 3>& x) const;
            double ResidualNorm() const;

        public:
            // constructor
            PMGravitySerial(const size_t n, const double L, const double G = G_GRAV,
                            const double tol = 1e-5, const size_t maxIter = 20000);

            // getters
            [[nodiscard]] size_t getN()              const noexcept { return _n; }
            [[nodiscard]] double getL()              const noexcept { return _L; }
            [[nodiscard]] double getSpacing()        const noexcept { return _h; }
            [[nodiscard]] double getMassOnGrid()     const noexcept { return _massOnGrid; }
            [[nodiscard]] size_t getNumOutside()     const noexcept { return _nOutside; }
            [[nodiscard]] size_t getLastIterations() const noexcept { return _lastIterations; }
            [[nodiscard]] double getLastResidual()   const noexcept { return _lastResidual; }
            [[nodiscard]] const vector<double>& getDensity()   const noexcept { return _density; }
            [[nodiscard]] const vector<double>& getPotential() const noexcept { return _potential; }

            // PM steps
            void CalculateDensity(const vector<DMparticle>& particles);
            void CalculatePotential();

            // AccelerationField interface
            void CalculateAccelerations(const vector<DMparticle>& particles,
                                        vector<array<double, 3>>& acc) override;
            double PotentialEnergy(const vector<DMparticle>& particles) const override;
            [[nodiscard]] string Status() const override;
    };


    /**
     * @brief Serial kick-drift-kick (KDK) leapfrog integrator
     *
     * One step of size dt:
     *     v_{i+1/2} = v_i       + a_i     dt/2     (half kick)
     *     x_{i+1}   = x_i       + v_{i+1/2} dt     (drift)
     *     a_{i+1}   = a(x_{i+1})                   (density -> potential -> force)
     *     v_{i+1}   = v_{i+1/2} + a_{i+1} dt/2     (half kick)
     *
     * Merging the last half kick of one step with the first half kick of the
     * next gives exactly the proposal's v_{i+1/2} = v_{i-1/2} + a_i dt,
     * x_{i+1} = x_i + v_{i+1/2} dt. Writing it as KDK keeps x and v at the same
     * time, so energies/snapshots are synchronized. Only one force evaluation
     * per step (a_{i+1} is reused as a_i of the next step).
     */
    class LeapfrogSerial
    {
        private:
            vector<DMparticle>       _particles;          // particles being evolved
            AccelerationField&       _field;              // where accelerations come from
            double                   _dt;                 // timestep [kpc/(km/s)]
            double                   _time = 0;           // current simulation time
            size_t                   _step = 0;           // number of completed steps
            vector<array<double, 3>> _acc;                // accelerations at the current positions
            bool                     _accCurrent = false; // are _acc up to date with the positions?
            double                   _E0 = 0;             // energy at the start of Run()
            bool                     _E0set = false;
            vector<double>           _profileEdges;       // radial bin edges for the density-convergence check

        public:
            // constructor: takes its own copy of the particles
            LeapfrogSerial(vector<DMparticle> particles, AccelerationField& field, const double dt);

            // copy and move constructors (assignment deleted: the reference to the field can't be reseated)
            LeapfrogSerial(const LeapfrogSerial& other) = default;
            LeapfrogSerial& operator=(const LeapfrogSerial& other) = delete;
            LeapfrogSerial(LeapfrogSerial&& other) noexcept = default;
            LeapfrogSerial& operator=(LeapfrogSerial&& other) noexcept = delete;

            // destructor
            ~LeapfrogSerial() = default;

            // getters
            [[nodiscard]] const vector<DMparticle>&       getParticles()     const noexcept { return _particles; }
            [[nodiscard]] const vector<array<double, 3>>& getAccelerations() const noexcept { return _acc; }
            [[nodiscard]] double getTime() const noexcept { return _time; }
            [[nodiscard]] double getDt()   const noexcept { return _dt; }
            [[nodiscard]] size_t getStep() const noexcept { return _step; }

            // setters
            void setDt(const double dt);
            void setProfileBins(const double rmin, const double rmax, const size_t nbins);

            // leapfrog pieces
            void UpdateAccelerations();
            void Kick(const double dt);
            void Drift(const double dt);
            void Step();
            void Run(const double tEnd, const size_t outputEvery = 100, const string& snapshotPrefix = "test_result_runs",
                     const double densityTol = 0.0, ostream& logStream = cout);

            // diagnostics
            [[nodiscard]] double KineticEnergy() const;
            double PotentialEnergy();
            double TotalEnergy();
            [[nodiscard]] array<double, 3> AngularMomentum() const;
            [[nodiscard]] array<double, 3> CenterOfMass() const;
            [[nodiscard]] vector<double> RadialDensityProfile(const vector<double>& edges) const;
            void PrintDiagnosticsHeader(ostream& out) const;
            void PrintDiagnostics(ostream& out);
            void WriteSnapshot(const string& filename) const;
    };


    // =====================================================================
    // AnalyticNFWField
    // =====================================================================

    /**
     * @brief Constructs an analytic NFW potential
     *
     * @param rho_s: characteristic density [Msun/kpc^3] (double)
     * @param r_s:   scale radius [kpc] (double)
     * @param G:     gravitational constant (double), default kpc (km/s)^2 / Msun
     */
    inline AnalyticNFWField::AnalyticNFWField(const double rho_s, const double r_s, const double G) :
        _rho_s(rho_s), _r_s(r_s), _G(G)
    {
        assert(this->_rho_s > 0);
        assert(this->_r_s > 0);
        assert(this->_G > 0);
    }

    /**
     * @brief NFW density rho(r) = rho_s / [(r/r_s)(1 + r/r_s)^2]
     *
     * @param r: radius, must be > 0 (double)
     * @return density [Msun/kpc^3] (double)
     */
    inline double AnalyticNFWField::Density(const double r) const
    {
        assert(r > 0);
        const double x = r / this->_r_s;
        return this->_rho_s / (x * (1.0 + x) * (1.0 + x));
    }

    /**
     * @brief Mass enclosed within radius r
     *
     * M(<r) = 4 pi rho_s r_s^3 [ln(1 + x) - x/(1 + x)],  x = r/r_s
     *
     * @param r: radius (double)
     * @return enclosed mass [Msun] (double)
     */
    inline double AnalyticNFWField::EnclosedMass(const double r) const
    {
        assert(r >= 0);
        const double x = r / this->_r_s;
        double mu;
        // ln(1+x) - x/(1+x) subtracts two nearly equal numbers for tiny x,
        // so use its Taylor series there: x^2/2 - 2x^3/3 + 3x^4/4
        if (x < 1e-3) { mu = x*x*(0.5 - x*(2.0/3.0 - 0.75*x)); }
        else          { mu = log1p(x) - x/(1.0 + x); }
        return 4.0*PI*this->_rho_s*this->_r_s*this->_r_s*this->_r_s*mu;
    }

    /**
     * @brief NFW potential phi(r) = -4 pi G rho_s r_s^3 ln(1 + r/r_s) / r
     *
     * @param r: radius (double)
     * @return potential [(km/s)^2] (double)
     */
    inline double AnalyticNFWField::Potential(const double r) const
    {
        assert(r >= 0);
        const double x = r / this->_r_s;
        // ln(1+x)/x -> 1 as x -> 0, so the potential is finite at the center
        const double lnterm = (x < 1e-8) ? 1.0 - 0.5*x : log1p(x)/x;
        return -4.0*PI*this->_G*this->_rho_s*this->_r_s*this->_r_s*lnterm;
    }

    /**
     * @brief Circular velocity V_c(r) = sqrt(G M(<r) / r)  (Eq. 2 of the proposal)
     *
     * @param r: radius (double)
     * @return circular velocity [km/s] (double)
     */
    inline double AnalyticNFWField::CircularVelocity(const double r) const
    {
        if (r <= 0) { return 0.0; }
        return sqrt(this->_G * this->EnclosedMass(r) / r);
    }

    /**
     * @brief Acceleration of each particle in the NFW potential, a = -G M(<r) r_vec / r^3
     *
     * @param particles: particles (vector<DMparticle>)
     * @param acc:       output accelerations, resized to particles.size() (vector<array<double,3>>)
     */
    inline void AnalyticNFWField::CalculateAccelerations(const vector<DMparticle>& particles,
                                                         vector<array<double, 3>>& acc)
    {
        acc.assign(particles.size(), {0.0, 0.0, 0.0});
        for (size_t p = 0; p < particles.size(); p++)
        {
            const vector<double> pos = particles[p].getPos();
            const double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
            if (r == 0.0) { continue; } // force vanishes at the exact center by symmetry
            const double factor = -this->_G * this->EnclosedMass(r) / (r*r*r);
            for (size_t d = 0; d < 3; d++) { acc[p][d] = factor * pos[d]; }
        }
    }

    /**
     * @brief Potential energy of the particles in the external NFW field, W = sum m phi(r)
     *
     * @param particles: particles (vector<DMparticle>)
     * @return potential energy [Msun (km/s)^2] (double)
     */
    inline double AnalyticNFWField::PotentialEnergy(const vector<DMparticle>& particles) const
    {
        double W = 0.0;
        for (const DMparticle& p : particles)
        {
            const vector<double> pos = p.getPos();
            const double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
            W += p.getMass() * this->Potential(r);
        }
        return W;
    }


    // =====================================================================
    // PMGravitySerial
    // =====================================================================

    /**
     * @brief Constructs a serial particle-mesh gravity solver
     *
     * @param n:       grid points along each axis, >= 5 (size_t)
     * @param L:       side length of the cubic box centered on the origin [kpc] (double)
     * @param G:       gravitational constant (double)
     * @param tol:     relative residual tolerance for the Poisson solve (double)
     * @param maxIter: max SOR iterations per solve (size_t)
     */
    inline PMGravitySerial::PMGravitySerial(const size_t n, const double L, const double G,
                                            const double tol, const size_t maxIter) :
        _n(n), _L(L), _h(L / static_cast<double>(n - 1)), _xmin(-0.5*L),
        _G(G), _tol(tol), _maxIter(maxIter)
    {
        assert(n >= 5);  // need a few interior points
        assert(L > 0);
        assert(tol > 0);
        // optimal over-relaxation factor for the 7-point Laplacian with Dirichlet boundaries
        this->_omega = 2.0 / (1.0 + sin(PI / static_cast<double>(n - 1)));
        this->_density.assign(n*n*n, 0.0);
        this->_potential.assign(n*n*n, 0.0);
    }

    /**
     * @brief Finds the grid cell a position falls in
     *
     * The CIC cloud touches nodes i and i+1 along each axis, and both must be
     * interior nodes (1..n-2) because the boundary nodes hold the fixed potential.
     *
     * @param x:    position (array<double,3>)
     * @param cell: output lower-corner node indices (array<size_t,3>)
     * @param frac: output fractional offset from that node, in [0,1) (array<double,3>)
     * @return true if the particle is in the usable part of the grid (bool)
     */
    inline bool PMGravitySerial::Locate(const array<double, 3>& x, array<size_t, 3>& cell,
                                        array<double, 3>& frac) const
    {
        for (size_t d = 0; d < 3; d++)
        {
            // position in units of the grid spacing, measured from the first grid point
            const double s = (x[d] - this->_xmin) / this->_h;
            // (written this way so NaN positions also count as outside)
            if (!(s >= 1.0 && s < static_cast<double>(this->_n - 2))) { return false; }
            cell[d] = static_cast<size_t>(s);  // floor, since s > 0
            frac[d] = s - static_cast<double>(cell[d]);
        }
        return true;
    }

    /**
     * @brief CIC weight of corner (di, dj, dk) of a cell
     *
     * @param frac: fractional offset within the cell (array<double,3>)
     * @param di, dj, dk: 0 or 1, which corner (size_t)
     * @return weight (double); the 8 weights sum to 1
     */
    inline double PMGravitySerial::CICWeight(const array<double, 3>& frac,
                                             const size_t di, const size_t dj, const size_t dk)
    {
        return (di ? frac[0] : 1.0 - frac[0])
             * (dj ? frac[1] : 1.0 - frac[1])
             * (dk ? frac[2] : 1.0 - frac[2]);
    }

    /**
     * @brief CIC-interpolates the mesh potential to a point inside the grid
     */
    inline double PMGravitySerial::InterpolatePotential(const array<size_t, 3>& cell,
                                                        const array<double, 3>& frac) const
    {
        double phi = 0.0;
        for (size_t di = 0; di < 2; di++)
        for (size_t dj = 0; dj < 2; dj++)
        for (size_t dk = 0; dk < 2; dk++)
        {
            phi += CICWeight(frac, di, dj, dk) * this->_potential[idx(cell[0]+di, cell[1]+dj, cell[2]+dk)];
        }
        return phi;
    }

    /**
     * @brief Potential of the grid mass treated as a point at its center of mass
     *
     * Softened at one grid spacing so it never blows up.
     */
    inline double PMGravitySerial::MonopolePotential(const array<double, 3>& x) const
    {
        if (this->_massOnGrid <= 0.0) { return 0.0; }
        const double dx = x[0] - this->_comOnGrid[0];
        const double dy = x[1] - this->_comOnGrid[1];
        const double dz = x[2] - this->_comOnGrid[2];
        const double r  = max(sqrt(dx*dx + dy*dy + dz*dz), this->_h);
        return -this->_G * this->_massOnGrid / r;
    }

    /**
     * @brief Acceleration from the grid mass treated as a point at its center of mass
     *
     * Used for particles that have left the usable part of the grid.
     */
    inline array<double, 3> PMGravitySerial::MonopoleAcceleration(const array<double, 3>& x) const
    {
        array<double, 3> a = {0.0, 0.0, 0.0};
        if (this->_massOnGrid <= 0.0) { return a; }
        array<double, 3> d;
        for (size_t c = 0; c < 3; c++) { d[c] = x[c] - this->_comOnGrid[c]; }
        const double r = max(sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]), this->_h);
        const double factor = -this->_G * this->_massOnGrid / (r*r*r);
        for (size_t c = 0; c < 3; c++) { a[c] = factor * d[c]; }
        return a;
    }

    /**
     * @brief L2 norm of the Poisson residual, del^2 phi - 4 pi G rho, over the interior nodes
     */
    inline double PMGravitySerial::ResidualNorm() const
    {
        const size_t n = this->_n;
        const size_t sx = n*n, sy = n;   // index strides along x and y (z stride is 1)
        const double invh2   = 1.0 / (this->_h * this->_h);
        const double fourPiG = 4.0 * PI * this->_G;
        const vector<double>& P = this->_potential;

        double sum = 0.0;
        for (size_t i = 1; i < n-1; i++)
        for (size_t j = 1; j < n-1; j++)
        for (size_t k = 1; k < n-1; k++)
        {
            const size_t c = idx(i, j, k);
            const double lap = (P[c+sx] + P[c-sx] + P[c+sy] + P[c-sy] + P[c+1] + P[c-1] - 6.0*P[c]) * invh2;
            const double r = lap - fourPiG * this->_density[c];
            sum += r*r;
        }
        return sqrt(sum);
    }

    /**
     * @brief Cloud-in-cell assignment of particle mass to the grid
     *
     * @param particles: particles (vector<DMparticle>)
     */
    inline void PMGravitySerial::CalculateDensity(const vector<DMparticle>& particles)
    {
        fill(this->_density.begin(), this->_density.end(), 0.0);
        this->_massOnGrid = 0.0;
        this->_comOnGrid  = {0.0, 0.0, 0.0};
        this->_nOutside   = 0;

        const double invVol = 1.0 / (this->_h * this->_h * this->_h);
        array<size_t, 3> cell;
        array<double, 3> frac;

        for (const DMparticle& p : particles)
        {
            const vector<double> pv = p.getPos();
            const array<double, 3> x = {pv[0], pv[1], pv[2]};
            const double m = p.getMass();

            if (!this->Locate(x, cell, frac)) { this->_nOutside++; continue; }

            // spread the mass over the 8 surrounding grid points
            for (size_t di = 0; di < 2; di++)
            for (size_t dj = 0; dj < 2; dj++)
            for (size_t dk = 0; dk < 2; dk++)
            {
                this->_density[idx(cell[0]+di, cell[1]+dj, cell[2]+dk)] += m * CICWeight(frac, di, dj, dk) * invVol;
            }

            this->_massOnGrid += m;
            for (size_t d = 0; d < 3; d++) { this->_comOnGrid[d] += m * x[d]; }
        }

        if (this->_massOnGrid > 0.0)
        {
            for (size_t d = 0; d < 3; d++) { this->_comOnGrid[d] /= this->_massOnGrid; }
        }
    }

    /**
     * @brief Solves del^2 phi = 4 pi G rho on the grid with red-black SOR
     *
     * Boundary nodes are fixed to the monopole potential of the grid mass.
     * The interior starts from the previous step's solution (warm start), so
     * after the first step only a handful of iterations are usually needed.
     */
    inline void PMGravitySerial::CalculatePotential()
    {
        const size_t n = this->_n;
        const size_t sx = n*n, sy = n;
        vector<double>& P = this->_potential;

        // nothing on the grid -> potential is zero everywhere
        if (this->_massOnGrid <= 0.0)
        {
            fill(P.begin(), P.end(), 0.0);
            this->_lastIterations = 0;
            this->_lastResidual = 0.0;
            return;
        }

        // boundary conditions (and initial guess on the very first call)
        for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
        for (size_t k = 0; k < n; k++)
        {
            const bool boundary = (i == 0 || j == 0 || k == 0 || i == n-1 || j == n-1 || k == n-1);
            if (boundary || !this->_potentialInitialized)
            {
                const array<double, 3> x = {this->_xmin + i*this->_h, this->_xmin + j*this->_h, this->_xmin + k*this->_h};
                P[idx(i, j, k)] = this->MonopolePotential(x);
            }
        }
        this->_potentialInitialized = true;

        // size of the right-hand side, so the stopping test is relative
        const double fourPiG = 4.0 * PI * this->_G;
        double rhsNorm = 0.0;
        for (size_t i = 1; i < n-1; i++)
        for (size_t j = 1; j < n-1; j++)
        for (size_t k = 1; k < n-1; k++)
        {
            const double f = fourPiG * this->_density[idx(i, j, k)];
            rhsNorm += f*f;
        }
        rhsNorm = sqrt(rhsNorm);

        const double h2 = this->_h * this->_h;
        size_t iter = 0;
        double relRes = this->ResidualNorm() / rhsNorm;

        while (relRes > this->_tol && iter < this->_maxIter)
        {
            // red-black ordering: update all nodes with (i+j+k) even, then all with (i+j+k) odd
            for (size_t color = 0; color < 2; color++)
            {
                for (size_t i = 1; i < n-1; i++)
                for (size_t j = 1; j < n-1; j++)
                {
                    const size_t k0 = ((i + j + 1) % 2 == color) ? 1 : 2;
                    for (size_t k = k0; k < n-1; k += 2)
                    {
                        const size_t c = idx(i, j, k);
                        const double neighbors = P[c+sx] + P[c-sx] + P[c+sy] + P[c-sy] + P[c+1] + P[c-1];
                        const double gaussSeidel = (neighbors - h2 * fourPiG * this->_density[c]) / 6.0;
                        P[c] += this->_omega * (gaussSeidel - P[c]);
                    }
                }
            }
            iter++;
            // checking the residual costs about one sweep, so don't do it every iteration
            if (iter % 5 == 0 || iter == this->_maxIter) { relRes = this->ResidualNorm() / rhsNorm; }
        }

        this->_lastIterations = iter;
        this->_lastResidual = relRes;
        if (relRes > this->_tol)
        {
            cerr << "PMGravitySerial: Poisson solve hit maxIter = " << this->_maxIter
                 << " with relative residual " << relRes << endl;
        }
    }

    /**
     * @brief Particle-mesh acceleration of every particle
     *
     * density (CIC) -> potential (SOR) -> -grad phi at the nodes -> CIC interpolation
     *
     * @param particles: particles (vector<DMparticle>)
     * @param acc:       output accelerations (vector<array<double,3>>)
     */
    inline void PMGravitySerial::CalculateAccelerations(const vector<DMparticle>& particles,
                                                        vector<array<double, 3>>& acc)
    {
        this->CalculateDensity(particles);
        this->CalculatePotential();

        acc.assign(particles.size(), {0.0, 0.0, 0.0});
        const size_t sx = this->_n * this->_n, sy = this->_n;
        const double inv2h = 1.0 / (2.0 * this->_h);
        const vector<double>& P = this->_potential;
        array<size_t, 3> cell;
        array<double, 3> frac;

        for (size_t p = 0; p < particles.size(); p++)
        {
            const vector<double> pv = particles[p].getPos();
            const array<double, 3> x = {pv[0], pv[1], pv[2]};

            if (!this->Locate(x, cell, frac))
            {
                acc[p] = this->MonopoleAcceleration(x);
                continue;
            }

            // all 8 nodes are interior, so their neighbors (used by the central difference) exist
            for (size_t di = 0; di < 2; di++)
            for (size_t dj = 0; dj < 2; dj++)
            for (size_t dk = 0; dk < 2; dk++)
            {
                const double w = CICWeight(frac, di, dj, dk);
                const size_t q = idx(cell[0]+di, cell[1]+dj, cell[2]+dk);
                acc[p][0] -= w * (P[q+sx] - P[q-sx]) * inv2h;
                acc[p][1] -= w * (P[q+sy] - P[q-sy]) * inv2h;
                acc[p][2] -= w * (P[q+1]  - P[q-1])  * inv2h;
            }
        }
    }

    /**
     * @brief Self-gravitational potential energy from the current mesh potential
     *
     * W = 1/2 sum_(on grid) m phi(x)  +  sum_(off grid) m phi_monopole(x)
     *
     * Must be called after CalculateAccelerations with the same particle positions.
     */
    inline double PMGravitySerial::PotentialEnergy(const vector<DMparticle>& particles) const
    {
        double W = 0.0;
        array<size_t, 3> cell;
        array<double, 3> frac;
        for (const DMparticle& p : particles)
        {
            const vector<double> pv = p.getPos();
            const array<double, 3> x = {pv[0], pv[1], pv[2]};
            if (this->Locate(x, cell, frac)) { W += 0.5 * p.getMass() * this->InterpolatePotential(cell, frac); }
            else                             { W += p.getMass() * this->MonopolePotential(x); }
        }
        return W;
    }

    /**
     * @brief One-line summary of the last Poisson solve for the run log
     */
    inline string PMGravitySerial::Status() const
    {
        ostringstream s;
        s << "SOR iters=" << this->_lastIterations
          << " res=" << scientific << setprecision(1) << this->_lastResidual
          << " off-grid=" << this->_nOutside;
        return s.str();
    }


    // =====================================================================
    // LeapfrogSerial
    // =====================================================================

    /**
     * @brief Constructs the leapfrog integrator
     *
     * @param particles: initial particles, any distribution in 6D phase space (vector<DMparticle>)
     * @param field:     source of accelerations, must outlive the integrator (AccelerationField&)
     * @param dt:        timestep [kpc/(km/s)] (double)
     */
    inline LeapfrogSerial::LeapfrogSerial(vector<DMparticle> particles, AccelerationField& field, const double dt) :
        _particles(std::move(particles)), _field(field), _dt(dt)
    {
        assert(this->_dt > 0);
        assert(!this->_particles.empty());
        this->_acc.assign(this->_particles.size(), {0.0, 0.0, 0.0});
    }

    /**
     * @brief Changes the timestep
     *
     * @param dt: new timestep (double)
     */
    inline void LeapfrogSerial::setDt(const double dt)
    {
        assert(dt > 0);
        this->_dt = dt;
    }

    /**
     * @brief Sets log-spaced radial bins (around the center of mass) used by the
     *        density-convergence stopping test in Run()
     *
     * @param rmin:  inner edge of the first bin (double)
     * @param rmax:  outer edge of the last bin (double)
     * @param nbins: number of bins (size_t)
     */
    inline void LeapfrogSerial::setProfileBins(const double rmin, const double rmax, const size_t nbins)
    {
        assert(rmin > 0 && rmax > rmin && nbins > 0);
        this->_profileEdges.resize(nbins + 1);
        for (size_t b = 0; b <= nbins; b++)
        {
            this->_profileEdges[b] = rmin * pow(rmax / rmin, static_cast<double>(b) / static_cast<double>(nbins));
        }
    }

    /**
     * @brief Recomputes the acceleration of every particle at its current position
     */
    inline void LeapfrogSerial::UpdateAccelerations()
    {
        this->_field.CalculateAccelerations(this->_particles, this->_acc);
        this->_accCurrent = true;
    }

    /**
     * @brief Velocity update ("kick"): v <- v + a dt
     *
     * @param dt: time interval (double)
     */
    inline void LeapfrogSerial::Kick(const double dt)
    {
        assert(this->_accCurrent);
        for (size_t p = 0; p < this->_particles.size(); p++)
        {
            const vector<double> v = this->_particles[p].getVel();
            this->_particles[p].setVel(v[0] + this->_acc[p][0]*dt,
                                       v[1] + this->_acc[p][1]*dt,
                                       v[2] + this->_acc[p][2]*dt);
        }
    }

    /**
     * @brief Position update ("drift"): x <- x + v dt
     *
     * @param dt: time interval (double)
     */
    inline void LeapfrogSerial::Drift(const double dt)
    {
        for (DMparticle& p : this->_particles)
        {
            const vector<double> x = p.getPos();
            const vector<double> v = p.getVel();
            p.setPos(x[0] + v[0]*dt, x[1] + v[1]*dt, x[2] + v[2]*dt);
        }
        // particles moved, so the stored accelerations are stale
        this->_accCurrent = false;
    }

    /**
     * @brief Advances the system by one kick-drift-kick step of size dt
     */
    inline void LeapfrogSerial::Step()
    {
        // a_i at the current positions (only actually computed on the first step;
        // afterwards it is left over from the end of the previous step)
        if (!this->_accCurrent) { this->UpdateAccelerations(); }

        this->Kick(0.5 * this->_dt);   // v_{i+1/2} = v_i + a_i dt/2
        this->Drift(this->_dt);        // x_{i+1}   = x_i + v_{i+1/2} dt
        this->UpdateAccelerations();   // a_{i+1}   (density -> potential -> force)
        this->Kick(0.5 * this->_dt);   // v_{i+1}   = v_{i+1/2} + a_{i+1} dt/2

        this->_time += this->_dt;
        this->_step++;
    }

    /**
     * @brief Integrates until tEnd, or until the density profile stops changing
     *
     * @param tEnd:           time to stop at (double)
     * @param outputEvery:    print diagnostics (and snapshot / convergence check) every this many steps, 0 = never (size_t)
     * @param snapshotPrefix: if not empty, write <prefix>_00000.csv, ... at each output (string)
     * @param densityTol:     if > 0, stop once the RMS fractional change of the radial density
     *                        profile between outputs is below this; needs setProfileBins() (double)
     * @param logStream:      where diagnostics are printed (ostream)
     */
    inline void LeapfrogSerial::Run(const double tEnd, const size_t outputEvery, const string& snapshotPrefix,
                                    const double densityTol, ostream& logStream)
    {
        if (!this->_E0set) { this->_E0 = this->TotalEnergy(); this->_E0set = true; }

        vector<double> oldProfile;
        if (densityTol > 0)
        {
            assert(!this->_profileEdges.empty() && "call setProfileBins() before using densityTol");
            oldProfile = this->RadialDensityProfile(this->_profileEdges);
        }

        size_t snapshotCount = 0;
        auto snapshotName = [&](size_t count)
        {
            ostringstream name;
            name << snapshotPrefix << "_" << setw(5) << setfill('0') << count << ".csv";
            return name.str();
        };

        this->PrintDiagnosticsHeader(logStream);
        this->PrintDiagnostics(logStream);
        if (!snapshotPrefix.empty()) { this->WriteSnapshot(snapshotName(snapshotCount++)); }

        // count steps up front so round-off in _time can't add an extra step
        const double remaining = tEnd - this->_time;
        const size_t nSteps = (remaining > 0) ? static_cast<size_t>(ceil(remaining / this->_dt - 1e-9)) : 0;

        bool printedLast = true;
        for (size_t s = 0; s < nSteps; s++)
        {
            this->Step();
            printedLast = false;

            if (outputEvery > 0 && this->_step % outputEvery == 0)
            {
                this->PrintDiagnostics(logStream);
                printedLast = true;
                if (!snapshotPrefix.empty()) { this->WriteSnapshot(snapshotName(snapshotCount++)); }

                if (densityTol > 0)
                {
                    const vector<double> newProfile = this->RadialDensityProfile(this->_profileEdges);
                    double sum = 0.0;
                    size_t used = 0;
                    for (size_t b = 0; b < newProfile.size(); b++)
                    {
                        if (oldProfile[b] > 0)
                        {
                            const double frac = (newProfile[b] - oldProfile[b]) / oldProfile[b];
                            sum += frac*frac;
                            used++;
                        }
                    }
                    const double change = (used > 0) ? sqrt(sum / used) : 0.0;
                    logStream << "    density profile change = " << change << endl;
                    if (change < densityTol)
                    {
                        logStream << "    density profile converged (< " << densityTol << "), stopping at t = "
                            << this->_time << endl;
                        break;
                    }
                    oldProfile = newProfile;
                }
            }
        }

        if (!printedLast) { this->PrintDiagnostics(logStream); }
    }

    /**
     * @brief Total kinetic energy, K = sum 1/2 m v^2
     */
    inline double LeapfrogSerial::KineticEnergy() const
    {
        double K = 0.0;
        for (const DMparticle& p : this->_particles)
        {
            const vector<double> v = p.getVel();
            K += 0.5 * p.getMass() * (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        }
        return K;
    }

    /**
     * @brief Total potential energy (recomputes accelerations first if positions changed)
     */
    inline double LeapfrogSerial::PotentialEnergy()
    {
        if (!this->_accCurrent) { this->UpdateAccelerations(); }
        return this->_field.PotentialEnergy(this->_particles);
    }

    /**
     * @brief Total energy, K + W
     */
    inline double LeapfrogSerial::TotalEnergy()
    {
        return this->KineticEnergy() + this->PotentialEnergy();
    }

    /**
     * @brief Total angular momentum about the origin, L = sum m (x cross v)
     */
    inline array<double, 3> LeapfrogSerial::AngularMomentum() const
    {
        array<double, 3> L = {0.0, 0.0, 0.0};
        for (const DMparticle& p : this->_particles)
        {
            const vector<double> x = p.getPos();
            const vector<double> v = p.getVel();
            const double m = p.getMass();
            L[0] += m * (x[1]*v[2] - x[2]*v[1]);
            L[1] += m * (x[2]*v[0] - x[0]*v[2]);
            L[2] += m * (x[0]*v[1] - x[1]*v[0]);
        }
        return L;
    }

    /**
     * @brief Center of mass of all particles
     */
    inline array<double, 3> LeapfrogSerial::CenterOfMass() const
    {
        array<double, 3> com = {0.0, 0.0, 0.0};
        double M = 0.0;
        for (const DMparticle& p : this->_particles)
        {
            const vector<double> x = p.getPos();
            for (size_t d = 0; d < 3; d++) { com[d] += p.getMass() * x[d]; }
            M += p.getMass();
        }
        for (size_t d = 0; d < 3; d++) { com[d] /= M; }
        return com;
    }

    /**
     * @brief Spherically averaged density in radial shells around the center of mass
     *
     * Simple version for the stopping test; the full profile/NFW fit lives elsewhere.
     *
     * @param edges: shell edges, increasing (vector<double>)
     * @return mass / shell volume for each of the edges.size()-1 shells (vector<double>)
     */
    inline vector<double> LeapfrogSerial::RadialDensityProfile(const vector<double>& edges) const
    {
        assert(edges.size() >= 2);
        const size_t nbins = edges.size() - 1;
        vector<double> rho(nbins, 0.0);
        const array<double, 3> com = this->CenterOfMass();

        for (const DMparticle& p : this->_particles)
        {
            const vector<double> x = p.getPos();
            const double r = sqrt((x[0]-com[0])*(x[0]-com[0]) + (x[1]-com[1])*(x[1]-com[1]) + (x[2]-com[2])*(x[2]-com[2]));
            if (r < edges.front() || r >= edges.back()) { continue; }
            // first edge greater than r; the bin is the one before it
            const size_t b = static_cast<size_t>(upper_bound(edges.begin(), edges.end(), r) - edges.begin()) - 1;
            rho[b] += p.getMass();
        }
        for (size_t b = 0; b < nbins; b++)
        {
            const double vol = 4.0/3.0 * PI * (pow(edges[b+1], 3) - pow(edges[b], 3));
            rho[b] /= vol;
        }
        return rho;
    }

    /**
     * @brief Prints the column names for PrintDiagnostics
     */
    inline void LeapfrogSerial::PrintDiagnosticsHeader(ostream& out) const
    {
        out << setw(8)  << "step"
            << setw(12) << "time"
            << setw(14) << "K"
            << setw(14) << "W"
            << setw(14) << "E"
            << setw(12) << "dE/|E0|"
            << setw(12) << "|L|"
            << setw(9)  << "2K/|W|"
            << "   field" << endl;
    }

    /**
     * @brief Prints one line of conserved quantities for monitoring the run
     */
    inline void LeapfrogSerial::PrintDiagnostics(ostream& out)
    {
        const double K = this->KineticEnergy();
        const double W = this->PotentialEnergy();
        const double E = K + W;
        if (!this->_E0set) { this->_E0 = E; this->_E0set = true; }
        const double dE = (this->_E0 != 0.0) ? (E - this->_E0) / fabs(this->_E0) : 0.0;
        const array<double, 3> L = this->AngularMomentum();
        const double Lmag = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]);
        const double virial = (W != 0.0) ? 2.0*K/fabs(W) : 0.0;

        const ios_base::fmtflags oldFlags = out.flags();
        const streamsize oldPrecision = out.precision();
        out << setw(8)  << this->_step
            << fixed << setprecision(4) << setw(12) << this->_time
            << scientific << setprecision(5)
            << setw(14) << K
            << setw(14) << W
            << setw(14) << E
            << setprecision(2)
            << setw(12) << dE
            << setw(12) << Lmag
            << fixed << setprecision(3) << setw(9) << virial
            << "   " << this->_field.Status() << endl;
        out.flags(oldFlags);
        out.precision(oldPrecision);
    }

    /**
     * @brief Writes all particles to a CSV file: x,y,z,vx,vy,vz,mass
     *
     * @param filename: output file (string)
     */
    inline void LeapfrogSerial::WriteSnapshot(const string& filename) const
    {
        ofstream file(filename);
        if (!file) { cerr << "LeapfrogSerial: could not open " << filename << endl; return; }
        file << "x,y,z,vx,vy,vz,mass\n";
        file << setprecision(10);
        for (const DMparticle& p : this->_particles)
        {
            const vector<double> x = p.getPos();
            const vector<double> v = p.getVel();
            file << x[0] << "," << x[1] << "," << x[2] << ","
                 << v[0] << "," << v[1] << "," << v[2] << "," << p.getMass() << "\n";
        }
    }
}

#endif
