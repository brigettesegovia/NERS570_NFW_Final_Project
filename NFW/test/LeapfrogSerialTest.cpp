/**
 * @brief Runs the serial leapfrog integrator in NBodySerial.hpp
 *
 * 1. Circular orbit in a fixed NFW halo   (orbit radius and energy should stay constant)
 * 2. Eccentric orbit in a fixed NFW halo  (2nd-order convergence + time reversibility)
 * 3. PM force on a uniform sphere         (compare with the analytic -G M(<r) / r^2)
 * 4. Self-gravitating PM run of a halo, printing energy/virial diagnostics as it goes
 *
 * Build (from the NFW/ directory):
 *     g++ -std=c++17 -O2 -Iinclude test/LeapfrogSerialTest.cpp src/DMparticle.cpp -o leapfrog_serial
 *
 * Run:
 *     ./leapfrog_serial [N] [ngrid] [tEnd] [dt] [snapshotPrefix] [densityTol]
 *       N              number of particles             (default 20000)
 *       ngrid          grid points per side            (default 64)
 *       tEnd           end time in kpc/(km/s) ~ Gyr    (default 1.5)
 *       dt             timestep                        (default 0.002)
 *       snapshotPrefix write <prefix>_00000.csv, ...   (default: no snapshots)
 *       densityTol     stop early once the radial density profile changes by less
 *                      than this between outputs       (default 0 = run to tEnd)
 */

#include "NBodySerial.hpp"

#include <random>
#include <chrono>
#include <string>

using namespace std;
using namespace NBody;

// halo used for tests 1 & 2: roughly Milky-Way sized NFW (M200 ~ 1e12 Msun, c ~ 10)
constexpr double RHO_S = 6.7e6;  // Msun / kpc^3
constexpr double R_S   = 20.0;   // kpc


/**
 * @brief Uniform-density sphere with isotropic Gaussian velocities
 *
 * Placeholder initial conditions so the PM run has something to evolve;
 * swap in the real initialization when it's ready. The velocity dispersion is
 * set by the virial ratio Q = 2K/|W| (Q < 1 collapses, Q = 1 is in equilibrium).
 *
 * @param N: number of particles
 * @param M: total mass [Msun]
 * @param R: sphere radius [kpc]
 * @param Q: virial ratio
 * @param seed: random seed
 */
vector<DMparticle> UniformSphere(const size_t N, const double M, const double R, const double Q, const unsigned seed)
{
    mt19937_64 rng(seed);
    uniform_real_distribution<double> uniform(-1.0, 1.0);

    // uniform sphere: W = -3/5 G M^2 / R, K = 3/2 M sigma^2  ->  sigma^2 = Q G M / (5 R)
    const double sigma = sqrt(Q * G_GRAV * M / (5.0 * R));
    normal_distribution<double> gauss(0.0, sigma);

    vector<DMparticle> particles;
    particles.reserve(N);
    const double m = M / static_cast<double>(N);
    vector<double> comPos = {0, 0, 0}, comVel = {0, 0, 0};

    while (particles.size() < N)
    {
        // rejection-sample a point in the unit ball
        const double x = uniform(rng), y = uniform(rng), z = uniform(rng);
        if (x*x + y*y + z*z > 1.0) { continue; }
        const vector<double> pos = {R*x, R*y, R*z};
        const vector<double> vel = {gauss(rng), gauss(rng), gauss(rng)};
        for (size_t d = 0; d < 3; d++) { comPos[d] += pos[d] / N; comVel[d] += vel[d] / N; }
        particles.emplace_back(m, pos, vel);
    }

    // shift to the center-of-mass frame so the halo stays centered in the box
    for (DMparticle& p : particles)
    {
        const vector<double> x = p.getPos(), v = p.getVel();
        p.setPos(x[0] - comPos[0], x[1] - comPos[1], x[2] - comPos[2]);
        p.setVel(v[0] - comVel[0], v[1] - comVel[1], v[2] - comVel[2]);
    }
    return particles;
}


/**
 * @brief Test 1: a particle on a circular orbit at r = r_s, with V_c from Eq. (2)
 */
bool TestCircularOrbit()
{
    AnalyticNFWField nfw(RHO_S, R_S);
    const double r0 = R_S;
    const double vc = nfw.CircularVelocity(r0);
    const double period = 2.0 * PI * r0 / vc;
    const size_t stepsPerOrbit = 500;
    const size_t nOrbits = 10;

    vector<DMparticle> particles = { DMparticle(1.0, {r0, 0.0, 0.0}, {0.0, vc, 0.0}) };
    LeapfrogSerial leapfrog(particles, nfw, period / stepsPerOrbit);
    const double E0 = leapfrog.TotalEnergy();

    double maxRadiusErr = 0.0, maxEnergyErr = 0.0;
    for (size_t s = 0; s < nOrbits * stepsPerOrbit; s++)
    {
        leapfrog.Step();
        const vector<double> x = leapfrog.getParticles()[0].getPos();
        const double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        maxRadiusErr = max(maxRadiusErr, fabs(r - r0) / r0);
        maxEnergyErr = max(maxEnergyErr, fabs((leapfrog.TotalEnergy() - E0) / E0));
    }

    const bool pass = (maxRadiusErr < 1e-3) && (maxEnergyErr < 1e-6);
    cout << "Test 1: circular orbit at r_s = " << R_S << " kpc, V_c = " << vc << " km/s, period = "
         << period << " kpc/(km/s)\n"
         << "        " << nOrbits << " orbits: max |dr|/r0 = " << maxRadiusErr
         << ", max |dE/E0| = " << maxEnergyErr << (pass ? "   PASS" : "   FAIL") << "\n\n";
    return pass;
}


/**
 * @brief Max |dE/E0| over a few orbits of an eccentric orbit with timestep dt
 */
double EccentricEnergyError(const AnalyticNFWField& nfwIn, const double dt, const double tEnd)
{
    AnalyticNFWField nfw = nfwIn;
    const double r0 = 50.0;
    const double vc = nfw.CircularVelocity(r0);
    vector<DMparticle> particles = { DMparticle(1.0, {r0, 0.0, 0.0}, {0.2*vc, 0.5*vc, 0.1*vc}) };
    LeapfrogSerial leapfrog(particles, nfw, dt);
    const double E0 = leapfrog.TotalEnergy();
    double maxErr = 0.0;
    const size_t nSteps = static_cast<size_t>(round(tEnd / dt));
    for (size_t s = 0; s < nSteps; s++)
    {
        leapfrog.Step();
        maxErr = max(maxErr, fabs((leapfrog.TotalEnergy() - E0) / E0));
    }
    return maxErr;
}


/**
 * @brief Test 2: eccentric orbit - halving dt should cut the energy error ~4x (2nd order),
 *        and running forward then backward should return to the start (time reversible)
 */
bool TestEccentricOrbit()
{
    AnalyticNFWField nfw(RHO_S, R_S);
    const double tEnd = 5.0;  // ~ several orbits at r ~ 50 kpc

    const double err1 = EccentricEnergyError(nfw, 0.01,  tEnd);
    const double err2 = EccentricEnergyError(nfw, 0.005, tEnd);
    const double err3 = EccentricEnergyError(nfw, 0.0025, tEnd);
    const double order12 = log2(err1 / err2);
    const double order23 = log2(err2 / err3);
    const bool orderPass = (order12 > 1.8 && order12 < 2.2) && (order23 > 1.8 && order23 < 2.2);

    // time reversibility: forward N steps, flip velocities, forward N steps again
    const double r0 = 50.0, vc = nfw.CircularVelocity(r0);
    const vector<double> x0 = {r0, 0.0, 0.0};
    const vector<double> v0 = {0.2*vc, 0.5*vc, 0.1*vc};
    LeapfrogSerial forward({DMparticle(1.0, x0, v0)}, nfw, 0.01);
    for (size_t s = 0; s < 500; s++) { forward.Step(); }

    DMparticle flipped = forward.getParticles()[0];
    const vector<double> v = flipped.getVel();
    flipped.setVel(-v[0], -v[1], -v[2]);
    LeapfrogSerial backward({flipped}, nfw, 0.01);
    for (size_t s = 0; s < 500; s++) { backward.Step(); }

    const vector<double> xBack = backward.getParticles()[0].getPos();
    const double returnErr = sqrt(pow(xBack[0]-x0[0], 2) + pow(xBack[1]-x0[1], 2) + pow(xBack[2]-x0[2], 2)) / r0;
    const bool reversePass = returnErr < 1e-9;

    cout << "Test 2: eccentric orbit, max |dE/E0| for dt = 0.01, 0.005, 0.0025: "
         << err1 << ", " << err2 << ", " << err3 << "\n"
         << "        convergence order = " << order12 << ", " << order23 << " (expect 2)"
         << (orderPass ? "   PASS" : "   FAIL") << "\n"
         << "        forward 500 steps then backward 500 steps: |x - x0|/r0 = " << returnErr
         << (reversePass ? "   PASS" : "   FAIL") << "\n\n";
    return orderPass && reversePass;
}


/**
 * @brief Test 3: PM acceleration of a uniform sphere vs. the analytic -G M r / R^3
 */
bool TestPMForce()
{
    const size_t N = 100000;
    const double M = 1e12, R = 100.0, L = 400.0;
    const size_t ngrid = 64;
    vector<DMparticle> particles = UniformSphere(N, M, R, 0.0, 42);

    PMGravitySerial pm(ngrid, L);
    vector<array<double, 3>> acc;
    pm.CalculateAccelerations(particles, acc);

    // compare radial acceleration with the analytic one away from the center
    // (the PM force is softened on the scale of ~2 grid cells)
    double sumRelErr = 0.0;
    size_t count = 0;
    for (size_t p = 0; p < N; p++)
    {
        const vector<double> x = particles[p].getPos();
        const double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        if (r < 3.0 * pm.getSpacing()) { continue; }
        const double aRadial = (acc[p][0]*x[0] + acc[p][1]*x[1] + acc[p][2]*x[2]) / r;
        const double aExact  = -G_GRAV * M * r / (R*R*R);
        sumRelErr += fabs((aRadial - aExact) / aExact);
        count++;
    }
    const double meanRelErr = sumRelErr / count;
    const bool pass = meanRelErr < 0.05;

    cout << "Test 3: PM force on a uniform sphere (N = " << N << ", " << ngrid << "^3 grid, h = "
         << pm.getSpacing() << " kpc)\n"
         << "        " << pm.Status() << "\n"
         << "        mean |a_PM - a_exact| / |a_exact| for r > 3h: " << meanRelErr
         << (pass ? "   PASS" : "   FAIL") << "\n\n";
    return pass;
}


int main(int argc, char* argv[])
{
    // ------------------------------------------------------------------
    // command-line options
    // ------------------------------------------------------------------
    const size_t N          = (argc > 1) ? stoul(argv[1]) : 20000;
    const size_t ngrid      = (argc > 2) ? stoul(argv[2]) : 64;
    const double tEnd       = (argc > 3) ? stod(argv[3])  : 1.5;
    const double dt         = (argc > 4) ? stod(argv[4])  : 0.002;
    const string prefix     = (argc > 5) ? argv[5]        : "";
    const double densityTol = (argc > 6) ? stod(argv[6])  : 0.0;

    // ------------------------------------------------------------------
    // checks of the integrator and the PM force
    // ------------------------------------------------------------------
    cout << "===== Leapfrog checks =====\n\n";
    bool allPass = true;
    allPass &= TestCircularOrbit();
    allPass &= TestEccentricOrbit();
    allPass &= TestPMForce();

    // ------------------------------------------------------------------
    // self-gravitating PM run
    // ------------------------------------------------------------------
    const double M = 1e12;     // halo mass [Msun]
    const double R = 100.0;    // initial radius [kpc]
    const double Q = 0.5;      // initial virial ratio (< 1: halo collapses, then relaxes)
    const double L = 4.0 * R;  // box side [kpc]

    cout << "===== Serial PM leapfrog run =====\n"
         << "N = " << N << ", grid = " << ngrid << "^3, box = " << L << " kpc, dt = " << dt
         << ", tEnd = " << tEnd << " kpc/(km/s)\n"
         << "initial halo: uniform sphere, M = " << M << " Msun, R = " << R << " kpc, 2K/|W| = " << Q << "\n\n";

    vector<DMparticle> particles = UniformSphere(N, M, R, Q, 12345);
    PMGravitySerial pm(ngrid, L);
    LeapfrogSerial leapfrog(particles, pm, dt);
    leapfrog.setProfileBins(2.0 * pm.getSpacing(), R, 12);

    const auto start = chrono::steady_clock::now();
    leapfrog.Run(tEnd, 50, prefix, densityTol);
    const double seconds = chrono::duration<double>(chrono::steady_clock::now() - start).count();

    cout << "\n" << leapfrog.getStep() << " steps in " << seconds << " s ("
         << 1000.0 * seconds / max<size_t>(leapfrog.getStep(), 1) << " ms/step)\n";
    cout << (allPass ? "All checks passed.\n" : "Some checks FAILED.\n");
    return allPass ? 0 : 1;
}
