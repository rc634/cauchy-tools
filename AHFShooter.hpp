#pragma once
#include "Surface.hpp"

class Spacetime;

class AHFShooter : public Surface {
public:

    std::vector<double> g;   // df/dsigma, tracked during RK4 integration

    double rmax  = 8.;       // upper radius clamp
    double rmin  = 1e-5;     // lower radius clamp
    bool   OOB   = false;    // out-of-bounds flag set during shoot()
    double g_END = 0.;       // target value of g at the equator (boundary condition)

    AHFShooter(int npoints);

    void initialize(double f0);

    // ODE right-hand side: dg/dsigma given (f, g, sigma)
    double DGDS(const Spacetime& spacetime, double F, double G, double TH);

    // single RK4 integration pass from pole to equator with initial radius rH0
    void shoot(const Spacetime& spacetime, double rH0);

    // bisection search over initial radius to satisfy the equatorial BC g_END = 0
    void interval_bisection(const Spacetime& spacetime);

    void hello() const;
};
