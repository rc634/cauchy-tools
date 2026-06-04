#include "AHFShooter.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>
#include <cmath>
#include <algorithm>

AHFShooter::AHFShooter(int npoints) : Surface(npoints) {
    g.resize(npoints, 0.0);
}

void AHFShooter::initialize(double f0) {
    Surface::initialize(f0);
    for (int i = 0; i < num_points; ++i) g[i] = 0.;
}

double AHFShooter::DGDS(const Spacetime& spacetime, double F, double G, double TH) {
    if (F <= 0.) F = rmin;
    double x = F * sin(TH);
    double y = F * cos(TH);

    double ddx   = spacetime.get_ddx_interp(spacetime.psi, x, y);
    double ddy   = spacetime.get_ddy_interp(spacetime.psi, x, y);
    double PSI   = spacetime.get_val_interp(spacetime.psi, x, y);
    double DRPSI = ddx * sin(TH) + ddy * cos(TH);
    double PSIP  = ddx * cos(TH) - ddy * sin(TH);

    return - (G*G*G/F/F + G) * (4.*PSIP/PSI + cos(TH)/sin(TH))
           + G*G * (3./F + 4.*DRPSI/PSI)
           + 2.*F + 4.*F*F*DRPSI/PSI;
}

void AHFShooter::shoot(const Spacetime& spacetime, double rH0) {
    double F0, F1, F2, F3;
    double G0, G1, G2, G3;
    double f0, f1, f2, f3;
    double g0, g1, g2, g3;
    double sig0, sig1, sig2, sig3;
    OOB = false;

    // pole boundary condition: small-angle approximation for g[0]
    double x0    = rH0 * sin(sigma[0]);
    double y0    = rH0 * cos(sigma[0]);
    double psi0  = spacetime.get_val_interp(spacetime.psi, x0, y0);
    double dR0   = spacetime.get_ddx_interp(spacetime.psi, x0, y0) * sin(sigma[0])
                 + spacetime.get_ddy_interp(spacetime.psi, x0, y0) * cos(sigma[0]);
    f[0] = rH0;
    g[0] = (rH0 + 2.*rH0*rH0*dR0/psi0) * sigma[0];

    // RK4 integration from pole to equator
    for (int i = 0; i < num_points - 1; ++i) {
        f0 = f[i];  g0 = g[i];  sig0 = sigma[i];

        F0 = g0;               G0 = DGDS(spacetime, f0, g0, sig0);
        f1 = f0 + 0.5*ds*F0;  g1 = g0 + 0.5*ds*G0;  sig1 = sig0 + 0.5*ds;
        F1 = g1;               G1 = DGDS(spacetime, f1, g1, sig1);
        f2 = f0 + 0.5*ds*F1;  g2 = g0 + 0.5*ds*G1;  sig2 = sig0 + 0.5*ds;
        F2 = g2;               G2 = DGDS(spacetime, f2, g2, sig2);
        f3 = f0 + ds*F2;       g3 = g0 + ds*G2;       sig3 = sig0 + ds;
        F3 = g3;               G3 = DGDS(spacetime, f3, g3, sig3);

        f[i+1] = f[i] + (F0 + 2.*F1 + 2.*F2 + F3) * ds / 6.;
        g[i+1] = g[i] + (G0 + 2.*G1 + 2.*G2 + G3) * ds / 6.;

        f[i+1] = std::min(f[i+1], rmax);
        f[i+1] = std::max(f[i+1], rmin);
        if (f[i+1] >= rmax) OOB = true;
    }

    // equatorial boundary condition: small-angle approximation for g at equator
    int iE      = num_points - 1;
    double xE   = f[iE] * sin(sigma[iE]);
    double yE   = f[iE] * cos(sigma[iE]);
    double psiE = spacetime.get_val_interp(spacetime.psi, xE, yE);
    double dRE  = spacetime.get_ddx_interp(spacetime.psi, xE, yE) * sin(sigma[iE])
                + spacetime.get_ddy_interp(spacetime.psi, xE, yE) * cos(sigma[iE]);
    g_END = (f[iE] + 2.*f[iE]*f[iE]*dRE/psiE) * sigma[0];
}

void AHFShooter::interval_bisection(const Spacetime& spacetime) {
    int iter_a = 501;
    int iter_b = 51;
    double rat = 0.98;
    double fb  = rmax;
    double fa  = rmax * rat;
    double fc  = 0., gc = 0.;

    // find a sensible bracket
    for (int i = 0; i < iter_a; ++i) {
        fc = fb;
        shoot(spacetime, fc);
        gc = g[num_points-1];
        if (OOB || gc > g_END) {
            fb = rat * fb;
        } else {
            fa = fb;
            fb = fb / rat;
            break;
        }
    }

    // bisect within bracket
    for (int j = 0; j < iter_b; ++j) {
        fc = 0.5 * (fa + fb);
        shoot(spacetime, fc);
        gc = g[num_points-1];
        if (gc > g_END) fb = fc;
        else            fa = fc;
    }

    shoot(spacetime, fc);
}

void AHFShooter::hello() const {
    std::cout << "Hello from AHFShooter!\n"
              << " - Number of horizon points = " << num_points << "\n";
}

