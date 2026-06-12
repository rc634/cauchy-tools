#include "Surface.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

Surface::Surface(int npoints) : num_points(npoints) {
    ds = 0.5 * Params::pi / npoints;
    dt = Params::CFL * ds * ds;
    sigma.resize(npoints, 0.0);
    f.resize(npoints, 0.0);
    psi.resize(npoints, 0.0);
    W.resize(npoints, 0.0);
    dpsi_dR.resize(npoints, 0.0);
    dW_dR.resize(npoints, 0.0);
    dfdt.resize(npoints, 0.0);
}

void Surface::announce(const std::string& label) const {
    int w = label.size() + 4;
    std::string bar(w, '=');
    std::cout << "\n*" << bar << "*\n";
    std::cout << "| " << label << "   |\n";
    std::cout << "*" << bar << "*\n\n";
}

void Surface::initialize(double f0) {
    st = nullptr;
    for (int i = 0; i < num_points; ++i) {
        sigma[i]    = (i + 0.5) * ds;   // cell-centred angular grid
        f[i]        = f0;
        dfdt[i]     = 0.;
        psi[i]      = 1.;               // flat space conformal factor
        W[i]        = 0.;               // no rotation / extrinsic curvature
        dpsi_dR[i]  = 0.;
        dW_dR[i]    = 0.;
    }
    // std::cout << " - Surface::initialize, f0 = " << f0 << std::endl;
}

// ellipsoid: a = polar radius (height), b = equatorial radius (width)
// f(sigma) = ab / sqrt(a^2*sin^2 + b^2*cos^2)
void Surface::initialize(double a, double b) {
    initialize();   // set sigma grid and flat-space field defaults
    for (int i = 0; i < num_points; ++i)
        f[i] = a * b / sqrt(a*a*sin(sigma[i])*sin(sigma[i]) + b*b*cos(sigma[i])*cos(sigma[i]));
}

double Surface::d(const std::vector<double>& field, int i) {
    int i1 = i - 1;
    int i3 = i + 1;
    if (i == 0) i1 = 0;
    else if (i == num_points - 1) i3 = num_points - 1;
    return (field[i3] - field[i1]) / (2. * ds);
}

double Surface::d2(const std::vector<double>& field, int i) {
    int i1 = i - 1;
    int i3 = i + 1;
    if (i == 0) i1 = 0;
    else if (i == num_points - 1) i3 = num_points - 1;
    return (field[i3] - 2. * field[i] + field[i1]) / (ds * ds);
}

void Surface::refresh(const Spacetime& spacetime) {
    st = &spacetime;
    double x = 0., y = 0.;
    double TH = 0.;
    for (int i = 0; i < num_points; ++i) {
        TH = sigma[i];
        x = f[i] * sin(TH);
        y = f[i] * cos(TH);

        psi[i] = spacetime.get_val_bicubic(spacetime.psi, x, y);
        W[i]   = spacetime.get_val_bicubic(spacetime.W,   x, y);

        dpsi_dR[i] = cos(TH) * spacetime.get_ddy_bicubic(spacetime.psi, x, y)
                   + sin(TH) * spacetime.get_ddx_bicubic(spacetime.psi, x, y);

        dW_dR[i]   = cos(TH) * spacetime.get_ddy_bicubic(spacetime.W, x, y)
                   + sin(TH) * spacetime.get_ddx_bicubic(spacetime.W, x, y);
    }
}

double Surface::dA(int i) {
    return 2. * Params::pi * pow(psi[i], 4)
             * sqrt(f[i] * f[i] + d(f, i) * d(f, i))
             * f[i] * sin(sigma[i]) * ds;
}

double Surface::area() {
    double A = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
    }
    return 2. * A;
}

// surface average of psi
double Surface::psi_avg() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * psi[i];
    }
    return B / A;
}

double Surface::area_flat() {
    double A = 0.;
    for (int i = 0; i < num_points; ++i) {
        double geom = sqrt(f[i] * f[i] + d(f, i) * d(f, i));
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
    }
    return 2. * A;
}

double Surface::mass_irr() {
    return 0.25 * sqrt(area() / Params::pi);
}

double Surface::mass() {
    double M_irr = mass_irr();
    double J = J_H();
    return sqrt(M_irr * M_irr + 0.25 * J * J / M_irr / M_irr);
}

double Surface::J_H() {
    double J = 0.;
    for (int i = 0; i < num_points; ++i) {
        J += 0.5 * (
            f[i] * (W[i] + f[i] * dW_dR[i]) * sin(sigma[i])
            - d(f, i) * (W[i] * cos(sigma[i]) + sin(sigma[i]) * d(W, i))
        ) * f[i] * sin(sigma[i]) * ds;
    }
    return 2. * J;
}

double Surface::a_H() {
    return J_H() / mass();
}

double Surface::chi_H() {
    return a_H() / mass();
}

double Surface::mass_MS() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * (psi[i] - 1.) * 2. * f[i];
    }
    return B / A;
}

double Surface::mass_SC() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * (1. - pow(psi[i], -4)) * f[i] / 2.;
    }
    return B / A;
}

double Surface::mass_ADM() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * dpsi_dR[i];
    }
    return -2. * coord_r() * coord_r() * B / A;
}

double Surface::r() {
    // polar areal radius
    return sqrt(area()/(4.*M_PI));
}

double Surface::coord_r() {
    // average coordinate radius
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * f[i];
    }
    return B / A;
}

double Surface::psi_h() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * psi[i];
    }
    return B / A;
}

double Surface::eccentricity() {
    double f0 = proper_height();
    double fE = proper_width();
    if (fE == f0) return 0.;
    else if (fE > f0) return sqrt(1. - pow(f0 / fE, 2));
    else              return sqrt(1. - pow(fE / f0, 2));
}

double Surface::coord_eccentricity() {
    double f0 = f[0];
    double fE = f[num_points - 1];
    if (fE == f0) return 0.;
    else if (fE > f0) return sqrt(1. - pow(f0 / fE, 2));
    else              return sqrt(1. - pow(fE / f0, 2));
}

double Surface::res() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * dfdt[i];
    }
    return B / A;
}

// polar great circle: N pole -> equator -> S pole -> equator -> N pole
// dl = psi^2 * sqrt(f^2 + f'^2) dsigma, factor 4 for full loop
double Surface::polar_circumference() {
    double C = 0.;
    for (int i = 0; i < num_points; ++i)
        C += pow(psi[i], 2) * sqrt(f[i]*f[i] + d(f,i)*d(f,i)) * ds;
    return 4. * C;
}

// azimuthal circle at the equator (sigma=pi/2, x=f_eq)
// dl = psi^2 * f_eq dphi, integrated over 2pi
double Surface::equatorial_circumference() {
    return 2. * Params::pi * pow(psi[num_points-1], 2) * f[num_points-1];
}

// integrate psi^2 dz from origin to the pole along the z-axis (x=0)
// falls back to coordinate value f[0] if no spacetime is loaded
double Surface::proper_height() {
    if (st == nullptr) return f[0];
    double z_top = f[0];
    double dz = z_top / num_points;
    double L = 0.;
    for (int i = 0; i < num_points; ++i) {
        double z = (i + 0.5) * dz;
        double p = st->get_val_bicubic(st->psi, 0., z);
        L += p * p * dz;
    }
    return L;
}

// integrate psi^2 dx from origin to the equator along the x-axis (y=0)
// falls back to coordinate value f[end] if no spacetime is loaded
double Surface::proper_width() {
    if (st == nullptr) return f[num_points-1];
    double x_eq = f[num_points-1];
    double dx = x_eq / num_points;
    double L = 0.;
    for (int i = 0; i < num_points; ++i) {
        double x = (i + 0.5) * dx;
        double p = st->get_val_bicubic(st->psi, x, 0.);
        L += p * p * dx;
    }
    return L;
}

// Recover effective ellipse semi-axes from Ramanujan perimeter inversion,
// then return the eccentricity sqrt(1 - (min/max)^2)
double Surface::ramanujan_eccentricity() {
    double pol = polar_circumference();
    double eqa = equatorial_circumference();
    double a = (sqrt(-5.*eqa*eqa + 6.*eqa*pol + 3.*pol*pol) - 2.*eqa + 3.*pol) / (6. * Params::pi);
    double b = eqa / (2. * Params::pi);
    if (a == b) return 0.;
    double mn = std::min(a, b);
    double mx = std::max(a, b);
    return sqrt(1. - (mn*mn) / (mx*mx));
}

// returns true if the point (x,y) lies inside the surface f(sigma)
// sigma is bicubicolated linearly from the discrete f array
bool Surface::inside(double x, double y) const {
    double R = sqrt(x*x + y*y);
    if (R == 0.) return true;
    double sig = atan2(x, y);   // angle from z-axis: 0 at pole, pi/2 at equator
    if (sig < 0. || sig > 0.5 * Params::pi) return false;
    // continuous index into cell-centred grid: sigma[i] = (i+0.5)*ds
    double ci = sig / ds - 0.5;
    int i0 = (int)ci;
    if (i0 < 0) return R < f[0];
    if (i0 >= num_points - 1) return R < f[num_points - 1];
    double t = ci - i0;
    double f_sig = (1. - t) * f[i0] + t * f[i0 + 1];
    return R < f_sig;
}

// proper volume element at (x,y):
//   dV = psi^6 * 2*pi * x * dx * dy
// psi^6 = sqrt(det gamma) for the isotropic 3-metric psi^4 * delta_ij
// 2*pi comes from integrating the azimuthal angle
// falls back to psi=1 (flat space) if no spacetime is loaded
double Surface::dV(double x, double y) {
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double p = (st != nullptr) ? st->get_val_bicubic(st->psi, x, y) : 1.0;
    return pow(p, 6) * 2. * Params::pi * x * ddx * ddy;
}

// total proper volume inside the surface
// loops over the spacetime grid (first quadrant x>=0, y>=0),
// accumulates dV for interior points, then doubles for southern hemisphere
double Surface::V() {
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double vol = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y)) vol += dV(x, y);
        }
    }
    return 2. * vol;   // southern hemisphere by reflective symmetry
}

double Surface::integrate_src() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src, i, j) * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_src_enclosed() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src, i, j) * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_J() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src, i, j) * st->get_val(st->srcW, i, j) * x * x * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_J_enclosed() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src, i, j) * st->get_val(st->srcW, i, j) * x * x * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_src_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src, i, j) * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_src_enclosed_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src, i, j) * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_J_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src, i, j) * st->get_val(st->srcW, i, j) * x * x * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_J_enclosed_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src, i, j) * st->get_val(st->srcW, i, j) * x * x * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_src_eff() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src_eff, i, j) * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_src_enclosed_eff() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src_eff, i, j) * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_J_eff() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src_eff, i, j) * st->get_val(st->srcW_eff, i, j) * x * x * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_J_enclosed_eff() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src_eff, i, j) * st->get_val(st->srcW_eff, i, j) * x * x * dV(x, y);
        }
    }
    return 2. * total;
}

double Surface::integrate_src_eff_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src_eff, i, j) * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_src_enclosed_eff_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src_eff, i, j) * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_J_eff_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            total += st->get_val(st->src_eff, i, j) * st->get_val(st->srcW_eff, i, j) * x * x * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

double Surface::integrate_J_enclosed_eff_flat() {
    if (st == nullptr) return 0.;
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double total = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y))
                total += st->get_val(st->src_eff, i, j) * st->get_val(st->srcW_eff, i, j) * x * x * 2. * Params::pi * x * ddx * ddy;
        }
    }
    return 2. * total;
}

// the sqrt of theta theta component of the 2 metric
double Surface::A(int i) {
    double q = d(f, i) / f[i];
    return psi[i] * psi[i] * f[i] * sqrt(1. + q * q);
}

// the sqrt of phi phi component of the 2 metric
double Surface::B(int i) {
    return psi[i] * psi[i] * f[i] * sin(sigma[i]);
}

// dA/dσ = d(ψ²fN)/dσ, analytic chain rule — no nested stencils
double Surface::Ap(int i) {
    double p   = psi[i],  pp  = d(psi, i);
    double f_  = f[i],    fp  = d(f, i),  fpp = d2(f, i);
    double q   = fp / f_;
    double N   = sqrt(1. + q*q);
    double qp  = fpp/f_ - q*q;
    return 2.*p*pp*f_*N  +  p*p*fp*N  +  p*p*f_*(q*qp/N);
} // p p f N

// dB/dσ = d(ψ²f sinσ)/dσ, analytic chain rule — no nested stencils
double Surface::Bp(int i) {
    double p  = psi[i],  pp = d(psi, i);
    double f_ = f[i],    fp = d(f, i);
    double TH = sigma[i];
    return 2.*p*pp*f_*sin(TH)  +  p*p*fp*sin(TH)  +  p*p*f_*cos(TH);
} // p p f sin(th)

// d²B/dσ² = d²(ψ²f sinσ)/dσ², analytic — no nested stencils
double Surface::Bpp(int i) {
    double p   = psi[i],  pp  = d(psi, i),  ppp = d2(psi, i);
    double f_  = f[i],    fp  = d(f, i),    fpp = d2(f, i);
    double TH  = sigma[i];
    double sinTH = sin(TH), cosTH = cos(TH);
    double s_term = 2.*(pp*pp + p*ppp)*f_ + 4.*p*pp*fp + p*p*fpp - p*p*f_;
    double c_term = 4.*p*pp*f_ + 2.*p*p*fp;
    return s_term * sinTH  +  c_term * cosTH;
} // d2 ( p p f sin(th) )

// 2d ricci scalar at i
// metric A^2 d th^2 + B^2 d phi^2
// 2-Ricci is 2/A^2 * [A'/A B'/B - B''/B]
double Surface::Ricci(int i) {
    double A_C  = A(i);
    double B_C  = B(i);
    double A_p  = Ap(i);
    double B_p  = Bp(i);
    double B_pp = Bpp(i);

    return 2. / (A_C * A_C) * (A_p / A_C * B_p / B_C - B_pp / B_C);
}

// covariant divergence of the outward unit normal s^i at grid point i
// in the physical metric gamma_ij = psi^4 delta_ij
double Surface::div_s(int i) {
    double TH = sigma[i];
    double x  = f[i] * sin(TH);
    double y  = f[i] * cos(TH);
    double q =  d(f,i) / f[i];
    double dqdr = -q/f[i];
    double dqdth = d2(f,i)/f[i] - q*q;
    double N    = sqrt(1. + q*q);
    double N3   = pow(1. + q*q, 1.5);

    // sqrt of 3-metric determinant: sqrt(gamma) = psi^6 R^2 sinθ
    // kept for conventions
    // double det3 = pow(psi[i], 6) * f[i]*f[i] * sin(TH);

    // spherical coordinate derivatives of psi from cylindrical spacetime data,
    // following the same pattern as refresh()
    double dpsi_dr = 0., dpsi_dth = 0.;
    if (st != nullptr) {
        double dpsi_dx  = st->get_ddx_bicubic(st->psi, x, y);
        double dpsi_dy  = st->get_ddy_bicubic(st->psi, x, y);
        dpsi_dr  = sin(TH) * dpsi_dx + cos(TH) * dpsi_dy;
        dpsi_dth = f[i] * (cos(TH) * dpsi_dx - sin(TH) * dpsi_dy);
    }

    // unit vector 
    double s_r = pow(psi[i],-2)/sqrt(1+q*q);
    double s_th = - (q / f[i]) * s_r;

    // s^i d_i ln(det3)
    double s_grad_ln_det3 = s_r * (6. * dpsi_dr / psi[i] + 2. / f[i]);
    s_grad_ln_det3 += s_th * (6. * dpsi_dth / psi[i] + cos(TH)/sin(TH) );

    // now calculate the partial_i s^i part
    // ∂_R s^R = -2ψ⁻³ ∂_R ψ / N  +  ψ⁻² q² / (f N³)
    double d_r_s_r  = - 2.*dpsi_dr * pow(psi[i],-3) / N
                      - pow(psi[i],-2) * q * dqdr / N3;
    // ∂_θ s^θ = 2ψ⁻³ ∂_θψ q/(fN)  -  ψ⁻² f''/(f²N³)
    double d_th_s_th = 2.*pow(psi[i],-3) * dpsi_dth * q / (f[i] * N)
                       - pow(psi[i],-2) * dqdth / (f[i] * N )
                       + pow(psi[i],-2) * q * q * (f[i] + d2(f,i)) / (f[i]*f[i] * N3);
    // combine the two pieces
    double div_s_flat = d_r_s_r + d_th_s_th;

    return div_s_flat + s_grad_ln_det3;
}

double Surface::K0() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * div_s(i);
    }
    return B / A;
}

double Surface::K2() {
    constexpr int trim = 5;

    // 5-point local average of Ricci, computed over interior points only
    auto R_av = [&](int i) {
        return (Ricci(i-2) + Ricci(i-1) + Ricci(i) + Ricci(i+1) + Ricci(i+2)) / 5.0;
    };

    // area-weighted mean of R_av over trimmed range
    double A = 0., mean = 0.;
    for (int i = trim; i < num_points - trim - 1; ++i) {
        A    += dA(i);
        mean += dA(i) * R_av(i);
    }
    mean /= A;

    // area-weighted variance of relative deviation (matches Python vis script)
    double B = 0.;
    for (int i = trim; i < num_points - trim - 1; ++i) {
        double dR = 1.0 - R_av(i) / mean;
        B += dA(i) * dR * dR;
    }
    return B / A;
}

double Surface::Ricci_mean() {
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        A += dA(i);
        B += dA(i) * Ricci(i);
    }
    return B / A;
}

double Surface::RY20() {
    // Y_2^0(σ) = sqrt(5/4π) * (3cos²σ - 1)/2
    // returns area-weighted projection: (1/A) ∫ Ricci * Y_20 dA
    constexpr double prefactor = 0.31539156525; // sqrt(5/(4*pi))
    double A = 0., B = 0.;
    for (int i = 0; i < num_points; ++i) {
        double c = cos(sigma[i]);
        double Y = prefactor * (3.*c*c - 1.) * 0.5;
        A += dA(i);
        B += dA(i) * Ricci(i) * Y;
    }
    return (B / A) / Ricci_mean();
}

double Surface::Ricci_ecc() {
    // e^2 = (3*sqrt(5*pi)/2) * |RY20|, derived for small deformations from a sphere
    constexpr double coeff = 3. * 3.96332729760601 / 2.;  // 3*sqrt(5*pi)/2 ~ 5.945
    return sqrt(coeff * std::abs(RY20()));
}

// moment of inertia tensor (surface integral version)
// axisymmetry => principal axes, so I = diag(I_xx, I_xx, I_zz)
// I_zz = integral (psi^2 f sinσ)^2 dA
// I_xx = (1/2)*I_zz + integral (psi^2 f cosσ)^2 dA
void Surface::moi_surf() {
    double Izz = 0., Izz_extra = 0.;
    for (int i = 0; i < num_points; ++i) {
        double p2  = psi[i] * psi[i];
        double rp  = p2 * f[i] * sin(sigma[i]);   // proper equatorial radius
        double zp  = p2 * f[i] * cos(sigma[i]);   // proper axial coord
        Izz       += rp * rp * dA(i);
        Izz_extra += zp * zp * dA(i);
    }
    Izz       *= 2.;
    Izz_extra *= 2.;
    double Ixx = 0.5 * Izz + Izz_extra;
    std::cout << std::scientific << std::setprecision(Params::cout_precision);
    std::cout << "| MOI_surf Ixx  = " << Ixx << "\n";
    std::cout << "| MOI_surf Izz  = " << Izz << "\n";
    std::cout << std::defaultfloat;
}

// moment of inertia tensor (volume integral version)
// axisymmetry => principal axes, so I = diag(I_xx, I_xx, I_zz)
// I_zz = integral (psi^2 x)^2 dV
// I_xx = (1/2)*I_zz + integral (psi^2 y)^2 dV
void Surface::moi_vol() {
    constexpr double ddx = (Params::XU - Params::XL) / Params::NX;
    constexpr double ddy = (Params::YU - Params::YL) / Params::NY;
    double Izz = 0., Izz_extra = 0.;
    for (int i = 0; i < Params::NX; ++i) {
        double x = Params::XL + (i + 0.5) * ddx;
        for (int j = 0; j < Params::NY; ++j) {
            double y = Params::YL + (j + 0.5) * ddy;
            if (inside(x, y)) {
                double p  = (st != nullptr) ? st->get_val_bicubic(st->psi, x, y) : 1.0;
                double p2 = p * p;
                double dv = dV(x, y);
                Izz       += (p2 * x) * (p2 * x) * dv;
                Izz_extra += (p2 * y) * (p2 * y) * dv;
            }
        }
    }
    Izz       *= 2.;
    Izz_extra *= 2.;
    double Ixx = 0.5 * Izz + Izz_extra;
    std::cout << std::scientific << std::setprecision(Params::cout_precision);
    std::cout << "| MOI_vol  Ixx  = " << Ixx << "\n";
    std::cout << "| MOI_vol  Izz  = " << Izz << "\n";
    std::cout << std::defaultfloat;
}

void Surface::smooth_edges(int hw) {
    std::vector<double> fc = f;   // read from copy so points don't contaminate each other

    for (int i = 0; i < hw; ++i) {
        double w   = 1.0 - static_cast<double>(i) / hw;  // 1 at edge, 0 at hw
        int    im  = (i == 0) ? 0 : i - 1;               // ghost: f[-1] = f[0]
        double avg = (fc[im] + fc[i] + fc[i + 1]) / 3.0;
        f[i] = (1.0 - w) * fc[i] + w * avg;
    }

    for (int i = num_points - hw; i < num_points; ++i) {
        double w   = 1.0 - static_cast<double>(num_points - 1 - i) / hw;  // 0 at hw from end, 1 at last point
        int    ip  = (i == num_points - 1) ? num_points - 1 : i + 1;      // ghost: f[N] = f[N-1]
        double avg = (fc[i - 1] + fc[i] + fc[ip]) / 3.0;
        f[i] = (1.0 - w) * fc[i] + w * avg;
    }
}

void Surface::smooth(int hw) {
    std::vector<double> f_smooth = f;
    for (int i = hw; i < num_points - hw; ++i) {
        double s = 0.;
        for (int j = -hw; j <= hw; ++j) s += f[i + j];
        f_smooth[i] = s / (2 * hw + 1);
    }
    f = f_smooth;
    if (st != nullptr) refresh(*st);
}

void Surface::save(const std::string& filename) {
    std::ofstream out(save_path + "/" + filename + ".dat");
    if (!out.is_open())
        throw std::runtime_error("Could not open file: " + filename);
    out << std::scientific << std::setprecision(Params::save_precision);
    out << "# theta  f  psi  div_s  Ricci\n";
    for (int i = 0; i < num_points; ++i)
        out << sigma[i] << "  " << f[i] << "  " << psi[i] << "  " << div_s(i)  << "  " << Ricci(i) << "\n";
    out.close();
}

void Surface::cout_state() {
    //std::cout << std::scientific << std::setprecision(Params::cout_precision);
    // std::cout << std::setprecision(Params::cout_precision);
    // std::cout << "\n* ============ 2-Surface ============== *\n";
    std::cout << "\n* ========== Surface Geometry ========== *\n";
    std::cout << "| Area               = " << area() << "\n";
    std::cout << "| Polar circ         = " << polar_circumference() << "\n";
    std::cout << "| Equatorial circ    = " << equatorial_circumference() << "\n";
    std::cout << "| Ramanujan ecc      = " << ramanujan_eccentricity() << "\n";
    std::cout << "|\n* =========== Curvatures ============== *\n";
    std::cout << "| Mean Div s         = " << K0() << "\n";
    std::cout << "| Mean 2-Ricci       = " << Ricci_mean() << "\n";
    std::cout << "| Ricci Deviation    = " << K2() << "\n";
    std::cout << "| RY20               = " << RY20() << "\n";
    std::cout << "| Ricci ecc          = " << Ricci_ecc() << "\n";
    std::cout << "|\n* ============= Physics =============== *\n";
    std::cout << "| Hawking Mass       = " << mass_irr() << "\n";
    std::cout << "| Dimensionless Spin = " << chi_H() << "\n";
    std::cout << "| Psi Avg            = " << psi_avg() << "\n";
    std::cout << "| 2R<psi-1>          = " << 2. * coord_r() * (psi_avg() - 1.) << "\n";
    std::cout << "| ADM Mass (dpsi/dR) = " << mass_ADM() << "\n";
    std::cout << "|\n* ========= Internal Geometry ========= *\n";
    std::cout << "| Height             = " << proper_height() << "\n";
    std::cout << "| Width              = " << proper_width() << "\n";
    std::cout << "| Eccentricity       = " << eccentricity() << "\n";
    std::cout << "| Volume             = " << V() << "\n";
    std::cout << "|\n* =========== Coord Stuff ============= *\n";
    std::cout << "| PA Radius          = " << r() << "\n";
    std::cout << "| Coord Radius Avg   = " << coord_r() << "\n";
    std::cout << "| Coord Eccentricity = " << coord_eccentricity() << "\n";
    std::cout << "|\n* ========== Matter Source ============ *\n";
    std::cout << "| Total Mass         = " << integrate_src() << "\n";
    std::cout << "| Mass Contained     = " << integrate_src_enclosed() << "\n";
    std::cout << "| Total AngMom       = " << integrate_J() << "\n";
    std::cout << "| Contained AngMom   = " << integrate_J_enclosed() << "\n";
    std::cout << "|\n* ======== Flat Source Integrals ======= *\n";
    std::cout << "| Total Mass (flat)  = " << integrate_src_flat() << "\n";
    std::cout << "| Mass Cont  (flat)  = " << integrate_src_enclosed_flat() << "\n";
    std::cout << "| Total J    (flat)  = " << integrate_J_flat() << "\n";
    std::cout << "| Contained J(flat)  = " << integrate_J_enclosed_flat() << "\n";
    std::cout << "|\n* ======== Effective Source ============ *\n";
    std::cout << "| Total Mass (eff)   = " << integrate_src_eff() << "\n";
    std::cout << "| Mass Cont  (eff)   = " << integrate_src_enclosed_eff() << "\n";
    std::cout << "| Total J    (eff)   = " << integrate_J_eff() << "\n";
    std::cout << "| Contained J(eff)   = " << integrate_J_enclosed_eff() << "\n";
    std::cout << "|\n* ====== Flat Effective Source ========= *\n";
    std::cout << "| Total Mass (eff,f) = " << integrate_src_eff_flat() << "\n";
    std::cout << "| Mass Cont  (eff,f) = " << integrate_src_enclosed_eff_flat() << "\n";
    std::cout << "| Total J    (eff,f) = " << integrate_J_eff_flat() << "\n";
    std::cout << "| Contained J(eff,f) = " << integrate_J_enclosed_eff_flat() << "\n";
    std::cout << "* ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ *\n";
    std::cout << std::defaultfloat;

    // verbose ricci printing
    // for (int i = 0; i < num_points; ++i)
    // std::cout << "i=" << i << " Ricci=" << Ricci(i) << "\n";

}

void Surface::print() const {
    std::cout << "# sigma  f\n";
    for (int i = 0; i < num_points; ++i)
        std::cout << sigma[i] << "  " << f[i] << "\n";
}
