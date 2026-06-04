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

        psi[i] = spacetime.get_val_interp(spacetime.psi, x, y);
        W[i]   = spacetime.get_val_interp(spacetime.W,   x, y);

        dpsi_dR[i] = cos(TH) * spacetime.get_ddy_interp(spacetime.psi, x, y)
                   + sin(TH) * spacetime.get_ddx_interp(spacetime.psi, x, y);

        dW_dR[i]   = cos(TH) * spacetime.get_ddy_interp(spacetime.W, x, y)
                   + sin(TH) * spacetime.get_ddx_interp(spacetime.W, x, y);
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
        double p = st->get_val_interp(st->psi, 0., z);
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
        double p = st->get_val_interp(st->psi, x, 0.);
        L += p * p * dx;
    }
    return L;
}

// Ramanujan approximation for the perimeter of an ellipse with
// semi-axes H (polar-height) and W (equatorial-width), in curved space
double Surface::ramanujan_perimeter() {
    double H = proper_height();
    double W = proper_width();
    return Params::pi * (3.*(H+W) - sqrt((3.*H+W)*(H+3.*W)));
}

// returns true if the point (x,y) lies inside the surface f(sigma)
// sigma is interpolated linearly from the discrete f array
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
    double p = (st != nullptr) ? st->get_val_interp(st->psi, x, y) : 1.0;
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
        double dpsi_dx  = st->get_ddx_interp(st->psi, x, y);
        double dpsi_dy  = st->get_ddy_interp(st->psi, x, y);
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
                double p  = (st != nullptr) ? st->get_val_interp(st->psi, x, y) : 1.0;
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

void Surface::save(const std::string& filename) {
    std::ofstream out("data/" + filename + ".dat");
    if (!out.is_open())
        throw std::runtime_error("Could not open file: " + filename);
    out << std::scientific << std::setprecision(Params::save_precision);
    for (int i = 0; i < num_points; ++i)
        out << f[i] << "," << sigma[i] << "\n";
    out.close();
}

void Surface::cout_state() {
    //std::cout << std::scientific << std::setprecision(Params::cout_precision);
    // std::cout << std::setprecision(Params::cout_precision);
    // std::cout << "\n* ============ 2-Surface ============== *\n";
    std::cout << "\n* ========== Proper Geometry ========== *\n";
    std::cout << "| Area               = " << area() << "\n";
    std::cout << "| Height             = " << proper_height() << "\n";
    std::cout << "| Width              = " << proper_width() << "\n";
    std::cout << "| Eccentricity       = " << eccentricity() << "\n";
    std::cout << "| Polar circ         = " << polar_circumference() << "\n";
    std::cout << "| Equatorial circ    = " << equatorial_circumference() << "\n";
    std::cout << "| Volume             = " << V() << "\n";
    std::cout << "| Mean Curvature     = " << K0() << "\n";
    std::cout << "|\n* ============= Physics =============== *\n";
    std::cout << "| Irreducible Mass   = " << mass_irr() << "\n";
    std::cout << "| Dimensionless Spin = " << chi_H() << "\n";
    std::cout << "| Psi Avg            = " << psi_avg() << "\n";
    std::cout << "|\n* ========== Other Geometry =========== *\n";
    std::cout << "| PA Radius          = " << r() << "\n";
    std::cout << "| Coord Radius Avg   = " << coord_r() << "\n";
    std::cout << "| Coord Eccentricity = " << coord_eccentricity() << "\n";
    std::cout << "| Ramanujan perim    = " << ramanujan_perimeter() << "\n";
    std::cout << "* ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ *\n";
    std::cout << std::defaultfloat;
}

void Surface::print() const {
    std::cout << "# sigma  f\n";
    for (int i = 0; i < num_points; ++i)
        std::cout << sigma[i] << "  " << f[i] << "\n";
}
