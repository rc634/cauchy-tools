#pragma once
#include <vector>
#include <string>

class Spacetime;

class Surface {
public:
    int num_points;
    double ds;
    double dt;
    std::vector<double> sigma;    // angular coord [0, pi/2]
    std::vector<double> f;        // spherical radius at each sigma
    std::vector<double> psi;
    std::vector<double> W;
    std::vector<double> dpsi_dR;  // spherical radial derivative
    std::vector<double> dW_dR;    // spherical radial derivative
    std::vector<double> dfdt;

    const Spacetime* st = nullptr;  // set by refresh(); null in flat-space mode

    Surface(int npoints);

    void announce(const std::string& label) const;
    void initialize(double f0 = 1.0);           // sphere of radius f0
    void initialize(double a, double b);         // ellipsoid: a=polar, b=equatorial

    // stencil derivatives
    double d(const std::vector<double>& field, int i);
    double d2(const std::vector<double>& field, int i);

    // sample spacetime fields onto surface
    void refresh(const Spacetime& spacetime);

    // area element
    double dA(int i);

    // geometric/physical properties
    double area();
    double area_flat();
    double psi_avg();
    double mass_irr();
    double mass();
    double J_H();
    double a_H();
    double chi_H();
    double mass_MS();
    double mass_SC();
    double r();
    double coord_r();
    double psi_h();
    double eccentricity();
    double coord_eccentricity();
    double res();

    double polar_circumference();
    double equatorial_circumference();
    double proper_height();   // integrates psi^2 dz if spacetime set, else returns f[0]
    double proper_width();    // integrates psi^2 dx if spacetime set, else returns f[end]
    double ramanujan_perimeter(); // flat ellipse perimeter approximation

    bool   inside(double x, double y) const; // is (x,y) inside the surface?
    double dV(double x, double y);           // proper volume element at (x,y), includes 2pi and grid cell area
    double V();                              // total proper volume inside the surface

    double div_s(int i); // covariant divergence of the outward unit normal s at grid point i
    double K0();         // area-weighted mean of div_s — mean curvature of the 2-surface

    void moi_surf(); // moment of inertia: surface integral of (psi^2 f sin sigma)^2 dA
    void moi_vol();  // moment of inertia: volume integral of (psi^2 x)^2 dV

    void save(const std::string& filename);
    void cout_state();
    void print() const;
};
