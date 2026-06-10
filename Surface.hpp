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
    double mass_ADM();
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
    double ramanujan_eccentricity(); // eccentricity from Ramanujan perimeter inversion

    bool   inside(double x, double y) const; // is (x,y) inside the surface?
    double dV(double x, double y);           // proper volume element at (x,y), includes 2pi and grid cell area
    double V();                              // total proper volume inside the surface

    double integrate_src();           // integral of src * dV over full grid
    double integrate_src_enclosed();  // integral of src * dV only inside the surface

    double integrate_J();             // integral of src * srcW * x^2 * dV over full grid
    double integrate_J_enclosed();    // integral of src * srcW * x^2 * dV only inside the surface

    double integrate_src_flat();           // integral of src * dV_flat (no psi^6) over full grid
    double integrate_src_enclosed_flat();  // integral of src * dV_flat only inside the surface

    double integrate_J_flat();             // integral of src * srcW * x^2 * dV_flat over full grid
    double integrate_J_enclosed_flat();    // integral of src * srcW * x^2 * dV_flat only inside the surface

    double integrate_src_eff();           // integral of src_eff * dV over full grid
    double integrate_src_enclosed_eff();  // integral of src_eff * dV only inside the surface

    double integrate_J_eff();             // integral of src_eff * srcW_eff * x^2 * dV over full grid
    double integrate_J_enclosed_eff();    // integral of src_eff * srcW_eff * x^2 * dV only inside the surface

    double integrate_src_eff_flat();           // integral of src_eff * dV_flat over full grid
    double integrate_src_enclosed_eff_flat();  // integral of src_eff * dV_flat only inside the surface

    double integrate_J_eff_flat();             // integral of src_eff * srcW_eff * x^2 * dV_flat over full grid
    double integrate_J_enclosed_eff_flat();    // integral of src_eff * srcW_eff * x^2 * dV_flat only inside the surface

    // 2-metric theta theta
    double A(int i);
    double Ap(int i);   // dA/dσ
    // 2-metric phi phi
    double B(int i);
    double Bp(int i);   // dB/dσ
    double Bpp(int i);  // d²B/dσ²
    // 2d ricci scalar
    double Ricci(int i);
    // area weighted mean
    double Ricci_mean();
    // area-weighted Y_20 projection of Ricci scalar
    double RY20();
    // eccentricity estimate from RY20: e^2 ~ (3*sqrt(5*pi)/2) * |RY20|
    double Ricci_ecc();

    double div_s(int i); // covariant divergence of the outward unit normal s at grid point i
    double K0();         // area-weighted mean of div_s — mean curvature of the 2-surface
    double K2();

    void moi_surf(); // moment of inertia: surface integral of (psi^2 f sin sigma)^2 dA
    void moi_vol();  // moment of inertia: volume integral of (psi^2 x)^2 dV

    void smooth(int hw = 5);        // moving-average smooth of f, then refresh
    void smooth_edges(int hw = 5);  // tapered 3-point smooth near pole and equator only

    void save(const std::string& filename);
    void cout_state();
    void print() const;
};
