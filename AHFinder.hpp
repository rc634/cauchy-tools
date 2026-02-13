#ifndef AHFINDER_HPP
#define AHFINDER_HPP

#include <vector>

class Spacetime;

class AHFinder {
public:
    // 1D horizon representation
    std::vector<double> sigma;     // parameter / angular coordinate
    std::vector<double> f;         // horizon radius
    std::vector<double> psi; // psi evaluated on horizon
    std::vector<double> dpsi_dr; // spherical radial derivative (not cylindrical radius)
    std::vector<double> dfdt; // for relaxation time updates

    // points 
    int num_points;
    double dt; // relaxation time stepper
    double ds; // delta parameter/angle

    AHFinder(int npoints);

    // area infinitessimal on surface
    double dA(const int i);

    // returns geometric horizon area
    double area();

    // returns flat metric area
    double area_flat();

    // irreducible horizon mass
    double mass();

    // Misner Sharp mass isotropic
    double mass_MS();

    // Misner Sharp mass sc gauge
    double mass_SC();

    // ODE error / residual
    double res();

    // approx horizon radius
    double r();

    // approx horizon psi
    double psi_h();

    // currently simple a/b, not root(1-aa/bb)
    double eccentricity();

    // partial derivs
    double d(const std::vector<double> &field, int i);
    double d2(const std::vector<double> &field, int i);

    void initialize(const Spacetime& spacetime, const double f0);
    void update(const Spacetime& spacetime);
    void relax();
    void refresh(const Spacetime& spacetime);

    void save(const std::string &filename);

    void hello() const;
};

#endif

