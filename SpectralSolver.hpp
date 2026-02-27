#ifndef SPECTRALSOLVER_HPP
#define SPECTRALSOLVER_HPP

#include <vector>

class Spacetime;

class SpectralSolver {
public:
    // 1D horizon representation
    std::vector<double> sigma;     // parameter / angular coordinate
    std::vector<double> f;         // horizon radius
    std::vector<double> dfds;         // legendre calculated deriv
    std::vector<double> W;         // curvature
    std::vector<double> psi; // psi evaluated on horizon
    std::vector<double> dpsi_dR; // spherical radial derivative (not cylindrical radius)
    std::vector<double> dW_dR; // spherical radial derivative (not cylindrical radius)
    std::vector<double> dfdt; // for relaxation time updates

    // points 
    int num_points;
    int N; // number of overtones
    std::vector<double> coeff; // basis expansion coefficients
    std::vector<double> J; // jacobean for basis coeffs
    double dt; // relaxation time stepper
    double ds; // delta parameter/angle

    SpectralSolver(int npoints, int Nspec);


    void initialize(const Spacetime& spacetime, const double f0);

    // best way to find coeffs
    void gradient_descent(const Spacetime& spacetime);

    // area infinitessimal on surface
    double dA(const int i);

    // returns geometric horizon area
    double area();

    // returns flat metric area
    double area_flat();

    // irreducible horizon mass
    double mass_irr();
    // full mass 
    double mass();

    // horizon angular momentum
    double J_H();
    double a_H(); // normalised by M
    double chi_H(); // normalised by M^2

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

    void residual();
    void refresh(const Spacetime& spacetime);
    void refresh_pureBH();

    void save(const std::string &filename);

    void hello() const;
};

#endif

