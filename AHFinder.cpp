#include "AHFinder.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>

AHFinder::AHFinder(int npoints) {
    sigma.resize(npoints);
    f.resize(npoints);
    psi.resize(npoints);
    num_points = npoints;
}

void AHFinder::initialize(const Spacetime& spacetime) {
    double x = 0., y = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        double d_sigma = 0.5 * Params::pi / num_points;
        // cell centred
        sigma[i] = (i + 1./2.) * d_sigma;
        x = f[i] * sin(sigma[i]);
        y = f[i] * cos(sigma[i]);
        f[i] = 1.0;           // initial guess for horizon radius
        psi[i] = spacetime.get_val_interp(spacetime.psi,x,y);
    }
}

// rough area
double AHFinder::area() {
    double d_sigma = 0.5 * Params::pi / num_points;
    double A = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = psi[i] * sqrt(f[i]*f[i] + dfds*dfds);
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * d_sigma;
    }
    // two because of reflective symmetry
    return 2.*A;
}

void AHFinder::update(const Spacetime& spacetime) {
    // Placeholder for expansion/minimization logic
    // code here 
}

void AHFinder::hello() const {
    std::cout << "Hello from AHFinder!\n"
              << " - Number of horizon points = "
              << num_points << std::endl;
}

