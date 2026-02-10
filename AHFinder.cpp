#include "AHFinder.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>

AHFinder::AHFinder(int npoints) {
    sigma.resize(npoints);
    f.resize(npoints);
    psi_sigma.resize(npoints);
    num_points = npoints;
}

void AHFinder::initialize() {
    for (size_t i = 0; i < num_points; ++i) {
        sigma[i] = 2.0 * Params::pi * i / num_points;
        f[i] = 1.0;           // initial guess for horizon radius
        psi_sigma[i] = 0.0;
    }
}

void AHFinder::update(const Spacetime& spacetime) {
    // Placeholder for expansion/minimization logic
    // code here 
}

void AHFinder::hello() const {
    std::cout << "Hello from AHFinder! "
              << "Number of horizon points = "
              << num_points << std::endl;
}

