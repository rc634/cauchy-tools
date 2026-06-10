#include "EXSurf.hpp"
#include <iostream>

EXSurf::EXSurf(int npoints) : Surface(npoints) {}

void EXSurf::set_sphere(double r) {
    initialize(r);
}

void EXSurf::set_ellipsoid(double a, double b) {
    initialize(a, b);
}

void EXSurf::extraction_output() {
    std::cout << "\n* ========= Extraction Surface ======== *\n";
    std::cout << "| Proper Area        = " << area() << "\n";
    std::cout << "| Coord Radius Avg   = " << coord_r() << "\n";
    std::cout << "| Isotropic Mass     = " << 2. * coord_r() * (psi_avg() - 1.) << "\n";
    std::cout << "| ADM Mass (dpsi/dR) = " << mass_ADM() << "\n";
    std::cout << "* ====================================== *\n";
}
