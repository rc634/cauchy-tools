#include "AHFinder.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>
#include <cmath>

AHFinder::AHFinder(int npoints) {
    sigma.resize(npoints);
    f.resize(npoints);
    psi.resize(npoints);
    dfdt.resize(npoints);
    num_points = npoints;

    double dx = (Params::XU-Params::XL)/Params::NX;
    double dy = (Params::YU-Params::YL)/Params::NY;
    dt = Params::CFL * std::min(dx,dy);
    ds = 0.5 * Params::pi / num_points;
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
        dfdt[i] = 0.;
        psi[i] = spacetime.get_val_interp(spacetime.psi,x,y);
    }
}

// 1st derivative
double AHFinder::d(const std::vector<double> &field, int i) {
    // integer sample points
    int i1 = i-2;
    int i2 = i-1;
    int i3 = i;
    int i4 = i+1;
    int i5 = i+2;

    // symmetric boundaries
    if (i==0) {
        i2 = 0;
        i1 = 1;
    } 
    else if (i==1) {
        i1 = 0;
    }
    if (i==num_points-1) {
        i5 = num_points-2;
        i4 = num_points-1;
    } 
    else if (i==num_points-2) {
        i5 = num_points-1;
    }

    // sample field f_i(s_i)
    double f1 = field[i1];
    double f2 = field[i2];
    double f3 = field[i3];
    double f4 = field[i4];
    double f5 = field[i5];

    return (-f5 + 8.*f4 - 8.*f2 + f1)/(12. * ds);
}

// 1st derivative
double AHFinder::d2(const std::vector<double> &field, int i) {
    // integer sample points
    int i1 = i-2;
    int i2 = i-1;
    int i3 = i;
    int i4 = i+1;
    int i5 = i+2;

    // symmetric boundaries
    if (i==0) {
        i2 = 0;
        i1 = 1;
    } 
    else if (i==1) {
        i1 = 0;
    }
    if (i==num_points-1) {
        i5 = num_points-2;
        i4 = num_points-1;
    } 
    else if (i==num_points-2) {
        i5 = num_points-1;
    }

    // sample field f_i(s_i)
    double f1 = field[i1];
    double f2 = field[i2];
    double f3 = field[i3];
    double f4 = field[i4];
    double f5 = field[i5];

    return (-f5 + 16.*f4 + 30.*f3 - 16.*f2 - f1)/(12. * ds * ds);
}

// area
double AHFinder::area() {
    double d_sigma = 0.5 * Params::pi / num_points;
    double A = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = pow(psi[i],8) * sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * d_sigma;
    }
    // two because of reflective symmetry
    return 2.*A;
}

// rough mass
double AHFinder::mass() {
    double d_sigma = 0.5 * Params::pi / num_points;
    double A = 0.;
    double MA = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = pow(psi[i],8) * sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * d_sigma;
        MA += A * (psi[i]-1.)*2.*f[i];
    }
    // two because of reflective symmetry
    return MA/A;
}

void AHFinder::update(const Spacetime& spacetime) {
    // Placeholder for expansion/minimization logic
    // code here 
}

void AHFinder::relax() {
    // Poisson solver towards correct surface
    for (size_t i = 0; i < num_points; ++i) {
        dfdt[i] = - 0.01 * f[i];
    }
    for (size_t i = 0; i < num_points; ++i) {
        f[i] = f[i] + dfdt[i] * dt;
    }
}

void AHFinder::hello() const {
    std::cout << "Hello from AHFinder!\n"
              << " - Number of horizon points = "
              << num_points << std::endl;
}

