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

        // cell centred
        sigma[i] = (i + 1./2.) * ds;
        f[i] = 1.7;           // initial guess for horizon radius
        dfdt[i] = 0.;

        // coords for interpolation
        x = f[i] * sin(sigma[i]);
        y = f[i] * cos(sigma[i]);
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

    return (-f5 + 16.*f4 + -30.*f3 + 16.*f2 - f1)/(12. * ds * ds);
}

// area
double AHFinder::area() {
    double A = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = pow(psi[i],8) * sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
    }
    // two because of reflective symmetry
    return 2.*A;
}

// rough mass
double AHFinder::mass() {
    double A = 0.;
    double dA = 0.;
    double MA = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = pow(psi[i],8) * sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        dA = 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
        A += dA;
        MA += dA * (psi[i]-1.)*2.*f[i];
    }
    // two because of reflective symmetry
    return MA/A;
}

// rough psi
double AHFinder::psi_h() {
    double A = 0.;
    double dA = 0.;
    double psiA = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = pow(psi[i],8) * sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        dA = 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
        A += dA;
        psiA += dA * psi[i];
    }
    // two because of reflective symmetry
    return psiA/A;
}

// rough r
double AHFinder::r() {
    double A = 0.;
    double dA = 0.;
    double rA = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = pow(psi[i],8) * sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        dA = 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
        A += dA;
        rA += dA * f[i];
    }
    // two because of reflective symmetry
    return rA/A;
}

void AHFinder::update(const Spacetime& spacetime) {
    // Placeholder for expansion/minimization logic
    // code here 
}

// to be called after we update the surface
// update external fields such as psi[f,s] 
void AHFinder::refresh(const Spacetime& spacetime) {
    double x = 0., y = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // update psi with new coords f,s
        x = f[i] * sin(sigma[i]);
        y = f[i] * cos(sigma[i]);
        psi[i] = spacetime.get_val_interp(spacetime.psi,x,y);
    }

}

void AHFinder::relax() {
    // Poisson solver towards correct surface
    double term1 = 0.;
    double term2 = 0.;
    double term3 = 0.;
    double dfds = 0.;
    double dpsi_dr = 0.;
    double r = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // dfdt[i] = - 0.2 * f[i];
        dfds = d(f,i);
        r = f[i];
        //
        term1 = (r*r + dfds*dfds)/(r*r + 2.*dfds*dfds);
        //
        term2 = (8. * dfds * dpsi_dr / psi[i] 
                    - 8. *d (psi,i) / psi[i] 
                    - dfds * cos(sigma[i])/sin(sigma[i]));
        //
        term3 = (1./r)*(2.*r*r + dfds*dfds)/(r*r + dfds*dfds);
        //
        dfdt[i] = d2(f,i) - (term1*term2 + term3);
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

