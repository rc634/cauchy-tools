#include "AHFinder.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

AHFinder::AHFinder(int npoints) {
    sigma.resize(npoints);
    f.resize(npoints);
    psi.resize(npoints);
    W.resize(npoints);
    dpsi_dR.resize(npoints);
    dW_dR.resize(npoints);
    dfdt.resize(npoints);
    num_points = npoints;

    // double dx = (Params::XU-Params::XL)/Params::NX;
    // double dy = (Params::YU-Params::YL)/Params::NY;
    //dt = Params::CFL * std::min(dx,dy);
    ds = 0.5 * Params::pi / num_points;
    dt = Params::CFL * ds * ds;
}

void AHFinder::initialize(const Spacetime& spacetime, const double f0) {
    double x = 0., y = 0.;
    for (size_t i = 0; i < num_points; ++i) {

        // cell centred
        sigma[i] = (i + 1./2.) * ds;
        f[i] = f0;           // initial guess for horizon radius
        dfdt[i] = 0.;

        // update interpolated values
        // e.g. psi, dpsi_dR, dW_dR
        refresh(spacetime); 

    }

    std::cout << " - setting f[i] = " << f0 << std::endl;
}

// // 1st derivative // 4th order
// double AHFinder::d(const std::vector<double> &field, int i) {
//     // integer sample points
//     int i1 = i-2;
//     int i2 = i-1;
//     int i3 = i;
//     int i4 = i+1;
//     int i5 = i+2;

//     // symmetric boundaries
//     if (i==0) {
//         i2 = 0;
//         i1 = 1;
//     } 
//     else if (i==1) {
//         i1 = 0;
//     }
//     if (i==num_points-1) {
//         i5 = num_points-2;
//         i4 = num_points-1;
//     } 
//     else if (i==num_points-2) {
//         i5 = num_points-1;
//     }

//     // sample field f_i(s_i)
//     double f1 = field[i1];
//     double f2 = field[i2];
//     double f3 = field[i3];
//     double f4 = field[i4];
//     double f5 = field[i5];

//     return (-f5 + 8.*f4 - 8.*f2 + f1)/(12. * ds);
// }

// 1st derivative
double AHFinder::d(const std::vector<double> &field, int i) {
    // integer sample points
    int i1 = i-1;
    int i2 = i;
    int i3 = i+1;

    // symmetric boundaries i1 = 1; 
    if (i==0) {
        i1 = 0;
    }
    else if (i==num_points-1) {
        i3 = num_points-1;
    }

    // sample field f_i(s_i)
    double f1 = field[i1];
    double f2 = field[i2];
    double f3 = field[i3];

    return (f3-f1)/(2. * ds);
}

// // 2nd derivative - 4th order
// double AHFinder::d2(const std::vector<double> &field, int i) {
//     // integer sample points
//     int i1 = i-2;
//     int i2 = i-1;
//     int i3 = i;
//     int i4 = i+1;
//     int i5 = i+2;

//     // symmetric boundaries
//     if (i==0) {
//         i2 = 0;
//         i1 = 1;
//     } 
//     else if (i==1) {
//         i1 = 0;
//     }
//     if (i==num_points-1) {
//         i5 = num_points-2;
//         i4 = num_points-1;
//     } 
//     else if (i==num_points-2) {
//         i5 = num_points-1;
//     }

//     // sample field f_i(s_i)
//     double f1 = field[i1];
//     double f2 = field[i2];
//     double f3 = field[i3];
//     double f4 = field[i4];
//     double f5 = field[i5];

//     return (-f5 + 16.*f4 + -30.*f3 + 16.*f2 - f1)/(12. * ds * ds);
// }

// 1st derivative
double AHFinder::d2(const std::vector<double> &field, int i) {
    // integer sample points
    int i1 = i-1;
    int i2 = i;
    int i3 = i+1;

    // symmetric boundaries i1 = 1; 
    if (i==0) {
        i1 = 0;
    }
    else if (i==num_points-1) {
        i3 = num_points-1;
    }

    // sample field f_i(s_i)
    double f1 = field[i1];
    double f2 = field[i2];
    double f3 = field[i3];

    return (f3 - 2.*f2 + f1)/pow(ds,2);
}


// area infinitessimal
// inline this later?
double AHFinder::dA(const int i) {
    return 2. * Params::pi * pow(psi[i],4) 
            * sqrt((f[i]*f[i]) + d(f,i)*d(f,i))
            * f[i] * sin(sigma[i]) * ds;
}

// proper area
double AHFinder::area() {
    double A = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        A += dA(i);
    }
    // two because of reflective symmetry
    return 2.*A;
}

// flat (conformal) area
double AHFinder::area_flat() {
    double A = 0.;
    double dfds = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = sqrt((f[i]*f[i]) + d(f,i)*d(f,i));
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
    }
    // two because of reflective symmetry
    return 2.*A;
}

// irreducible horizon mass
double AHFinder::mass_irr() {
    return 0.25*sqrt(area()/(Params::pi));
}

// full horizon mass M^2 = M_irr^2 + J^2/(4*M_irr^2)
double AHFinder::mass() {
    double M_irr = mass_irr();
    double J = J_H();
    return sqrt(M_irr*M_irr + 0.25*J*J/M_irr/M_irr);
}

// J
double AHFinder::J_H() {
    double J = 0;
    for (size_t i = 0; i < num_points; ++i) {
        // calculate J_H, formula in notes/paper
        J +=  0.5 * 
        ( 
            f[i]*(W[i] + f[i]*dW_dR[i])*sin(sigma[i])
            - d(f,i)*(W[i]*cos(sigma[i]) + sin(sigma[i])*d(W,i))
        ) 
        * f[i] * sin(sigma[i]) * ds;
    }
    // factor of two for reflection symmetry
    return 2. * J;
}

// J/M, from 0 to M
double AHFinder::a_H() {
    return J_H()/mass();
}

// J/M^2, from 0 to 1
double AHFinder::chi_H() {
    return a_H()/mass();
}



double AHFinder::mass_MS() {
    double A = 0.;
    double B = 0.;
    double integrand = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // area 
        A += dA(i);
        // function integral B = int integrand dA
        integrand = (psi[i]-1.)*2.*f[i];
        B += dA(i) * integrand;
    }
    // two because of reflective symmetry
    return B/A;
}

double AHFinder::mass_SC() {
    double A = 0.;
    double B = 0.;
    double integrand = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // area 
        A += dA(i);
        // function integral B = int integrand dA
        integrand = (1.-pow(psi[i],-4))*f[i]/2.;
        B += dA(i) * integrand;
    }
    // two because of reflective symmetry
    return B/A;
}


double AHFinder::res() {
    double A = 0.;
    double B = 0.;
    double integrand = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // area 
        A += dA(i);
        B += dA(i) * dfdt[i];
    }
    // two because of reflective symmetry
    return B/A;
}


double AHFinder::psi_h() {
    double A = 0.;
    double B = 0.;
    double integrand = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // area 
        A += dA(i);
        B += dA(i) * psi[i];
    }
    // two because of reflective symmetry
    return B/A;
}

// average r
double AHFinder::r() {
    double A = 0.;
    double B = 0.;
    double integrand = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // area 
        A += dA(i);
        B += dA(i) * f[i];
    }
    // two because of reflective symmetry
    return B/A;
}

double AHFinder::eccentricity() {
    return f[0]/f[num_points-1];
}


void AHFinder::update(const Spacetime& spacetime) {
    // Placeholder for expansion/minimization logic
    // code here 
}

// to be called after we update the surface
// update external fields such as psi[f,s] 
void AHFinder::refresh(const Spacetime& spacetime) {
    double x = 0., y = 0.;
    double df = 0., dpsi = 0., dW=0.;
    double test_bh_mass = 1.0;
    for (size_t i = 0; i < num_points; ++i) {
        // numerical 
        df = ds * f[i];
        // update psi with new coords f,s
        x = f[i] * sin(sigma[i]);
        y = f[i] * cos(sigma[i]);
        psi[i] = spacetime.get_val_interp(spacetime.psi,x,y);
        W[i] = spacetime.get_val_interp(spacetime.W,x,y);
        // calculate d_psi/d_r
        x = (f[i] + df) * sin(sigma[i]);
        y = (f[i] + df) * cos(sigma[i]);
        dpsi = spacetime.get_val_interp(spacetime.psi,x,y);
        dW = spacetime.get_val_interp(spacetime.W,x,y);
        x = (f[i] - df) * sin(sigma[i]);
        y = (f[i] - df) * cos(sigma[i]);
        dpsi -= spacetime.get_val_interp(spacetime.psi,x,y);
        dW -= spacetime.get_val_interp(spacetime.W,x,y);
        // set values for spheical radial deriv
        dpsi_dR[i]=dpsi/(2.*df);
        dW_dR[i]=dW/(2.*df);

        // // algebraic schwarzschild
        // psi[i] = 1. + test_bh_mass/f[i]/2.;
        // dpsi_dR[i] = - test_bh_mass/f[i]/f[i]/2.;
    }

}

void AHFinder::relax() {
    // Poisson solver towards correct surface
    double term1 = 0.;
    double term2 = 0.;
    double term3 = 0.;
    double dfds = 0.;
    double R = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // dfdt[i] = - 0.2 * f[i];
        dfds = d(f,i);
        R = f[i];
        //
        term1 = (dfds + pow(dfds,3)/R/R)*( cos(sigma[i])/sin(sigma[i]) + 4.*d(psi,i)/psi[i]);
        //
        term2 = - dfds * dfds * (4.*dpsi_dR[i]/psi[i] + 3./R);
        //
        term3 = -2.*R - 4.*R*R*dpsi_dR[i]/psi[i];
        //
        dfdt[i] = d2(f,i) + term1 + term2 + term3;
    }
    for (size_t i = 0; i < num_points; ++i) {
        f[i] = f[i] + dfdt[i] * dt;
    }
}

void AHFinder::save(const std::string &filename) {
    std::ofstream outFile("data/" + filename + ".dat");
    if (!outFile.is_open()) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }

    // sig figs 
    outFile << std::scientific << std::setprecision(Params::save_precision);

    // for (size_t i = 0; i < num_points; ++i) {
    //     outFile << f[i];
    //     if (i != num_points - 1) {
    //         outFile << ","; // comma separator
    //     }
    // }

    for (size_t i = 0; i < num_points; ++i) {
        outFile << f[i] << "," << sigma[i] <<"\n";
    }

    outFile.close();
}

void AHFinder::hello() const {
    std::cout << "Hello from AHFinder!\n"
              << " - Number of horizon points = "
              << num_points << std::endl;
}

