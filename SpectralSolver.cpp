#include "SpectralSolver.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include "mymaths.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

SpectralSolver::SpectralSolver(int npoints, int Nspec) {
    sigma.resize(npoints);
    f.resize(npoints);
    dfds.resize(npoints);
    psi.resize(npoints);
    W.resize(npoints);
    dpsi_dR.resize(npoints);
    dW_dR.resize(npoints);
    dfdt.resize(npoints);
    coeff.resize(Nspec);
    J.resize(Nspec);
    num_points = npoints;
    N = Nspec;

    // double dx = (Params::XU-Params::XL)/Params::NX;
    // double dy = (Params::YU-Params::YL)/Params::NY;
    //dt = Params::CFL * std::min(dx,dy);
    ds = 0.5 * Params::pi / num_points;
    dt = Params::CFL * ds * ds;
}

void SpectralSolver::initialize(const Spacetime& spacetime, const double f0) {
    double x = 0., y = 0.;

    // initial Legendre coefficients
    for (size_t j = 0; j < N; j++) coeff[j] = 0.;
    for (size_t j = 0; j < N; j++) J[j] = 0.;
    for (size_t j = 0; j < N; j++)  if (j%2==0)  {
        coeff[j] = 5./(1.+j*j*pow(-1,1+j/2.));
    }

    for (size_t i = 0; i < num_points; ++i) {

        // cell centred
        
        
        sigma[i] = (i + 1./2.) * ds;
        f[i] = mymaths::LegendreSeries(coeff,cos(sigma[i]));
        dfds[i] = mymaths::LegendreSeriesDiff1(coeff,cos(sigma[i]));
        dfdt[i] = 0.;
        psi[i] = 0.;
        dpsi_dR[i] = 0.;
        W[i] = 0.;
        dW_dR[i] = 0.;

    }

    std::cout << " - setting f[i] = " << f0 << std::endl;
}

void SpectralSolver::gradient_descent(const Spacetime& spacetime) {
    //
    double RES0 = 0.;
    double delta = 0.001;
    double delta_min = 0.00001;
    double modJ = 0.;
    double alpha = 0.01;
    double alpha_min = 0.001;
    int vk = 500; // verbosity in k
    double halflife = 5000.; // alpha gets smaller
    for (size_t k = 0; k < 15001; k++) {


        refresh(spacetime);
        residual();
        RES0 = res();
        alpha *= (1. - 1./halflife);
        if (alpha<alpha_min) alpha = alpha_min;
        if (delta>delta_min) delta *= (1.-1./halflife);

        if(k%vk==0) std::cout << " * - * descent step = " << k << " * initial res : " << RES0 << "\n";

        for (size_t j = 0; j < N; j++)  if (j%2==0)  {

            coeff[j] += delta;

            for (size_t i = 0; i < num_points; ++i) {
                f[i] = mymaths::LegendreSeries(coeff,cos(sigma[i]));
                dfds[i] = mymaths::LegendreSeriesDiff1(coeff,cos(sigma[i]));
            }

            refresh(spacetime);
            residual();

            J[j] = (res()-RES0)/delta;
            modJ += J[j]*J[j];
            // if(k%vk==0) std::cout << " * j = " << j << " : cj " << coeff[j] << ", J : " << J[j] << "\n";

            coeff[j] -= delta; // undo teh change for next deriv
        }

        modJ = sqrt(modJ);

        for (size_t j = 0; j < N; j++)  if (j%2==0)  {
            coeff[j] -= alpha*J[j]/modJ;
        }

        for (size_t i = 0; i < num_points; ++i) {
                f[i] = mymaths::LegendreSeries(coeff,cos(sigma[i]));
                dfds[i] = mymaths::LegendreSeriesDiff1(coeff,cos(sigma[i]));
            }

        refresh(spacetime);
        residual();
        modJ = 0.;
    }
}

// // 1st derivative // 4th order
// double SpectralSolver::d(const std::vector<double> &field, int i) {
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
double SpectralSolver::d(const std::vector<double> &field, int i) {
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
// double SpectralSolver::d2(const std::vector<double> &field, int i) {
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
double SpectralSolver::d2(const std::vector<double> &field, int i) {
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
double SpectralSolver::dA(const int i) {
    return 2. * Params::pi * pow(psi[i],4) 
            * sqrt((f[i]*f[i]) + d(f,i)*d(f,i))
            * f[i] * sin(sigma[i]) * ds;
}

// proper area
double SpectralSolver::area() {
    double A = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        A += dA(i);
    }
    // two because of reflective symmetry
    return 2.*A;
}

// flat (conformal) area
double SpectralSolver::area_flat() {
    double A = 0.;
    double geom = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        geom = sqrt((f[i]*f[i]) + dfds[i]*dfds[i]);
        A += 2. * Params::pi * f[i] * sin(sigma[i]) * geom * ds;
    }
    // two because of reflective symmetry
    return 2.*A;
}

// irreducible horizon mass
double SpectralSolver::mass_irr() {
    return 0.25*sqrt(area()/(Params::pi));
}

// full horizon mass M^2 = M_irr^2 + J^2/(4*M_irr^2)
double SpectralSolver::mass() {
    double M_irr = mass_irr();
    double J = J_H();
    return sqrt(M_irr*M_irr + 0.25*J*J/M_irr/M_irr);
}

// J
double SpectralSolver::J_H() {
    double J = 0.;
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
double SpectralSolver::a_H() {
    return J_H()/mass();
}

// J/M^2, from 0 to 1
double SpectralSolver::chi_H() {
    return a_H()/mass();
}



double SpectralSolver::mass_MS() {
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

double SpectralSolver::mass_SC() {
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


double SpectralSolver::res() {
    double A = 0.;
    double B = 0.;
    double C = 0.;
    double D = 0.;
    double integrand = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // area 
        A += dA(i);
        B += dA(i) * dfdt[i];
        C += dfdt[i]*dfdt[i] * sin(sigma[i]);
        D += dfdt[i]*dfdt[i];
    }
    // two because of reflective symmetry
    return C/num_points;
}


double SpectralSolver::psi_h() {
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
double SpectralSolver::r() {
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

double SpectralSolver::eccentricity() {
    double f0 = f[0];
    double fE = f[num_points-1];
    if (fE==f0) return 0;
    else if (fE>f0) {
        return sqrt(1-pow(f0/fE,2));
    }
    else {
        return sqrt(1-pow(fE/f0,2));
    }
}

// to be called after we update the surface
// update external fields such as psi[f,s] 
void SpectralSolver::refresh(const Spacetime& spacetime) {
    double x = 0., y = 0.;
    double df = 0., dpsi = 0., dW=0.;
    double test_bh_mass = 1.0;
    double TH = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // for numerical spherical radius R deriv
        df = ds * f[i];
        TH = sigma[i];
        // update psi with new coords f,s
        x = f[i] * sin(TH);

        y = f[i] * cos(TH);

        psi[i] = spacetime.get_val_interp(spacetime.psi,x,y);

        W[i] = spacetime.get_val_interp(spacetime.W,x,y);

        dpsi_dR[i] = cos(sigma[i]) * spacetime.get_ddy_interp(spacetime.psi,x,y) 
                    + sin(sigma[i]) * spacetime.get_ddx_interp(spacetime.psi,x,y);
                    
        dW_dR[i] = cos(sigma[i]) * spacetime.get_ddy_interp(spacetime.W,x,y) 
                    + sin(sigma[i]) * spacetime.get_ddx_interp(spacetime.W,x,y);
  
        // //  old method 
        // // calculate d_psi/d_r
        // x = (f[i] + df) * sin(sigma[i]);
        // y = (f[i] + df) * cos(sigma[i]);
        // dpsi = spacetime.get_val_interp(spacetime.psi,x,y);
        // dW = spacetime.get_val_interp(spacetime.W,x,y);
        // x = (f[i] - df) * sin(sigma[i]);
        // y = (f[i] - df) * cos(sigma[i]);
        // dpsi -= spacetime.get_val_interp(spacetime.psi,x,y);
        // dW -= spacetime.get_val_interp(spacetime.W,x,y);


        // // set values for spheical radial deriv
        // dpsi_dR[i]=dpsi/(2.*df);
        // dW_dR[i]=dW/(2.*df);
    }
}

// to be called after we update the surface
// update external fields such as psi[f,s] 
void SpectralSolver::refresh_pureBH() {
    double test_bh_mass = 1.0;
    for (size_t i = 0; i < num_points; ++i) {

        // algebraic schwarzschild
        psi[i] = 1. + test_bh_mass/f[i]/2.;
        dpsi_dR[i] = - test_bh_mass/f[i]/f[i]/2.;
        W[i] = 0.;
        dW_dR[i] = 0.;
    }

}


void SpectralSolver::residual() {
    // Poisson solver towards correct surface
    double term1 = 0.;
    double term2 = 0.;
    double term3 = 0.;
    double _dfds_ = 0.;
    double R = 0.;
    for (size_t i = 0; i < num_points; ++i) {
        // dfdt[i] = - 0.2 * f[i];
        _dfds_ = dfds[i];//d(f,i);
        R = f[i];
        //
        term1 = (_dfds_ + pow(_dfds_,3)/R/R)*( cos(sigma[i])/sin(sigma[i]) + 4.*d(psi,i)/psi[i]);
        //
        term2 = - _dfds_ * _dfds_ * (4.*dpsi_dR[i]/psi[i] + 3./R);
        //
        term3 = -2.*R - 4.*R*R*dpsi_dR[i]/psi[i];
        //
        dfdt[i] = abs(d2(f,i) + term1 + term2 + term3);
    }
    // removed moving the surface
    // for (size_t i = 0; i < num_points; ++i) {
    //     f[i] = f[i] + dfdt[i] * dt;
    // }
}



// stuff 

void SpectralSolver::save(const std::string &filename) {
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

void SpectralSolver::hello() const {
    std::cout << "Hello from SpectralSolver!\n"
              << " - Number of horizon points = "
              << num_points << std::endl;
}

