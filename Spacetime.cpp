#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>

Spacetime::Spacetime() {
    
    nx = Params::NX; 
    ny = Params::NY;
    xl = Params::XL;
    xu = Params::XU;
    yl = Params::YL;
    yu = Params::YU;
    dx = (xu-xl)/nx;
    dy = (yu-yl)/ny;

    psi.resize(nx * ny);
    W.resize(nx * ny);
}

void Spacetime::hello() const {
    std::cout << "Hello from Spacetime! \n"
              << " - Grid size = " << nx << " x " << ny << "\n"
              << " - Boundaries at [XL:" << xl << ",XU:" << xu 
              << ",YL:" << yl << ",YU:" << yu
              << ",DX:" << dx << ",DY:" << dy << "]\n";

    // debugging / writing
    // for (size_t i = 0; i < nx; i++)
    // {
    //     std::cout << "Spacetime : " << get_val(psi,i,i) << "\n";
    // }
    
}

const double Spacetime::get_val(const std::vector<double>& field, int i, int j) const {
// #ifndef DEBUG
    // if (field.size() != static_cast<size_t>(nx * ny)) {
    //     throw std::runtime_error("Field size does not match spacetime grid");
    // }
// #endif
    return field[index(i, j)];
}


// ---- Bilinear interpolation in physical coordinates ----
const double Spacetime::get_val_interp(const std::vector<double>& field, double x, double y) const {
    // Convert physical coordinates to fractional grid indices

    // CAREFUL NEED TO IMPLEMENT CELL CENTERED STUFF!
    double i = (x - 0.5*dx) / dx;
    double j = (y - 0.5*dy) / dy;

    // Floor indices
    int i0 = static_cast<int>(std::floor(i));
    int j0 = static_cast<int>(std::floor(j));
    int i1 = i0 + 1;
    int j1 = j0 + 1;

    // Clamp indices to grid
    i0 = std::clamp(i0, 0, nx - 1);
    i1 = std::clamp(i1, 0, nx - 1);
    j0 = std::clamp(j0, 0, ny - 1);
    j1 = std::clamp(j1, 0, ny - 1);

    // Fractional distances
    double di = i - i0;
    double dj = j - j0;

    // Four surrounding points
    double f00 = field[index(i0, j0)];
    double f10 = field[index(i1, j0)];
    double f01 = field[index(i0, j1)];
    double f11 = field[index(i1, j1)];

    // Bilinear interpolation
    double val = f00 * (1 - di) * (1 - dj)
               + f10 * di * (1 - dj)
               + f01 * (1 - di) * dj
               + f11 * di * dj;

    return val;
}

// ---- Bilinear interpolation in physical coordinates ----
const double Spacetime::get_ddx_interp(const std::vector<double>& field, double x, double y) const {
    // Convert physical coordinates to fractional grid indices

    // CAREFUL NEED TO IMPLEMENT CELL CENTERED STUFF!
    double i = (x - 0.5*dx) / dx;
    double j = (y - 0.5*dy) / dy;

    // Floor indices
    int i0 = static_cast<int>(std::floor(i));
    int j0 = static_cast<int>(std::floor(j));
    int i1 = i0 + 1;
    int j1 = j0 + 1;

    // Clamp indices to grid
    i0 = std::clamp(i0, 0, nx - 1);
    i1 = std::clamp(i1, 0, nx - 1);
    j0 = std::clamp(j0, 0, ny - 1);
    j1 = std::clamp(j1, 0, ny - 1);

    // Fractional distances
    double di = i - i0;
    double dj = j - j0;

    // Four surrounding points
    double f00 = field[index(i0, j0)];
    double f10 = field[index(i1, j0)];
    double f01 = field[index(i0, j1)];
    double f11 = field[index(i1, j1)];

    // Bilinear interpolation
    double val = (f10 - f00) * (1. - dj) +
                 (f11 - f01) * dj;

    return val / dx; // dx needed as di dj are fractional coordinates 
}

// ---- Bilinear interpolation in physical coordinates ----
const double Spacetime::get_ddy_interp(const std::vector<double>& field, double x, double y) const {
    // Convert physical coordinates to fractional grid indices

    // CAREFUL NEED TO IMPLEMENT CELL CENTERED STUFF!
    double i = (x - 0.5*dx) / dx;
    double j = (y - 0.5*dy) / dy;

    // Floor indices
    int i0 = static_cast<int>(std::floor(i));
    int j0 = static_cast<int>(std::floor(j));
    int i1 = i0 + 1;
    int j1 = j0 + 1;

    // Clamp indices to grid
    i0 = std::clamp(i0, 0, nx - 1);
    i1 = std::clamp(i1, 0, nx - 1);
    j0 = std::clamp(j0, 0, ny - 1);
    j1 = std::clamp(j1, 0, ny - 1);

    // Fractional distances
    double di = i - i0;
    double dj = j - j0;

    // Four surrounding points
    double f00 = field[index(i0, j0)];
    double f10 = field[index(i1, j0)];
    double f01 = field[index(i0, j1)];
    double f11 = field[index(i1, j1)];

    // Bilinear interpolation
    double val = (f11 - f10) * di +
                 (f01 - f00) * (1. - di);

    return val / dy; // dy needed as di dj are fractional coordinates 
}





inline int Spacetime::index(int i, int j) const {
// #ifndef DEBUG
//     if (i < 0 || i >= nx || j < 0 || j >= ny) {
//         throw std::out_of_range("Spacetime index out of bounds");
//     }
// #endif

    // CAREFUL, CHECK THIS IS CORRECT
    return i + nx * j;
}

// set nakamura data 
void Spacetime::set_data_nakamura_prolate() {

    // coords 
    double x=0., y=0.;

    // zero W 
    for (size_t k = 0; k < nx*ny-1; k++)
    {
        W[k] = 0.;
    }

    // nakamura const
    double M_N = 1.; // assumed implicitly for now
    double M = 2. * M_N;
    double e = 0.99;
    double c = M * 0.4;
    double a = c*sqrt(1-e*e);
    double k = 1.5 * pow(c*e,-3);
    double alpha = 6. / (5.*c*e) * log((1+e)/(1-e));
    double M_0 = 2. + alpha;
    double beta_internal = std::asinh(c*e/a); // internal beta
    double pho = -1.5/(c*e);
    // nakamura vars
    double phiN = 0.;
    double Kr = 0;
    double Kz = 0.;
    double beta = 0.;
    // test
    bool internal = false;

    // loop i and j and set psi
    for (size_t i = 0; i < nx-1; i++)
    {
        x = (i+0.5)*dx;
        for (size_t j = 0; j < ny-1; j++) {
            y = (j+0.5)*dy;

            // are we inside the nakamura ellipsoid?
            if (x*x/a/a + y*y/c/c < 1.) {
                beta = beta_internal;
                internal = true;
            }
            else {
                beta = beta_external_prolate(x,y,c,e);
            }

            // brrrrrrt!
            Kr = k * x * (beta - sinh(beta)*cosh(beta)); 
            Kz = 2. * k * y * (tanh(beta) - beta); 
            phiN = pho*beta - 0.5*(x*Kr+y*Kz);

            // output psi!
            psi[index(i,j)] = 1. - phiN;
        }
    }
}

double Spacetime::beta_external_prolate(double x, double y, double c, double e) {
    double A = x*x;
    double B = x*x + y*y - c*c*e*e;
    double C = -c*c*e*e;
    // quadratic formula
    double ss = (- B + sqrt(B*B - 4.*A*C))/(2.*A);
    return std::asinh(sqrt(ss));
}

// set nakamura data 
void Spacetime::set_data_nakamura_oblate() {

    // coords 
    double x=0., y=0.;

    // zero W 
    for (size_t k = 0; k < nx*ny-1; k++)
    {
        W[k] = 0.;
    }

    // nakamura const
    double M_N = 1.; // assumed implicitly for now
    double M = 2. * M_N;
    double e = 0.99;
    double a = M * 0.54;
    double c = a*sqrt(1-e*e);
    double k = pow(a*e,-3);
    double beta_internal = std::asin(e); // internal beta
    double pho = -1.5/(a*e);
    // nakamura vars
    double phiN = 0.;
    double Kr = 0;
    double Kz = 0.;
    double beta = 0.;
    // test
    bool internal = false;

    // loop i and j and set psi
    for (size_t i = 0; i < nx-1; i++)
    {
        x = (i+0.5)*dx;
        for (size_t j = 0; j < ny-1; j++) {
            y = (j+0.5)*dy;

            // are we inside the nakamura ellipsoid?
            if (x*x/a/a + y*y/c/c < 1.) {
                beta = beta_internal;
                internal = true;
            }
            else {
                beta = beta_external_oblate(x,y,a,e);
            }

            // brrrrrrt!
            Kr = - 1.5 * k * x * (beta - sin(beta)*cos(beta)); 
            Kz = - 3.0 * k * y * (tan(beta) - beta); 
            phiN = pho*beta - 0.5*(x*Kr+y*Kz);

            // output psi!
            psi[index(i,j)] = 1. - phiN;
        }
    }
}

// this one is subtle, in the limit x and y are large the solution is 
// 0 or (r^2 + z^2) / r^2. we dont want bigger than one so we take negative root
double Spacetime::beta_external_oblate(double x, double y, double a, double e) {
    double A = x*x;
    double B = - x*x - y*y - a*a*e*e;
    double C = a*a*e*e;
    // quadratic formula MINUS ROOT!
    double ss = (- B - sqrt(B*B - 4.*A*C))/(2.*A);
    return std::asin(sqrt(ss));
}

// set BH data 
void Spacetime::set_data_BH() {

    // coords 
    double x=0., y=0.;

    // zero W 
    for (size_t k = 0; k < nx*ny-1; k++)
    {
        W[k] = 0.;
    }

    // BH mass const
    double M = 1.;

    // loop i and j and set psi
    for (size_t i = 0; i < nx-1; i++)
    {
        x = (i+0.5)*dx;

        for (size_t j = 0; j < ny-1; j++) {

            y = (j+0.5)*dy;

            psi[index(i,j)] = 1. + 0.5*M/sqrt(x*x+y*y);
        }
    }
}




