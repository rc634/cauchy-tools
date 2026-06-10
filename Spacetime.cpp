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
    src.resize(nx * ny, 0.0);
    srcW.resize(nx * ny, 0.0);
    src_eff.resize(nx * ny, 0.0);
    srcW_eff.resize(nx * ny, 0.0);
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



// BICUBIC SPLINE STUFF


// ---- Bicubic interpolation in physical coordinates ----
const double Spacetime::get_val_bicubic(const std::vector<double>& field, double x, double y) const {
    // cell-centred grid: x_i = (i + 0.5)*dx  →  i = floor(x/dx - 0.5)
    int i = std::clamp(static_cast<int>(std::floor(x / dx - 0.5)), 0, nx - 2);
    int j = std::clamp(static_cast<int>(std::floor(y / dy - 0.5)), 0, ny - 2);

    // local coordinates in [0, 1]: u=0 at grid point i, u=1 at i+1
    double u = x / dx - 0.5 - i;
    double v = y / dy - 0.5 - j;

    fill_a_mn(field, i, j);

    double val = 0.;
    double um = 1.;
    for (int m = 0; m < 4; ++m) {
        double vn = 1.;
        for (int n = 0; n < 4; ++n) {
            val += a_mn[m][n] * um * vn;
            vn  *= v;
        }
        um *= u;
    }
    return val;
}

const double Spacetime::get_ddx_bicubic(const std::vector<double>& field, double x, double y) const {
    int i = std::clamp(static_cast<int>(std::floor(x / dx - 0.5)), 0, nx - 2);
    int j = std::clamp(static_cast<int>(std::floor(y / dy - 0.5)), 0, ny - 2);

    double u = x / dx - 0.5 - i;
    double v = y / dy - 0.5 - j;

    fill_a_mn(field, i, j);

    double val = 0.;
    double um = 1.;                  // u^(m-1), starting at m=1
    for (int m = 1; m < 4; ++m) {
        double vn = 1.;
        for (int n = 0; n < 4; ++n) {
            val += m * a_mn[m][n] * um * vn;
            vn  *= v;
        }
        um *= u;
    }
    return val / dx;
}

const double Spacetime::get_ddy_bicubic(const std::vector<double>& field, double x, double y) const {
    int i = std::clamp(static_cast<int>(std::floor(x / dx - 0.5)), 0, nx - 2);
    int j = std::clamp(static_cast<int>(std::floor(y / dy - 0.5)), 0, ny - 2);

    double u = x / dx - 0.5 - i;
    double v = y / dy - 0.5 - j;

    fill_a_mn(field, i, j);

    double val = 0.;
    double um = 1.;
    for (int m = 0; m < 4; ++m) {
        double vn = 1.;                  // v^(n-1), starting at n=1
        for (int n = 1; n < 4; ++n) {
            val += n * a_mn[m][n] * um * vn;
            vn  *= v;
        }
        um *= u;
    }
    return val / dy;
}

void Spacetime::fill_f_ij(const std::vector<double>& field, int i, int j) const {
    //
    for (int di = 0; di < 4; ++di) {

        int ii = i - 1 + di;

        if (ii < 0)       ii = -ii;       // axial symmetry: reflect across x=0
        else if (ii >= nx) ii = nx - 1;  // clamp to last valid index

        for (int dj = 0; dj < 4; ++dj) {
            int jj = j - 1 + dj;
            if (jj < 0)        jj = -jj;       // equatorial symmetry: reflect across y=0
            else if (jj >= ny) jj = ny - 1;    // clamp to last valid index
            f_ij[di][dj] = field[index(ii, jj)];
        }
        
    }
}

void Spacetime::fill_a_mn(const std::vector<double>& field, int i, int j) const {
    fill_f_ij(field, i, j);

    const double f00=f_ij[0][0], f01=f_ij[0][1], f02=f_ij[0][2], f03=f_ij[0][3];
    const double f10=f_ij[1][0], f11=f_ij[1][1], f12=f_ij[1][2], f13=f_ij[1][3];
    const double f20=f_ij[2][0], f21=f_ij[2][1], f22=f_ij[2][2], f23=f_ij[2][3];
    const double f30=f_ij[3][0], f31=f_ij[3][1], f32=f_ij[3][2], f33=f_ij[3][3];

    a_mn[0][0] = f11;
    a_mn[0][1] = (1./6.) * (-2.*f10 - 3.*f11 + 6.*f12 - f13);
    a_mn[0][2] = 0.5 * (f10 - 2.*f11 + f12);
    a_mn[0][3] = (1./6.) * (-f10 + 3.*f11 - 3.*f12 + f13);

    a_mn[1][0] = (1./6.) * (-2.*f01 - 3.*f11 + 6.*f21 - f31);
    a_mn[1][1] = (1./36.) * ( 4.*f00 +  6.*f01 - 12.*f02 +  2.*f03
                             + 6.*f10 +  9.*f11 - 18.*f12 +  3.*f13
                             -12.*f20 - 18.*f21 + 36.*f22 -  6.*f23
                             + 2.*f30 +  3.*f31 -  6.*f32 +     f33);
    a_mn[1][2] = (1./12.) * (-2.*f00 +  4.*f01 -  2.*f02
                             - 3.*f10 +  6.*f11 -  3.*f12
                             + 6.*f20 - 12.*f21 +  6.*f22
                             -    f30 +  2.*f31 -     f32);
    a_mn[1][3] = (1./36.) * ( 2.*f00 -  6.*f01 +  6.*f02 -  2.*f03
                             + 3.*f10 -  9.*f11 +  9.*f12 -  3.*f13
                             - 6.*f20 + 18.*f21 - 18.*f22 +  6.*f23
                             +    f30 -  3.*f31 +  3.*f32 -     f33);

    a_mn[2][0] = 0.5 * (f01 - 2.*f11 + f21);
    a_mn[2][1] = (1./12.) * (-2.*f00 -  3.*f01 +  6.*f02 -     f03
                             + 4.*f10 +  6.*f11 - 12.*f12 +  2.*f13
                             - 2.*f20 -  3.*f21 +  6.*f22 -     f23);
    a_mn[2][2] = 0.25 * (   f00 - 2.*f01 +    f02
                         - 2.*f10 + 4.*f11 - 2.*f12
                         +    f20 - 2.*f21 +    f22);
    a_mn[2][3] = (1./12.) * (-   f00 +  3.*f01 -  3.*f02 +     f03
                             + 2.*f10 -  6.*f11 +  6.*f12 -  2.*f13
                             -    f20 +  3.*f21 -  3.*f22 +     f23);

    a_mn[3][0] = (1./6.) * (-f01 + 3.*f11 - 3.*f21 + f31);
    a_mn[3][1] = (1./36.) * ( 2.*f00 +  3.*f01 -  6.*f02 +     f03
                             - 6.*f10 -  9.*f11 + 18.*f12 -  3.*f13
                             + 6.*f20 +  9.*f21 - 18.*f22 +  3.*f23
                             - 2.*f30 -  3.*f31 +  6.*f32 -     f33);
    a_mn[3][2] = (1./12.) * (-   f00 +  2.*f01 -     f02
                             + 3.*f10 -  6.*f11 +  3.*f12
                             - 3.*f20 +  6.*f21 -  3.*f22
                             +    f30 -  2.*f31 +     f32);
    a_mn[3][3] = (1./36.) * (    f00 -  3.*f01 +  3.*f02 -     f03
                             - 3.*f10 +  9.*f11 -  9.*f12 +  3.*f13
                             + 3.*f20 -  9.*f21 +  9.*f22 -  3.*f23
                             -    f30 +  3.*f31 -  3.*f32 +     f33);
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







// NAKAMURA SOLS




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
    double x=0., y=0., r=0., cos_theta=0.;

    // zero W 
    for (size_t k = 0; k < nx*ny-1; k++)
    {
        W[k] = 0.;
    }

    // BH mass const
    double M = 1.;

    // perturbative ellipse
    double eps = 0.003;
    double PL = 0.; // legendre polynomial

    // loop i and j and set psi
    for (size_t i = 0; i <= nx-1; i++)
    {
        x = (i+0.5)*dx;

        for (size_t j = 0; j <= ny-1; j++) {

            y = (j+0.5)*dy;

            r = sqrt(x*x+y*y);

            cos_theta = y/r;

            PL = 0.5*(3.*cos_theta*cos_theta - 1.); // L=2

            psi[index(i,j)] = 1. + 0.5*M/r + eps * pow(0.5*M/r,3) * PL;
        }
    }
}




