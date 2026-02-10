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
              << ",YL:" << yl << ",YU:" << yu << "]\n";

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
    double i = x / dx;
    double j = y / dy;

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




inline int Spacetime::index(int i, int j) const {
// #ifndef DEBUG
//     if (i < 0 || i >= nx || j < 0 || j >= ny) {
//         throw std::out_of_range("Spacetime index out of bounds");
//     }
// #endif
    return i + nx * j;
}




