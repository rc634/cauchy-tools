#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>

Spacetime::Spacetime(int nx_, int ny_)
    : nx(nx_), ny(ny_)
{
    psi.resize(nx * ny);
    W.resize(nx * ny);
}

void Spacetime::hello() const {
    std::cout << "Hello from Spacetime! "
              << "Grid size = " << nx << " x " << ny << std::endl;

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

inline int Spacetime::index(int i, int j) const {
// #ifndef DEBUG
//     if (i < 0 || i >= nx || j < 0 || j >= ny) {
//         throw std::out_of_range("Spacetime index out of bounds");
//     }
// #endif
    return i + nx * j;
}




