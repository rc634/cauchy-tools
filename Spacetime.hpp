#ifndef SPACETIME_HPP
#define SPACETIME_HPP

#include <vector>
#include <string>

class Spacetime {
public:
    // Raw 2D fields
    std::vector<double> psi;
    std::vector<double> W;

    int nx, ny;

    Spacetime(int nx_, int ny_);

    // on-gridpoint accessor
    const double get_val(const std::vector<double>& field, int i, int j) const;

    void hello() const;

private:

    inline int index(int i, int j) const;
};

#endif

