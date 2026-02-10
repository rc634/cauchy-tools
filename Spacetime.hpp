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
    double dx, dy;

    double xl, xu, yl, yu; 

    Spacetime();

    // on-gridpoint accessor
    const double get_val(const std::vector<double>& field, int i, int j) const;
    // inbetween gridpoint interpolator
    const double get_val_interp(const std::vector<double>& field, double x, double y) const;


    void hello() const;

private:

    inline int index(int i, int j) const;
};

#endif

