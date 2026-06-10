#ifndef SPACETIME_HPP
#define SPACETIME_HPP

#include <vector>
#include <string>

class Spacetime {
public:
    // Raw 2D fields
    // can be filled analytically (functinos in this file)
    // can also be read from file using the data loader (from main)
    std::vector<double> psi;
    std::vector<double> W;
    std::vector<double> src;
    std::vector<double> srcW;
    std::vector<double> src_eff;
    std::vector<double> srcW_eff;

    int nx, ny;
    double dx, dy;

    double xl, xu, yl, yu; 

    Spacetime();

    // on-gridpoint accessor
    const double get_val(const std::vector<double>& field, int i, int j) const;

    // inbetween gridpoint interpolator
    const double get_val_interp(const std::vector<double>& field, double x, double y) const;

    // bicubic interpolator
    const double get_val_bicubic(const std::vector<double>& field, double x, double y) const;
    const double get_ddx_bicubic(const std::vector<double>& field, double x, double y) const;
    const double get_ddy_bicubic(const std::vector<double>& field, double x, double y) const;

    // get derivatives
    const double get_ddx_interp(const std::vector<double>& field, double x, double y) const;
    const double get_ddy_interp(const std::vector<double>& field, double x, double y) const;

    // set initial data like Nakamuras solutions 
    void set_data_nakamura_prolate();
    double beta_external_prolate(double x, double y, double c, double e);
    void set_data_nakamura_oblate();
    double beta_external_oblate(double x, double y, double a, double e);

    // BH or perturbed BH to test code convergence
    void set_data_BH();

    void hello() const;

private:

    mutable double f_ij[4][4] = {};   // 4x4 grid values of current patch
    mutable double a_mn[4][4] = {};   // bicubic coefficients for current patch

    void fill_a_mn(const std::vector<double>& field, int i, int j) const;
    void fill_f_ij(const std::vector<double>& field, int i, int j) const;

    inline int index(int i, int j) const;
};

#endif

