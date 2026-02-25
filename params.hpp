#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

namespace Params {

    constexpr int NX = 10000;
    constexpr int NY = 10000;

    // coords of x/y axes, may not be cartesian
    constexpr double XL = 0.;
    constexpr double XU = 12.;
    constexpr double YL = 0.;
    constexpr double YU = 12.;

    // courant friedrich lewis thingy 
    constexpr double CFL = 0.25; 

    // sample points
    constexpr int AH_NPOINTS = 128; 
    constexpr int EX_NPOINTS = 2048;
    constexpr int SH_NPOINTS = 40000;

    // initial radii
    constexpr double RH = 5.;
    constexpr double RX = 5.;

    // consts 
    constexpr double pi = M_PI;

    // surface resolution
    // constexpr int SURFACE_N_POINTS = 16;

    // save precision
    constexpr int save_precision = 8;

}

#endif

