#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

namespace Params {

    constexpr int NX = 256;
    constexpr int NY = 256;

    // coords of x/y axes, may not be cartesian
    constexpr double XL = 0.;
    constexpr double XU = 50.;
    constexpr double YL = 0.;
    constexpr double YU = 50.;

    // courant friedrich lewis thingy 
    constexpr double CFL = 0.25; 

    // apparent horizon sample points
    constexpr int AH_NPOINTS = 64; 
    constexpr double RH = 5.;
    constexpr double RX = 40.;

    constexpr double pi = M_PI;

    // surface resolution
    // constexpr int SURFACE_N_POINTS = 16;

    // save precision
    constexpr int save_precision = 8;

}

#endif

