#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

namespace Params {

    constexpr int NX = 128;
    constexpr int NY = 128;

    // coords of x/y axes, may not be cartesian
    constexpr double XL = 0.;
    constexpr double XU = 2.;
    constexpr double YL = 0.;
    constexpr double YU = 2.;

    // courant friedrich lewis thingy 
    constexpr double CFL = 0.25; 

    // apparent horizon sample points
    constexpr int AH_NPOINTS = 400; 

    constexpr double pi = M_PI;

}

#endif

