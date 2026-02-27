#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

namespace Params {

    constexpr int NX = 1024;
    constexpr int NY = 1024;

    // coords of x/y axes, may not be cartesian
    constexpr double XL = 0.;
    constexpr double XU = 20.;
    constexpr double YL = 0.;
    constexpr double YU = 20.;

    // courant friedrich lewis thingy 
    constexpr double CFL = 0.25; 
    
    // legendre coefficient order max
    constexpr int PL_LMAX = 24; 

    // sample points
    constexpr int AH_NPOINTS = 256; 
    constexpr int SP_NPOINTS = 512; 
    constexpr int EX_NPOINTS = 2048;
    constexpr int SH_NPOINTS = 40000;

    // initial radii
    constexpr double RH = 5.;
    constexpr double RX = 5.;
    constexpr double RS = 5.;

    // consts 
    constexpr double pi = M_PI;

    // surface resolution
    // constexpr int SURFACE_N_POINTS = 16;

    // save precision
    constexpr int save_precision = 8;

}

#endif

