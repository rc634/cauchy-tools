#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

namespace Params {

    constexpr int NX = 1024;
    constexpr int NY = 1024;

    // coords of x/y axes, may not be cartesian
    constexpr double XL = 0.;
    constexpr double XU = 40.;
    constexpr double YL = 0.;
    constexpr double YU = 40.;

    // courant friedrich lewis thingy 
    constexpr double CFL = 0.25; 
    
    // legendre coefficient order max
    constexpr int PL_LMAX = 16; 

    // sample points
    constexpr int AH_NPOINTS = 128; // relax
    constexpr int SP_NPOINTS = 512;  // spectral
    constexpr int EX_NPOINTS = 2048; // extraction
    constexpr int SH_NPOINTS = 16000; //shoot

    // initial radii
    constexpr double RH = 5.;
    constexpr double RX = 5.;
    constexpr double RS = 5.;

    // consts 
    constexpr double pi = M_PI;

    // surface resolution
    // constexpr int SURFACE_N_POINTS = 16;

    // save precision
    constexpr int cout_precision = 8;
    constexpr int save_precision = 12;

}

#endif

