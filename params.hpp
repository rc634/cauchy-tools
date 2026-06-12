#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <cmath>

namespace Params {

    constexpr int NX = 4096; //  8192;
    constexpr int NY = 4096; //  8192;

    // coords of x/y axes, may not be cartesian
    constexpr double XL = 0.;
    constexpr double XU = 40. ; // 80. ;
    constexpr double YL = 0.;
    constexpr double YU = 40. ; // 80. ;

    // courant friedrich lewis thingy 
    constexpr double CFL = 0.25; 
    
    // legendre coefficient order max
    constexpr int PL_LMAX = 16; 

    // sample points
    constexpr int AH_NPOINTS = 256; // relax
    constexpr int SP_NPOINTS = 512;  // spectral
    constexpr int CH_NPOINTS = 32;  // chase
    constexpr int EX_NPOINTS = 2048; // extraction
    constexpr int SH_NPOINTS = 2000; //shoot

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

