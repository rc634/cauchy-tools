#include "mymaths.hpp"

namespace mymaths
{

double LegendreP(int l, double x)
{
    if (l == 0)
        return 1.0;

    if (l == 1)
        return x;

    double P_lm2 = 1.0;  // P0, P Legendre l minus 2
    double P_lm1 = x;    // P1, P Legendre l minus 1
    double P_l  = 0.0;  // output

    for (int n = 2; n <= l; ++n)
    {
        P_l = ((2.0 * n - 1.0) * x * P_lm1 - (n - 1.0) * P_lm2) / n;
        P_lm2 = P_lm1;
        P_lm1 = P_l;
    }

    return P_l;
}

double LegendreSeries(std::vector<double> &coeff, double x) {
    // first two expansions points 
    double out = coeff[0]  + x*coeff[1];

    double P_lm2 = 1.0;  // P0, P Legendre l minus 2
    double P_lm1 = x;    // P1, P Legendre l minus 1
    double P_l  = 0.0;  // output

    for (int n = 2; n < coeff.size(); ++n)
    {
        P_l = ((2.0 * n - 1.0) * x * P_lm1 - (n - 1.0) * P_lm2) / n;
        P_lm2 = P_lm1;
        P_lm1 = P_l;
        out += coeff[n] * P_l;
    }

    return out;
}

// dirst derivative
// d/dtheta P(cos(theta)) = dx/dtheta d/dx P(x)
double LegendreSeriesDiff1(std::vector<double> &coeff, double x) {
    // first two expansions points 
    double out = coeff[1];

    double P_lm2 = 1.0;  // P0, P Legendre l minus 2
    double P_lm1 = x;    // P1, P Legendre l minus 1
    double P_l  = 0.0;  // output

    for (int n = 2; n < coeff.size(); ++n)
    {
        P_l = ((2.0 * n - 1.0) * x * P_lm1 - (n - 1.0) * P_lm2) / n;

        out += coeff[n] * n * ( x * P_l - P_lm1) / sqrt(1. - x*x);

        P_lm2 = P_lm1;
        P_lm1 = P_l;
    }

    return out;
}


}