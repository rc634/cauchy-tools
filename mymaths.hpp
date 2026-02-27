#ifndef MYMATHS_HPP
#define MYMATHS_HPP

#include <vector>

namespace mymaths
{
    double LegendreP(int l, double x);

    double LegendreSeries(std::vector<double> &coeff, double x);
    double LegendreSeriesDiff1(std::vector<double> &coeff, double x);
}

#endif