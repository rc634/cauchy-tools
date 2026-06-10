#pragma once
#include "Surface.hpp"

class EXSurf : public Surface {
public:
    EXSurf(int npoints);

    void set_sphere(double r);
    void set_ellipsoid(double a, double b);

    void extraction_output();
};
