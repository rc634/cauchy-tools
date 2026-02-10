#ifndef AHFINDER_HPP
#define AHFINDER_HPP

#include <vector>

class Spacetime;

class AHFinder {
public:
    // 1D horizon representation
    std::vector<double> sigma;     // parameter / angular coordinate
    std::vector<double> f;         // horizon radius
    std::vector<double> psi; // psi evaluated on horizon
    std::vector<double> dfdt; // for relaxation time updates

    // points 
    int num_points;
    double dt; // relaxation time stepper
    double ds; // delta parameter/angle

    AHFinder(int npoints);

    // returns horizon area
    double area();
    
    // approx horizon mass
    double mass();

    // partial derivs
    double d(const std::vector<double> &field, int i);
    double d2(const std::vector<double> &field, int i);

    void initialize(const Spacetime& spacetime);
    void update(const Spacetime& spacetime);
    void relax();

    void hello() const;
};

#endif

