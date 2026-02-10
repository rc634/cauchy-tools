#ifndef AHFINDER_HPP
#define AHFINDER_HPP

#include <vector>

class Spacetime;

class AHFinder {
public:
    // 1D horizon representation
    std::vector<double> sigma;     // parameter / angular coordinate
    std::vector<double> f;         // horizon radius
    std::vector<double> psi_sigma; // psi evaluated on horizon

    // points 
    int num_points;

    AHFinder(int npoints);

    void initialize();
    void update(const Spacetime& spacetime);

    void hello() const;
};

#endif

