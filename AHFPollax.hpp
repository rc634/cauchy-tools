#pragma once
#include "Surface.hpp"
class Spacetime;

class AHFPollax : public Surface {
public:
    AHFPollax(int npoints);

    void relax();
    void refine(const Spacetime& spacetime);

    void hello() const;
};
