#pragma once
#include "Surface.hpp"

class AHFRelax : public Surface {
public:
    AHFRelax(int npoints);

    // one relaxation step: evaluates the AH residual into dfdt, then steps f
    void relax();

    void hello() const;
};
