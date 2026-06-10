#pragma once
#include "Surface.hpp"

class AHFChase : public Surface {
public:
    AHFChase(int npoints);

    void chase();
    void mutate();

    void hello() const;
};
