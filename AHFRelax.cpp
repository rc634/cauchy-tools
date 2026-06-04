#include "AHFRelax.hpp"
#include "params.hpp"
#include <iostream>
#include <cmath>

AHFRelax::AHFRelax(int npoints) : Surface(npoints) {}

void AHFRelax::relax() {
    double term1, term2, term3, dfds, R;
    for (int i = 0; i < num_points; ++i) {
        dfds  = d(f, i);
        R     = f[i];
        term1 = (dfds + pow(dfds,3)/R/R) * (cos(sigma[i])/sin(sigma[i]) + 4.*d(psi,i)/psi[i]);
        term2 = -dfds * dfds * (4.*dpsi_dR[i]/psi[i] + 3./R);
        term3 = -2.*R - 4.*R*R*dpsi_dR[i]/psi[i];
        dfdt[i] = d2(f,i) + term1 + term2 + term3;
    }
    for (int i = 0; i < num_points; ++i)
        f[i] += dfdt[i] * dt;
}

void AHFRelax::hello() const {
    std::cout << "Hello from AHFRelax!\n"
              << " - Number of horizon points = " << num_points << "\n";
}
