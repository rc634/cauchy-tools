#include "AHFPollax.hpp"
#include "Spacetime.hpp"
#include "params.hpp"
#include <iostream>
#include <cmath>
#include <vector>

AHFPollax::AHFPollax(int npoints) : Surface(npoints) {}

void AHFPollax::relax() {
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

void AHFPollax::refine(const Spacetime& spacetime) {
    int    N_old = num_points;
    double ds_old = ds;
    std::vector<double> f_old    = f;
    std::vector<double> sigma_old = sigma;

    int    N_new  = 2 * N_old;
    double ds_new = ds_old / 2.0;
    double dt_new = Params::CFL * ds_new * ds_new;

    // rebuild sigma and interpolate f onto the finer grid
    sigma.resize(N_new);
    f.resize(N_new);
    for (int i = 0; i < N_new; ++i) {
        double s = (i + 0.5) * ds_new;
        sigma[i] = s;
        // find bracket in old cell-centred grid: s = (j+0.5)*ds_old => j = s/ds_old - 0.5
        double j  = s / ds_old - 0.5;
        int    j0 = static_cast<int>(j);
        if (j0 < 0)           j0 = 0;
        if (j0 > N_old - 2)   j0 = N_old - 2;
        double t  = j - j0;
        f[i] = (1.0 - t) * f_old[j0] + t * f_old[j0 + 1];
    }

    // resize scratch arrays (refresh will fill them)
    psi.resize(N_new,     0.0);
    W.resize(N_new,       0.0);
    dpsi_dR.resize(N_new, 0.0);
    dW_dR.resize(N_new,   0.0);
    dfdt.resize(N_new,    0.0);

    num_points = N_new;
    ds         = ds_new;
    dt         = dt_new;

    refresh(spacetime);
    std::cout << " - * - Pollax refine : N " << N_old << " -> " << N_new
              << ",  ds " << ds_old << " -> " << ds_new << "\n";
}

void AHFPollax::hello() const {
    std::cout << "Hello from AHFPollax!\n"
              << " - Number of horizon points = " << num_points << "\n";
}
