#include "AHFChase.hpp"
#include "params.hpp"
#include <iostream>
#include <vector>
#include <cmath>

AHFChase::AHFChase(int npoints) : Surface(npoints) {}

void AHFChase::chase() {
    for (int i = 0; i < num_points; ++i)
        dfdt[i] = div_s(i);
    for (int i = 0; i < num_points; ++i)
        f[i] -= 0.25 * f[i]*f[i] * ds*ds * dfdt[i];
}

void AHFChase::mutate() {
    int    N_old  = num_points;
    double ds_old = ds;
    std::vector<double> f_old    = f;
    std::vector<double> sigma_old = sigma;

    int    N_new  = 2 * N_old;
    double ds_new = ds_old / 2.0;

    sigma.resize(N_new);
    f.resize(N_new);
    for (int i = 0; i < N_new; ++i) {
        double s  = (i + 0.5) * ds_new;
        sigma[i]  = s;
        double j  = s / ds_old - 0.5;
        int    j0 = static_cast<int>(j);
        if (j0 < 0)          j0 = 0;
        if (j0 > N_old - 2)  j0 = N_old - 2;
        double t  = j - j0;
        f[i] = (1.0 - t) * f_old[j0] + t * f_old[j0 + 1];
    }

    psi.resize(N_new,     0.0);
    W.resize(N_new,       0.0);
    dpsi_dR.resize(N_new, 0.0);
    dW_dR.resize(N_new,   0.0);
    dfdt.resize(N_new,    0.0);

    num_points = N_new;
    ds         = ds_new;
    dt         = Params::CFL * ds_new * ds_new;  // dt ~ ds^2 for stability

    smooth(5);  // smooth f then refresh(*st) internally

    std::cout << " - * - Chase mutate : N " << N_old << " -> " << N_new
              << ",  ds " << ds_old << " -> " << ds_new << "\n";
}

void AHFChase::hello() const {
    std::cout << "Hello from AHFChase!\n"
              << " - Number of horizon points = " << num_points << "\n";
}
