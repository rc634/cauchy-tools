#include <iostream>

#include "params.hpp"
#include "Spacetime.hpp"
#include "DataLoader.hpp"
#include "Surface.hpp"
#include "AHFShooter.hpp"
#include "AHFRelax.hpp"

int main() {

    ///////////////////////////
    // Hello World!!
    
    std::cout << "\n";
    std::cout << "*========================*\n";
    std::cout << "| Starting Cauchy Tools! |\n";
    std::cout << "*========================*\n\n";

    ////////////////////////////
    // Objects

    // spacetime holds 2d fields and 2d manifold information
    Spacetime spacetime;

    // nakamura comparisons
    // spacetime.set_data_nakamura_prolate();
    // spacetime.set_data_nakamura_oblate();

    // perturbed BH test (controlled by small param epsilon internally)
    // spacetime.set_data_BH();

    // loads initial data files e.g. .dat
    DataLoader loader;
    loader.loadCSV(spacetime, "data/psi.dat", "data/W.dat");

    /////////////////////////////////////////////////
    // initial radii
    double r_horizon = Params::RH;
    double r_extraction = Params::RX;
    double r_spectral = Params::RS;

    //////////////////////////////////////////////////////////////
    // tests of pure ellipsoidal surfaces
    // raw test surface
    Surface surface(Params::EX_NPOINTS);

    // flat space (no spacetime loaded)
    surface.announce("Spherical Surface");
    surface.initialize(1.0);
    surface.cout_state();

    surface.announce("Spheroidal Surface");
    surface.initialize(2.0, 1.0);
    surface.cout_state();

    // same surfaces with spacetime metric
    surface.announce("Spherical Surface (curved)");
    surface.initialize(1.0);
    surface.refresh(spacetime);
    surface.cout_state();

    surface.announce("Spheroidal Surface (curved)");
    surface.initialize(2.0, 1.0);
    surface.refresh(spacetime);
    surface.cout_state();

    //////////////////////////////////////////////////////////////

    // AHF Shooter Method
    AHFShooter shooter(Params::SH_NPOINTS);
    shooter.announce("AHFShooter");
    shooter.initialize(r_horizon);
    shooter.interval_bisection(spacetime);
    shooter.refresh(spacetime);
    shooter.cout_state();
    shooter.save("shot_new");

    //////////////////////////////////////////////////////////////

    // AHF Relaxation Method
    AHFRelax relaxation(Params::AH_NPOINTS);
    relaxation.announce("AHFRelax");
    relaxation.initialize(r_horizon);

    for (size_t i = 0; i < 500001; i++) {
        relaxation.refresh(spacetime);
        relaxation.relax();
        if (i % 50000 == 0) {
            std::cout << " - * - Relax iter : " << i
                      << ",  r ~ "   << relaxation.r()
                      << ",  A ~ "   << relaxation.area()
                      << ",  M ~ "   << relaxation.mass_irr()
                      << ",  e ~ "   << relaxation.eccentricity()
                      << ",  res ~ " << relaxation.res() << "\n";
        }
    }

    relaxation.cout_state();
    relaxation.save("relaxed_new");



    ///////////////////////////
    // Close Program
    
    std::cout << "\n";
    std::cout << "*=================*\n";
    std::cout << "| Shutting Down ~ |\n";
    std::cout << "*=================*\n\n";

    return 0;
}

