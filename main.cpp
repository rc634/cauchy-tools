#include <iostream>

#include "params.hpp"
#include "Spacetime.hpp"
#include "DataLoader.hpp"
#include "Surface.hpp"
#include "AHFShooter.hpp"
#include "AHFRelax.hpp"
#include "AHFChase.hpp"
#include "AHFPollax.hpp"
#include "EXSurf.hpp"

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
    loader.loadCSV(spacetime, "data/psi.dat", "data/W.dat", "data/rho-raw.dat", "data/v-raw.dat", "data/rho-eff.dat", "data/v-eff.dat");

    /////////////////////////////////////////////////
    // initial radii
    double r_horizon = Params::RH;
    double r_extraction = Params::RX;
    double r_spectral = Params::RS;

    //////////////////////////////////////////////////////////////
    // extraction surfaces
    EXSurf exsurf(Params::EX_NPOINTS);

    // // flat space (no spacetime loaded)
    // exsurf.announce("Spherical Surface");
    // exsurf.set_sphere(1.0);
    // exsurf.cout_state();

    // exsurf.announce("Spheroidal Surface");
    // exsurf.set_ellipsoid(2.0, 1.0);
    // exsurf.cout_state();

    // same surfaces with spacetime metric
    exsurf.announce("Spherical Surface (curved)");
    exsurf.set_sphere(1.0);
    exsurf.refresh(spacetime);
    exsurf.extraction_output();
    exsurf.set_sphere(5.0);
    exsurf.refresh(spacetime);
    exsurf.extraction_output();
    exsurf.set_sphere(10.0);
    exsurf.refresh(spacetime);
    exsurf.extraction_output();
    exsurf.set_sphere(20.0);
    exsurf.refresh(spacetime);
    exsurf.extraction_output();
    exsurf.set_sphere(30.0);
    exsurf.refresh(spacetime);
    exsurf.extraction_output();

    // exsurf.announce("Spheroidal Surface (curved)");
    // exsurf.set_ellipsoid(2.0, 1.0);
    // exsurf.refresh(spacetime);
    // exsurf.cout_state();

    // //////////////////////////////////////////////////////////////

    // // AHF Shooter Method
    // AHFShooter shooter(Params::SH_NPOINTS);
    // shooter.announce("AHFShooter");
    // shooter.initialize(r_horizon);
    // shooter.interval_bisection(spacetime);
    // shooter.refresh(spacetime);
    // // for (size_t i = 0; i < 100000; ++i) {
    // //     shooter.refresh(spacetime);
    // //     shooter.relax();
    // // }
    // shooter.cout_state();
    // shooter.save("shot_surf");

    // //////////////////////////////////////////////////////////////

    // // AHF Relaxation Method
    // AHFRelax relaxation(Params::AH_NPOINTS);
    // relaxation.announce("AHFRelax");
    // relaxation.initialize(r_horizon);

    // for (size_t i = 0; i < 2400001; i++) {
    //     relaxation.refresh(spacetime);
    //     relaxation.relax();
    //     if (i % 240000 == 0) {
    //         std::cout << " - * - Relax iter : " << i
    //                   << ",  r ~ "   << relaxation.r()
    //                   << ",  A ~ "   << relaxation.area()
    //                   << ",  M ~ "   << relaxation.mass_irr()
    //                   << ",  e ~ "   << relaxation.eccentricity()
    //                   << ",  res ~ " << relaxation.res() << "\n";
    //     }
    // }

    // relaxation.cout_state();
    // relaxation.save("relaxed_surf");

    // //////////////////////////////////////////////////////////////

    // AHF Chase Method
    AHFChase chaser(Params::CH_NPOINTS);
    chaser.announce("AHFChase");
    chaser.initialize(r_horizon);

    for (size_t i = 0; i < 100001; i++) {
        chaser.refresh(spacetime);
        chaser.chase();
        if (i % 1000 == 0) chaser.smooth_edges(5);
        if (i % 10000 == 0) {
            std::cout << " - * - Chase iter : " << i
                      << ",  r ~ "   << chaser.r()
                      << ",  A ~ "   << chaser.area()
                      << ",  M ~ "   << chaser.mass_irr()
                      << ",  e ~ "   << chaser.eccentricity()
                      << ",  res ~ " << chaser.res() << "\n";
        }
    }

    chaser.mutate();

    for (size_t i = 0; i < 200001; i++) {
        chaser.refresh(spacetime);
        chaser.chase();
        if (i % 1000 == 0) chaser.smooth_edges(5);
        if (i % 20000 == 0) {
            std::cout << " - * - Chase (fine) iter : " << i
                      << ",  r ~ "   << chaser.r()
                      << ",  A ~ "   << chaser.area()
                      << ",  M ~ "   << chaser.mass_irr()
                      << ",  e ~ "   << chaser.eccentricity()
                      << ",  res ~ " << chaser.res() << "\n";
        }
    }

    chaser.mutate();
    
    for (size_t i = 0; i < 500001; i++) {
        chaser.refresh(spacetime);
        chaser.chase();
        if (i % 1000 == 0) chaser.smooth_edges(5);
        if (i % 50000 == 0) {
            std::cout << " - * - Chase (fine) iter : " << i
                      << ",  r ~ "   << chaser.r()
                      << ",  A ~ "   << chaser.area()
                      << ",  M ~ "   << chaser.mass_irr()
                      << ",  e ~ "   << chaser.eccentricity()
                      << ",  res ~ " << chaser.res() << "\n";
        }
    }

    chaser.cout_state();
    chaser.save("chased_surf");

    //////////////////////////////////////////////////////////////

    // AHF Pollax Method (polygrid relaxation)
    AHFPollax pollax(32);
    pollax.announce("AHFPollax");
    pollax.initialize(r_horizon);

    for (size_t i = 0; i < 100001; i++) {
        pollax.refresh(spacetime);
        pollax.relax();
        if (i % 10000 == 0) {
            std::cout << " - * - Pollax iter : " << i
                      << ",  r ~ "   << pollax.r()
                      << ",  A ~ "   << pollax.area()
                      << ",  M ~ "   << pollax.mass_irr()
                      << ",  e ~ "   << pollax.eccentricity()
                      << ",  res ~ " << pollax.res() << "\n";
        }
    }
    pollax.refine(spacetime);
    for (size_t i = 0; i < 500001; i++) {
        pollax.refresh(spacetime);
        pollax.relax();
        if (i % 100000 == 0) {
            std::cout << " - * - Pollax iter : " << i
                      << ",  r ~ "   << pollax.r()
                      << ",  A ~ "   << pollax.area()
                      << ",  M ~ "   << pollax.mass_irr()
                      << ",  e ~ "   << pollax.eccentricity()
                      << ",  res ~ " << pollax.res() << "\n";
        }
    }
    pollax.refine(spacetime);
    for (size_t i = 0; i < 1200001; i++) {
        pollax.refresh(spacetime);
        pollax.relax();
        if (i % 400000 == 0) {
            std::cout << " - * - Pollax iter : " << i
                      << ",  r ~ "   << pollax.r()
                      << ",  A ~ "   << pollax.area()
                      << ",  M ~ "   << pollax.mass_irr()
                      << ",  e ~ "   << pollax.eccentricity()
                      << ",  res ~ " << pollax.res() << "\n";
        }
    }
    // pollax.refine(spacetime);
    // for (size_t i = 0; i < 800001; i++) {
    //     pollax.refresh(spacetime);
    //     pollax.relax();
    //     if (i % 80000 == 0) {
    //         std::cout << " - * - Pollax iter : " << i
    //                   << ",  r ~ "   << pollax.r()
    //                   << ",  A ~ "   << pollax.area()
    //                   << ",  M ~ "   << pollax.mass_irr()
    //                   << ",  e ~ "   << pollax.eccentricity()
    //                   << ",  res ~ " << pollax.res() << "\n";
    //     }
    // }
    
    pollax.cout_state();
    pollax.save("pollax_surf");

    ///////////////////////////
    // Close Program
    
    std::cout << "\n";
    std::cout << "*=================*\n";
    std::cout << "| Shutting Down ~ |\n";
    std::cout << "*=================*\n\n";

    return 0;
}

