#include <iostream>

#include "params.hpp"
#include "Spacetime.hpp"
#include "DataLoader.hpp"
#include "AHFinder.hpp"
#include "Shooter.hpp"
#include "SpectralSolver.hpp"

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
    // spacetime.set_data_nakamura_prolate();
    // spacetime.set_data_nakamura_oblate();
    spacetime.set_data_BH();

    // loads initial data files e.g. .dat
    DataLoader loader;

    // loader.loadCSV(spacetime, "data/psi.dat", "data/W.dat");

    // apparent horizon finder object
    AHFinder ahfinder(Params::AH_NPOINTS);

    // pure integrating surface, too heavy for relaxation!
    AHFinder surface(Params::EX_NPOINTS);

    // RK4 shooter
    Shooter shooter(Params::SH_NPOINTS);

    // spectral solver
    SpectralSolver spectral(Params::SP_NPOINTS, Params::PL_LMAX);

   
    
    //////////////////////////
    // shouts
    spacetime.hello();
    loader.hello();
    ahfinder.hello();
    surface.hello();
    shooter.hello();
    spectral.hello();


    ///////////////////////////
    // Run Code
    
    std::cout << "\n";
    std::cout << "*======================*\n";
    std::cout << "| Test Surface <> ! |\n";
    std::cout << "*======================*\n\n";

    //////////////////////
    // dev 

    // initial radii
    double r_horizon = Params::RH;
    double r_extraction = Params::RX;
    double r_spectral = Params::RS;

    // data
    ahfinder.initialize(spacetime,r_horizon);
    surface.initialize(spacetime,r_extraction);
    shooter.initialize(r_horizon);
    spectral.initialize(spacetime, r_spectral);
    ahfinder.save("before");

    /////////
    // initial output from surface 
    surface.refresh(spacetime);
    // output from AHFinder
    std::cout << "\n* ============ Surface ============= *\n";
    std::cout << "| Surf R0 = " << surface.f[0] << "\n";
    std::cout << "| Surf RE = " << surface.f[ahfinder.num_points-1] << "\n";
    std::cout << "| Surf Area = " << surface.area() << "\n";
    std::cout << "| Surf Irreducible Mass = " << surface.mass_irr() << "\n";
    std::cout << "| Surf Radius = " << surface.r() << "\n";
    std::cout << "| Surf Eccentricity = " << surface.eccentricity() << "\n";
    std::cout << "| Surf Spin = " << surface.a_H() << "\n";
    std::cout << "| Surf Dimensionless Spin = " << surface.chi_H() << "\n";
    std::cout << "* ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ * \n";


    std::cout << "\n";
    std::cout << "*======================*\n";
    std::cout << "| Relaxation Method <> ! |\n";
    std::cout << "*======================*\n\n";

    ///////////////////////
    // relaxation iterations for AH solver          
    for (size_t i = 0; i < 1600001; i++)      
    {
        // ahfinder.refresh_Nakamura();
        ahfinder.refresh(spacetime);
        ahfinder.relax();
        if (i%100000==0) {
            // a breather step, we can do more expensice things here
            std::cout << " - * - Relax iter : " << i << ", ";
            std::cout << " r ~ " << ahfinder.r() << ", ";
            std::cout << " A ~ " << ahfinder.area() << ", ";
            std::cout << " psi ~ " << ahfinder.psi_h() << ", ";
            std::cout << " M ~ " << ahfinder.mass() << ", ";
            std::cout << " e ~ " << ahfinder.eccentricity() << ", ";
            std::cout << " res ~ " << ahfinder.res() << "\n";
        }
    }

    // output from AHFinder
    std::cout << "\n* ============ A H Finder ============= *\n";
    std::cout << "| Horizon R0 = " << ahfinder.f[0] << "\n";
    std::cout << "| Horizon RE = " << ahfinder.f[ahfinder.num_points-1] << "\n";
    std::cout << "| Horizon Area = " << ahfinder.area() << "\n";
    std::cout << "| Horizon Irreducible Mass = " << ahfinder.mass_irr() << "\n";
    std::cout << "| Horizon Radius = " << ahfinder.r() << "\n";
    std::cout << "| Horizon Eccentricity = " << ahfinder.eccentricity() << "\n";
    std::cout << "| Horizon Spin = " << ahfinder.a_H() << "\n";
    std::cout << "| Horizon Dimensionless Spin = " << ahfinder.chi_H() << "\n";
    std::cout << "* ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ * \n";

    // save relaxation solution
    ahfinder.save("relaxed");


    std::cout << "\n";
    std::cout << "*======================*\n";
    std::cout << "| Shooting Method <> ! |\n";
    std::cout << "*======================*\n\n";

    /////////////////
    // shooter 
    shooter.interval_bisection(spacetime);
    shooter.refresh(spacetime);

    // output from shooter
    std::cout << "\n* ============ Shooter ================ *\n";
    std::cout << "| Horizon R0 = " << shooter.f[0] << "\n";
    std::cout << "| Horizon RE = " << shooter.f[shooter.num_points-1] << "\n";
    std::cout << "| Horizon Area = " << shooter.area() << "\n";
    std::cout << "| Horizon Irreducible Mass = " << shooter.mass_irr() << "\n";
    std::cout << "| Horizon Radius = " << shooter.r() << "\n";
    std::cout << "| Horizon Eccentricity = " << shooter.eccentricity() << "\n";
    std::cout << "| Horizon Spin = " << shooter.a_H() << "\n";
    std::cout << "| Horizon Dimensionless Spin = " << shooter.chi_H() << "\n";
    std::cout << "* ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ * \n";

    // save shooting solution
    shooter.save("shot");

    // spectral stuff
    std::cout << "\n";
    std::cout << "*======================*\n";
    std::cout << "| Spectral Solver <> ! |\n";
    std::cout << "*======================*\n\n";
    spectral.refresh(spacetime);
    spectral.gradient_descent(spacetime);
    spectral.residual(); // to update the residual
    spectral.save("spectral");
    std::cout << "Spec res : " << spectral.res() << "\n";

    // output from shooter
    std::cout << "\n* ============ Spectral ================ *\n";
    std::cout << "| Horizon R0 = " << spectral.f[0] << "\n";
    std::cout << "| Horizon RE = " << spectral.f[spectral.num_points-1] << "\n";
    std::cout << "| Horizon Area = " << spectral.area() << "\n";
    std::cout << "| Horizon Res = " << spectral.res() << "\n";
    std::cout << "| Horizon Irreducible Mass = " << spectral.mass_irr() << "\n";
    std::cout << "| Horizon Radius = " << spectral.r() << "\n";
    std::cout << "| Horizon Eccentricity = " << spectral.eccentricity() << "\n";
    std::cout << "| Horizon Spin = " << spectral.a_H() << "\n";
    std::cout << "| Horizon Dimensionless Spin = " << spectral.chi_H() << "\n";
    std::cout << "* ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ * \n";

    

    ///////////////////////////
    // Close Program
    
    std::cout << "\n";
    std::cout << "*=================*\n";
    std::cout << "| Shutting Down ~ |\n";
    std::cout << "*=================*\n\n";

    return 0;
}

