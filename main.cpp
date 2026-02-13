#include <iostream>

#include "params.hpp"
#include "Spacetime.hpp"
#include "DataLoader.hpp"
#include "AHFinder.hpp"

int main() {

    ///////////////////////////
    // Hello World!!
    
    std::cout << "*========================*\n";
    std::cout << "| Starting Cauchy Tools! |\n";
    std::cout << "*========================*\n\n";

    ////////////////////////////
    // Objects

    // spacetime holds 2d fields and 2d manifold information
    Spacetime spacetime;

    // loads initial data files e.g. .dat
    DataLoader loader;
    loader.loadCSV("data/psi.dat", spacetime);

    // apparent horizon finder object
    AHFinder ahfinder(Params::AH_NPOINTS);

    // pure integrating surface, too heavy for relaxation!
    AHFinder surface(Params::EX_NPOINTS);

   
    
    //////////////////////////
    // shouts
    spacetime.hello();
    loader.hello();
    ahfinder.hello();
    surface.hello();

    ///////////////////////////
    // Run Code
    
    std::cout << "\n";
    std::cout << "*======================*\n";
    std::cout << "| Finding Horizon <> ! |\n";
    std::cout << "*======================*\n\n";

    //////////////////////
    // dev 

    // initial radii
    double r_horizon = Params::RH;
    double r_extraction = Params::RX;

    // data
    ahfinder.initialize(spacetime,r_horizon);
    surface.initialize(spacetime,r_extraction);
    ahfinder.update(spacetime);
    ahfinder.save("before");

    // initial output from AHFinder
    std::cout << "Horizon Area ~ " << ahfinder.area() << "\n";
    std::cout << "Horizon Mass ~ " << ahfinder.mass() << "\n";
    std::cout << "Horizon Psi ~ " << ahfinder.psi_h() << "\n";
    std::cout << "Horizon Radius ~ " << ahfinder.r() << "\n";
    std::cout << "Horizon Res ~ " << ahfinder.res() << "\n";
    std::cout << "Horizon Eccentricity ~ " << ahfinder.eccentricity() << "\n";
    std::cout << "* ~ ~ * \n";

    // initial output from surface 
    std::cout << "Surface Area = " << surface.area_flat() << "\n";
    std::cout << "Surface psi = " << surface.psi_h() << "\n";
    std::cout << "Surface Misner-Sharp Mass = " << surface.mass_MS() << "\n";
    std::cout << "Surface Schwarzschild Mass = " << surface.mass_SC() << "\n";
    std::cout << "Surface Hawking Mass = " << surface.mass() << "\n";
    std::cout << "* ~ ~ * \n";

    // relaxation iterations for AH solver          
    for (size_t i = 0; i < 2000001; i++)
    {
        ahfinder.relax();
        ahfinder.refresh(spacetime);
        if (i%200000==0) {
            // a breather step, we can do more expensice things here
            std::cout << " - iter : " << i << ", ";
            std::cout << " r ~ " << ahfinder.r() << ", ";
            std::cout << " A ~ " << ahfinder.area() << ", ";
            std::cout << " psi ~ " << ahfinder.psi_h() << ", ";
            std::cout << " M ~ " << ahfinder.mass() << ", ";
            std::cout << " e ~ " << ahfinder.eccentricity() << ", ";
            std::cout << " res ~ " << ahfinder.res() << "\n";
            std::cout << " - ~ surf : Area ~ " << surface.area() 
              << " : Misner-Sharp Mass ~ " << surface.mass_MS() << "\n";
            std::cout << " * \n";
        }
    }

    // initial output from AHFinder
    std::cout << "Horizon Area ~ " << ahfinder.area() << "\n";
    std::cout << "Horizon Mass ~ " << ahfinder.mass() << "\n";
    std::cout << "Horizon Psi ~ " << ahfinder.psi_h() << "\n";
    std::cout << "Horizon Radius ~ " << ahfinder.r() << "\n";
    std::cout << "Horizon Res ~ " << ahfinder.res() << "\n";
    std::cout << "Horizon Eccentricity ~ " << ahfinder.eccentricity() << "\n";
    std::cout << "* ~ ~ * \n";

    // initial output from surface 
    std::cout << "Surface Area = " << surface.area_flat() << "\n";
    std::cout << "Surface psi = " << surface.psi_h() << "\n";
    std::cout << "Surface Misner-Sharp Mass = " << surface.mass_MS() << "\n";
    std::cout << "Surface Schwarzschild Mass = " << surface.mass_SC() << "\n";
    std::cout << "Surface Hawking Mass = " << surface.mass() << "\n";
    std::cout << "* ~ ~ * \n";

    // save after
    ahfinder.save("after");
    

    ///////////////////////////
    // Close Program
    
    std::cout << "\n";
    std::cout << "*=================*\n";
    std::cout << "| Shutting Down ~ |\n";
    std::cout << "*=================*\n\n";

    return 0;
}

