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

   
    
    //////////////////////////
    // shouts
    spacetime.hello();
    loader.hello();
    ahfinder.hello();

    ///////////////////////////
    // Run Code
    
    std::cout << "\n";
    std::cout << "*======================*\n";
    std::cout << "| Finding Horizon <> ! |\n";
    std::cout << "*======================*\n\n";

    // dev 
    ahfinder.initialize(spacetime);
    ahfinder.update(spacetime);
    ahfinder.save("before");

    std::cout << "Horizon Area ~ " << ahfinder.area() << "\n";
    std::cout << "Horizon Mass ~ " << ahfinder.mass() << "\n";
    std::cout << "Horizon Psi ~ " << ahfinder.psi_h() << "\n";
    std::cout << "Horizon Radius ~ " << ahfinder.r() << "\n";
    std::cout << "Horizon Res ~ " << ahfinder.res() << "\n";
    for (size_t i = 0; i < 200001; i++)
    {
        ahfinder.relax();
        ahfinder.refresh(spacetime);
        if (i%10000==0) {
            std::cout << " - iter : " << i << ", ";
            std::cout << " r ~ " << ahfinder.r() << ", ";
            std::cout << " A ~ " << ahfinder.area() << ", ";
            std::cout << " psi ~ " << ahfinder.psi_h() << ", ";
            std::cout << " M ~ " << ahfinder.mass() << ", ";
            std::cout << " res ~ " << ahfinder.res() << "\n";
        }
    }

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

