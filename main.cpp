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
    Spacetime spacetime(Params::NX, Params::NY);

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

    // dev 
    ahfinder.initialize();
    ahfinder.update(spacetime);

    ///////////////////////////
    // Close Program
    
    std::cout << "\n*=================*\n";
    std::cout << "| Shutting Down ~ |\n";
    std::cout << "*=================*\n\n";

    return 0;
}

