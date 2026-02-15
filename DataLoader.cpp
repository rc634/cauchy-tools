#include "DataLoader.hpp"
#include "Spacetime.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

void DataLoader::loadCSV(Spacetime& spacetime, 
                        const std::string& filename1, 
                        const std::string& filename2) {

    std::ifstream file1(filename1);
    if (!file1.is_open()) {
        throw std::runtime_error("Could not open file: " + filename1);
    }

    std::string line;
    size_t index = 0;

    while (std::getline(file1, line)) {
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            if (index >= spacetime.psi.size()) {
                throw std::runtime_error(
                    "Too much data in psi file (exceeds grid size)"
                );
            }

            spacetime.psi[index] = std::stod(value);
            //std::cout << "DataLoader : " << std::stod(value) << "\n";
            ++index;
        }
    }

    if (index != spacetime.psi.size()) {
        throw std::runtime_error(
            "Not enough data in psi file: expected " +
            std::to_string(spacetime.psi.size()) +
            ", got " + std::to_string(index)
        );
    }

    std::cout << "Loaded psi field from " << filename1
              << " (" << index << " values)" << std::endl;


    std::ifstream file2(filename2);

    if (!file2.is_open()) {
        throw std::runtime_error("Could not open file: " + filename2);
    }

    index = 0;

    while (std::getline(file2, line)) {
        std::stringstream ss(line);
        std::string value;

        while (std::getline(ss, value, ',')) {
            if (index >= spacetime.W.size()) {
                throw std::runtime_error(
                    "Too much data in W file (exceeds grid size)"
                );
            }

            spacetime.W[index] = std::stod(value);
            //std::cout << "DataLoader : " << std::stod(value) << "\n";
            ++index;
        }
    }

    if (index != spacetime.W.size()) {
        throw std::runtime_error(
            "Not enough data in W file: expected " +
            std::to_string(spacetime.W.size()) +
            ", got " + std::to_string(index)
        );
    }

    std::cout << "Loaded W field from " << filename2
              << " (" << index << " values)" << std::endl;



    file1.close();
    file2.close();
}

void DataLoader::hello() const {
    std::cout << "Hello from DataLoader!" << std::endl;
}

