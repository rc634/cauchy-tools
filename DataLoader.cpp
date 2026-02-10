#include "DataLoader.hpp"
#include "Spacetime.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

void DataLoader::loadCSV(const std::string& filename, Spacetime& spacetime) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    size_t index = 0;

    while (std::getline(file, line)) {
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

    std::cout << "Loaded psi field from " << filename
              << " (" << index << " values)" << std::endl;
}

void DataLoader::hello() const {
    std::cout << "Hello from DataLoader!" << std::endl;
}

