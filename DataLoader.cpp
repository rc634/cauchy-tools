#include "DataLoader.hpp"
#include "Spacetime.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

void DataLoader::loadCSV(Spacetime& spacetime,
                        const std::string& filename1,
                        const std::string& filename2,
                        const std::string& filename3,
                        const std::string& filename4,
                        const std::string& filename5,
                        const std::string& filename6) {

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



    if (!filename3.empty()) {
        std::ifstream file3(filename3);
        if (!file3.is_open())
            throw std::runtime_error("Could not open file: " + filename3);

        index = 0;
        while (std::getline(file3, line)) {
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                if (index >= spacetime.src.size())
                    throw std::runtime_error("Too much data in src file (exceeds grid size)");
                try {
                    spacetime.src[index] = std::stod(value);
                } catch (const std::out_of_range&) {
                    spacetime.src[index] = 0.0;  // subnormal underflow → treat as zero
                }
                ++index;
            }
        }
        if (index != spacetime.src.size())
            throw std::runtime_error(
                "Not enough data in src file: expected " +
                std::to_string(spacetime.src.size()) +
                ", got " + std::to_string(index));
        std::cout << "Loaded src field from " << filename3
                  << " (" << index << " values)" << std::endl;
        file3.close();
    }

    if (!filename4.empty()) {
        std::ifstream file4(filename4);
        if (!file4.is_open())
            throw std::runtime_error("Could not open file: " + filename4);

        index = 0;
        while (std::getline(file4, line)) {
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                if (index >= spacetime.srcW.size())
                    throw std::runtime_error("Too much data in srcW file (exceeds grid size)");
                try {
                    spacetime.srcW[index] = std::stod(value);
                } catch (const std::out_of_range&) {
                    spacetime.srcW[index] = 0.0;
                }
                ++index;
            }
        }
        if (index != spacetime.srcW.size())
            throw std::runtime_error(
                "Not enough data in srcW file: expected " +
                std::to_string(spacetime.srcW.size()) +
                ", got " + std::to_string(index));
        std::cout << "Loaded srcW field from " << filename4
                  << " (" << index << " values)" << std::endl;
        file4.close();
    }

    if (!filename5.empty()) {
        std::ifstream file5(filename5);
        if (!file5.is_open())
            throw std::runtime_error("Could not open file: " + filename5);

        index = 0;
        while (std::getline(file5, line)) {
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                if (index >= spacetime.src_eff.size())
                    throw std::runtime_error("Too much data in src_eff file (exceeds grid size)");
                try {
                    spacetime.src_eff[index] = std::stod(value);
                } catch (const std::out_of_range&) {
                    spacetime.src_eff[index] = 0.0;
                }
                ++index;
            }
        }
        if (index != spacetime.src_eff.size())
            throw std::runtime_error(
                "Not enough data in src_eff file: expected " +
                std::to_string(spacetime.src_eff.size()) +
                ", got " + std::to_string(index));
        std::cout << "Loaded src_eff field from " << filename5
                  << " (" << index << " values)" << std::endl;
        file5.close();
    }

    if (!filename6.empty()) {
        std::ifstream file6(filename6);
        if (!file6.is_open())
            throw std::runtime_error("Could not open file: " + filename6);

        index = 0;
        while (std::getline(file6, line)) {
            std::stringstream ss(line);
            std::string value;
            while (std::getline(ss, value, ',')) {
                if (index >= spacetime.srcW_eff.size())
                    throw std::runtime_error("Too much data in srcW_eff file (exceeds grid size)");
                try {
                    spacetime.srcW_eff[index] = std::stod(value);
                } catch (const std::out_of_range&) {
                    spacetime.srcW_eff[index] = 0.0;
                }
                ++index;
            }
        }
        if (index != spacetime.srcW_eff.size())
            throw std::runtime_error(
                "Not enough data in srcW_eff file: expected " +
                std::to_string(spacetime.srcW_eff.size()) +
                ", got " + std::to_string(index));
        std::cout << "Loaded srcW_eff field from " << filename6
                  << " (" << index << " values)" << std::endl;
        file6.close();
    }

    file1.close();
    file2.close();
}

void DataLoader::hello() const {
    std::cout << "Hello from DataLoader!" << std::endl;
}

