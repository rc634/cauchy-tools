#ifndef DATALOADER_HPP
#define DATALOADER_HPP

#include <string>

class Spacetime;

class DataLoader {
public:
    DataLoader() = default;

    void loadCSV(const std::string& filename, Spacetime& spacetime);

    void hello() const;
};

#endif

