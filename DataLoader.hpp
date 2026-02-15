#ifndef DATALOADER_HPP
#define DATALOADER_HPP

#include <string>

class Spacetime;

class DataLoader {
public:
    DataLoader() = default;

    void loadCSV(Spacetime& spacetime, 
                const std::string& filename1, 
                const std::string& filename2);

    void hello() const;
};

#endif

