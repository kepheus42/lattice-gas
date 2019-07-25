#include "output.hpp"

void vector_to_file(std::vector<int> vec, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        for(int n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(int));
        }
        ofs.close();
}

void vector_to_file(std::vector<double> vec, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        for(double n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(double));
        }
        ofs.close();
}
