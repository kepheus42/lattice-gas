#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <vector>
#include <string>
#include <ios>

template <typename T>
void vector_to_file(std::vector<T> vec, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        for(T n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}

template <typename T>
void vector_to_file(std::vector<T> vec, std::string header, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        // check for newline after header maybe?
        // ofs.write(header);
        for(T n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}

#endif
