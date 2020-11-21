#include "output.hpp"

// Debugging code
#ifdef DEBUG
#define D(x) (x)
#else
#define D(x) do {} while(0)
#endif
// To add debugging messages, use D(std::cerr << "Debugging message 1 2 3!" << std::endl; )

// - - - - - - - - - - - - - - - - - - - - - - - -
void vector_to_file(std::vector<int> vec, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        for(int n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void vector_to_file(std::vector<unsigned long> vec, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        for(unsigned long n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void vector_to_file(std::vector<double> vec, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        for(double n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
// W I T H  H E A D E R  S T R I N G
// - - - - - - - - - - - - - - - - - - - - - - - -
void vector_to_file(std::vector<int> vec, std::string header, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        // check for newline after header maybe?
        // ofs.write(header);
        for(int n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}
// - - - - - - - - - - - - - - - - - - - - - - - -
void vector_to_file(std::vector<double> vec, std::string header, std::string fname){
        std::ofstream ofs;
        ofs.open(fname, std::ios::out | std::ios::binary);
        // check for newline after header maybe?
        // ofs.write(header);
        for(double n : vec) {
                ofs.write(reinterpret_cast<const char*>(&n), sizeof(n));
        }
        ofs.close();
}
