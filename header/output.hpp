#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <vector>
#include <string>
#include <ios>

// without header string
void vector_to_file(std::vector<int>, std::string);
void vector_to_file(std::vector<double>, std::string);
// with header string
void vector_to_file(std::vector<int>, std::string, std::string);
void vector_to_file(std::vector<double>, std::string, std::string);

#endif
