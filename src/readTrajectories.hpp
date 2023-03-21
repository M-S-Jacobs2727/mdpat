#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

template<typename T>
struct Dump {
    long int nAtoms;
    int nCols;
    T box[6];
    std::vector<T> atoms;
};

std::string getFilename(std::string directory, int step);

template<typename T>
void readAtomsSection(std::ifstream &dumpfile, int nAtoms, std::vector<int> columnFlag, Dump<T> & dump);

Dump<double> readDumpBinary(std::string filepath);

template<typename T>
int readDumpText(const std::string & filepath, int nAtoms, std::vector<int> columnFlag, Dump<T> & dump);

template<typename T>
Dump<T> readDumpFiles(
    int firstStep, int nSteps, int dumpStep, int nAtoms, int totalCols,
    std::vector<int> columns, std::string directory
);
