#pragma once

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace MDPAT
{
template<typename T>
struct Frame {
    long timestep = 0;
    long nAtoms = 0;
    int nCols = 0;
    T box[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<T> atoms;
};

std::string getFilename(std::string directory, int step);

Frame<double> readDumpBinary(std::string filepath);

template<typename T>
int readAtomsSection(std::istream & dumpFile, std::vector<int> columnFlag, int nCols, long nAtoms, std::vector<T> & atoms);

template<typename T>
int readDumpTextFrame(std::istream & dumpFile, std::vector<int>columnFlag, int nCols, Frame<T> * dump);

template<typename T>
int readDumpTextFile(std::istream & dumpFile, std::vector<int> columnFlag, std::vector<Frame<T>*> & dumpFrames);

template<typename T>
Frame<T> readDumpFiles(
    int firstStep, int nSteps, int dumpStep, int nAtoms, int totalCols,
    std::vector<int> columns, std::string directory
);
}