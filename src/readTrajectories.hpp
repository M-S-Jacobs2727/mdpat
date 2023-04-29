#pragma once

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

namespace MDPAT
{
struct StepRange
{
    uint64_t initStep = 0UL;
    uint64_t endStep = 0UL;
    uint64_t dumpStep = 0UL;
    StepRange(uint64_t iStep, uint64_t eStep, uint64_t dStep) :
        initStep(iStep), endStep(eStep), dumpStep(dStep) {}
};

template<typename T>
struct Frame {
    uint64_t timestep = 0;
    uint64_t nAtoms = 0;
    uint32_t nCols = 0;
    T box[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<std::string> columnLabels;
    std::vector<T> atoms;
};

std::string getFilepath(std::string directory, int step);

Frame<double> readDumpBinary(std::string filepath);

int readAtomsSection(std::istream & dumpFile, std::vector<int> columnFlag, uint32_t nValidCols, uint64_t nAtoms, std::vector<float> & atoms);

int readAtomsSectionByType(std::istream & dumpFile, std::vector<int> columnFlag, uint32_t nValidCols, std::vector<float> & atoms, std::vector<uint32_t> typeflag, uint64_t nValidAtoms);

int readDumpTextFrame(std::istream & dumpFile, std::vector<int>columnFlag, int nCols, Frame<float> * dump);

int readDumpTextFile(std::istream & dumpFile, std::vector<int> columnFlag, std::vector<Frame<float>*> & dumpFrames);

std::vector<Frame<float>*> readDumpFiles(
    int firstStep, int nSteps, int dumpStep, int nAtoms, int totalCols,
    std::vector<int> columns, std::string directory
);

std::vector<float> getTrajectories(bool wrapped,
                                    StepRange stepRange,
                                    std::filesystem::path directory,
                                    uint32_t atomType=0,
                                    uint32_t dim=3);

}