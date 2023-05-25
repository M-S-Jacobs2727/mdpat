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
        uint64_t nSteps = 0UL;
        StepRange(uint64_t iStep, uint64_t eStep, uint64_t dStep) : initStep(iStep), endStep(eStep), dumpStep(dStep)
        {
            nSteps = (eStep - iStep) / dStep + 1UL;
        }
    };

    template <typename T>
    struct Frame
    {
        uint64_t timestep = 0;
        uint64_t nAtoms = 0;
        uint32_t nCols = 0;
        T box[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        std::vector<std::string> columnLabels;
        std::vector<T> atoms;
    };

    Frame<double> readDumpBinary(std::string filepath);

    std::string getFilepath(
        const std::string & directory,
        const int step);

    std::string getFilename(const uint64_t step);

    int32_t findCoordColumn(
        const Frame<float> & frame,
        const char coord,
        const bool wrapped);
    
    int readAtomsSection(
        std::istream &dumpFile,
        std::vector<int> columnFlag,
        uint32_t nValidCols,
        uint64_t nAtoms,
        std::vector<float> &atoms);

    int readAtomsSectionByType(
        std::istream &dumpFile,
        std::vector<int> columnFlag,
        uint32_t nValidCols,
        std::vector<float> &atoms,
        std::vector<uint32_t> typeflag,
        uint64_t nValidAtoms);

    Frame<float> readDumpTextHeader(std::istream &dumpFile);

    void skipDumpTextHeader(std::istream &dumpFile);

    std::vector<float> getTrajectories(
        std::filesystem::path directory,
        const StepRange stepRange,
        const bool wrapped,
        const uint32_t atomType = 0,
        const uint32_t dim = 3);

}