#pragma once

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <regex>
#include <sstream>
#include <map>
#include <string>
#include <vector>

#include <mpi.h>

#include "bcastContainers.hpp"
#include "error.hpp"
#include "msd.hpp"   // add other analysis files as we write them
#include "stepRange.hpp"
#include "trajectory.hpp"

namespace MDPAT
{
    class InputReader
    {
    public:
        InputReader(const std::string &);
        ~InputReader();
        void runFile();
        std::vector<std::string> parseLine(const std::string &);
        void executeCommand(const std::vector<std::string> &);

        typedef void(*CommandPtr)(MDPAT::Trajectory&, const std::vector<std::string> &);
        std::map<char *, CommandPtr> commandMap {
            {"msd", &meanSquaredDisplacement},
        };

    private:
        std::filesystem::path inputFile;
        std::ifstream input;
        std::string dumpfileString;
        StepRange stepRange;
        MDPAT::Trajectory trajectory;

        double timestep = 0.0;
        int me;

        void trajCmd();
        void timestepCmd();
        void incorrectArgs(
            const std::string & command,
            const int expected_nargs,
            const int found_nargs);
        StepRange parseRange(const std::string &);
        StepRange parseRange(const std::string &, const uint64_t);
    };

}
/*=============================================================================
# `InputValues readInput(std::istream& stream)`
Reads input from `stream` consisting of keyword-value pairs (or more than one
value) separated by whitespace. The results are returned as a struct of these
values. Valid keywords are listed below:

## Dump file definitions
These define which dump files/timesteps to read. For now, filenames are assumed
to be `dump.<timestep>.txt`, where <timestep> is a 9-digit integer left-padded
with 0's.
* `directory`: The directory containing the dump files.
* `initStep`: The first step.
* `endStep`: The final step.
* `dumpStep`: The increment between snapshots to be read. E.g., a value of 10
with `initStep 0` would read frames 0, 10, 20, etc.
* `breakStep`: The step furthest away from `initStep` to consider for dynamic
analysis. E.g., a value of 100 with `initStep 0` and `dumpStep 10` would
perform analysis between frames 0:10, 0:20, 0:30, ..., 0:100, 10:20, ...,
10:110, 20:30, etc.
* `nAtoms`: The number of atoms in the dump files.
* `totalCols`: The number of columns in the dump files
* `atomTypes`: Following this keyword is one whole number indicating the number
of atom types to read from the files followed by that many atom types (whole
numbers).
* `columns`: As with `atomTypes`; one number to indicate the number of columns
to be allocated, and one number for each column to read. The ID column (must
be first) is column 0.

## Simulation-specific definitions
* `NN`: The number of atoms in a molecule. Generally used for the DP of a
coarse-grained polymer chain.
* `dim`: The number of spatial dimensions in the simulation.
* `dt`: The timestep of the simulation. Used with `dumpStep` to determine the
amount of simulation time between dump files.

## Scattering definitions
* `binFactor`: Indicates linear scaling of the scattering vector. The
scattering vectors are scaled by this factor, starting from `2*binFactor*pi/L`.
* `numBins`: Indicates log-scaling of the scattering vector. The scattering
vectors are log-scaled from `2*pi/L`.
=============================================================================*/