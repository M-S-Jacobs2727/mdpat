#include "readInput.hpp"

using std::string;

namespace MDPAT
{
InputValues readInput(std::istream& stream) {
    InputValues input;
    string line;
    while (!stream.eof()) {
        stream >> line;
        if (line.empty())
            break;

        // Mandatory arguments
        if (line == "initStep")
            stream >> input.initStep;
        else if (line == "endStep")
            stream >> input.endStep;
        else if (line == "dumpStep")
            stream >> input.dumpStep;
        else if (line == "directory")
            stream >> input.directory;
        else if (line == "nAtoms")
            stream >> input.nAtoms;
        else if (line == "totalCols")
            stream >> input.totalCols;
        // Optional arguments
        else if (line == "breakStep")
            stream >> input.breakStep;
        else if (line == "atomTypes") {
            std::getline(stream, line);
            std::istringstream myline(line);
            int tmp;
            myline >> tmp;
            while (!myline && tmp > 0) {
                input.atomTypes.push_back(tmp);
                myline >> tmp;
            }
        }
        else if (line == "columns") {
            std::getline(stream, line);
            std::istringstream myline(line);
            int tmp;
            myline >> tmp;
            while (!myline && tmp > 0) {
                input.columns.push_back(tmp);
                myline >> tmp;
            }
        }
        else if (line == "breakStep")
            stream >> input.breakStep;
        else if (line == "dt")
            stream >> input.dt;
        else if (line == "dim")
            stream >> input.dim;
        else if (line == "NN")
            stream >> input.NN;
        else if (line == "binFactor")
            stream >> input.binFactor;
        else if (line == "numBins")
            stream >> input.numBins;
        else if (line == "vectorFile")
            stream >> input.vectorFile;
    }
    if (input.initStep < 0 || input.endStep <= 0 || input.dumpStep <= 0 || input.directory == "" || input.nAtoms <= 0 || input.totalCols <= 0)
        std::cerr << "Invalid or missing settings in input file!" << std::endl;
    return input;
}
}