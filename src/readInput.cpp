#include "readInput.hpp"

namespace fs = std::filesystem;
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

InputReader::InputReader()
{
}

void InputReader::parse(std::istream input)
{
    string line;
    std::istringstream linestream;
    while (input)
    {
        std::getline(input, line);
        if (line.size() == 0)
            continue;
        linestream.str(line);
        linestream >> line;
        if (linestream.eof())
            continue;
        
        if (line == "timestep")
        {
            linestream >> timestep;
            if (timestep <= 0)
            {
                std::cerr << "Invalid value of timestep: " << timestep << '\n';
                throw timestep;
            }
        }
        else if (line == "dim")
        {
            linestream >> dim;
            if (dim == 0)
            {
                std::cerr << "dim=0 is not allowed.\n";
                throw -3;
            }
        }
        else if (line == "directory")
        {
            linestream >> directory;
            if (!fs::exists(directory))
            {
                std::cerr << "Could not find directory " << directory << '\n';
                throw directory.string();
            }
            directorySet = true;
        }
        else if (line == "stepRange")
        {
            std::getline(linestream, line, '-');
            stepRange.initStep = std::stoul(line);
            std::getline(linestream, line, ':');
            stepRange.endStep = std::stoul(line);
            linestream >> stepRange.dumpStep;

            if (stepRange.endStep <= stepRange.initStep || stepRange.dumpStep == 0)
            {
                std::cerr << "Invalid values of stepRange: " << stepRange.initStep << '-' << stepRange.endStep << ':' << stepRange.dumpStep;
                throw -1;
            }
            stepRangeSet = true;
        }
        else if (line == "msd")
        {
            if (!directorySet || !stepRangeSet || timestep == 0.0)
            {
                std::cerr << "Command 'msd' was called, but either 'timestep'," <<
                    " 'directory', or 'stepRange' was not set!\n";
                throw -2;
            }

            std::size_t beginPos;
            std::size_t commaPos;

            // Atom type
            uint32_t atomType;
            linestream >> atomType;

            // std::vector<uint32_t> typelist;

            // beginPos = 0;
            // commaPos = line.find(',');
            // while (commaPos != std::string::npos)
            // {
            //     typelist.push_back(std::stoul(line.substr(beginPos, commaPos)));
            //     beginPos = commaPos + 1;
            //     commaPos = line.find(',', beginPos);
            // }
            // typelist.push_back(std::stoul(line.substr(beginPos, line.size() - beginPos)));

            // std::sort(typelist.begin(), typelist.end());
            // auto newend = std::unique(typelist.begin(), typelist.end());
            // typelist.resize(newend - typelist.begin());

            // Gap range in steps
            linestream >> line;
            beginPos = 0;
            commaPos = line.find('-');

            uint64_t minGap = std::stoul(line.substr(0, commaPos));
            beginPos = commaPos + 1;
            uint64_t maxGap = std::stoul(line.substr(beginPos, line.size() - beginPos));

            // Output file
            linestream >> line;
            fs::path outfile(line);
            if (!fs::exists(outfile.parent_path()))
            {
                std::cerr << "Can't write to " << outfile << '\n';
                throw outfile.string();
            }

            // Run
            msd(outfile,
                directory,
                stepRange,
                timestep,
                atomType,
                minGap,
                maxGap,
                dim);
        }
        
    }
}

InputReader::~InputReader()
{

}



}