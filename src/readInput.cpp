#include "readInput.hpp"

namespace fs = std::filesystem;
using std::string;
using std::vector;

namespace MDPAT
{
    InputReader::InputReader(const string & inputFile) : 
        inputFile(inputFile)
    {
        MPI_Comm_rank(MPI_COMM_WORLD, &me);
    }

    InputReader::~InputReader()
    {}

    void InputReader::runFile()
    {
        if (!fs::is_regular_file(inputFile))
            errorAll(Error::IOERROR, "Input file %s does not exist.", inputFile.c_str());

        if (me == 0)
        {
            input.open(inputFile);
            if (input.bad())
                errorOne(Error::IOERROR, "Couldn't open input file: %s", inputFile.c_str());
        }

        bool eof = false;
        string line;
        while (!eof)
        {
            if (me == 0)
            {
                std::getline(input, line);
                eof = input.eof();
            }
            bcast(line, 0, MPI_COMM_WORLD);
            const auto words = parseLine(line);
            if (words.size() > 0UL)
                executeCommand(words);
            MPI_Bcast(&eof, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
        }

        input.close();
    }

    vector<string> InputReader::parseLine(const string & line)
    {
        vector<string> words;
        if (line.find_first_not_of(" \t\n") == string::npos)
            return words;

        string word;
        std::istringstream iss(line);
        iss >> word;

        if (commandMap.find(word.data()) == commandMap.end())
        {
            errorAll(Error::SYNTAXERROR, "Command not recognized: %s", word.c_str());
        }

        words.push_back(word);
        while (!iss.eof())
        {
            iss >> word;

            auto pos = word.find('"');
            if (pos == 0UL)
            {
                pos = word.find('"', pos);
                if (pos == word.size() - 1)
                {
                    words.push_back(word.substr(1, word.size() - 2));
                }
                else if (pos != string::npos)
                {
                    errorAll(Error::SYNTAXERROR, "Invalid syntax: %s", line);
                }
                else
                {
                    std::ostringstream oss(word.substr(1));
                    while (true)
                    {
                        iss >> word;
                        pos = word.find('"');
                        if (pos == string::npos)
                        {
                            oss << word;
                        }
                        else if (pos == word.size() - 1)
                        {
                            oss << word.substr(0, word.size() - 1);
                            break;
                        }
                        else
                        {
                            errorAll(Error::SYNTAXERROR, "Invalid syntax: %s", line);
                        }

                        if (iss.eof())
                        {
                            errorAll(Error::SYNTAXERROR, "Unmatched quotation mark: %s", line);
                        }
                    }

                    words.push_back(oss.str());
                }
                continue;
            }
            else if (pos != string::npos)
            {
                errorAll(Error::SYNTAXERROR, "Invalid syntax: %s", line);
            }

            pos = word.find('#');
            if (pos == 0UL)
                break;
            else if (pos != string::npos)
            {
                words.push_back(word.substr(0, pos));
                break;
            }

            words.push_back(word);
        }
        return words;
    }

    void InputReader::executeCommand(const vector<string> & words)
    {
        const string command = words[0];
        if (command == "trajectory") trajCmd();
        else if (commandMap.find(command) != commandMap.end()) 
        {
            if (!trajectory.isLoaded)
                errorAll(Error::ARGUMENTERROR, "Command `%s` called without a loaded trajectory", command);
            const auto args = {words.begin()+1, words.end()};
            commandMap[command](trajectory, args);
            // trajectory->reset();  // undo any permutation of the data
        }
        else errorAll(Error::SYNTAXERROR, "Unknown command: %s", command);
    }

    void InputReader::trajCmd(const vector<string> & words)
    {
        if (words.size() != 3)
            incorrectArgs(words[0], 1, words.size() - 1);

        fs::path tmp(words[1]);
        fs::path directory = tmp.parent_path();
        if (!fs::is_directory(directory))
            errorAll(Error::IOERROR, "Could not find directory for dumpfiles: %s", directory.c_str());

        dumpfileString = tmp.string();

        // TODO: make general for single-file dumpfiles and per-timestep dumpfiles
        // For single-file dumpfiles, `dumpfileString` will be a single file with no % signs.
        // These will have to be parsed step-by-step in `readTrajectories.cpp` with a check on the recorded timestep.
        // For per-timestep dumpfiles, `dumpfileString` will contain a % sign (substitution point)
        // and should be able to reference multiple files within a directory. 
        // This is probably easiest to do in `readTrajectories.cpp` completely, 
        // so that `dumpfileString` and `stepRange` will be smaller and easier to pass around.

        stepRange = parseRange(words[2]);

        if (stepRange.endStep <= stepRange.initStep || stepRange.dumpStep == 0)
            errorAll(
                Error::ARGUMENTERROR,
                "Invalid range: %s\n"
                "End of range must be at least as large as beginning,"
                " and the step size must be greater than 0.",
                words[2]);
        
        trajectory.read(directory, dumpfileString, stepRange);
    }

    void InputReader::incorrectArgs(
        const string & command,
        const int expected_nargs,
        const int found_nargs)
    {
        errorAll(
            Error::ARGUMENTERROR, 
            "Incorrect number of args for command %s: expected %d, found %d",
            command.c_str(),
            expected_nargs,
            found_nargs);
    }
    
    StepRange InputReader::parseRange(
        const std::string & rangeString)
    {
        const char message[] = "Invalid range syntax: %s\nMust be of form <init>-<end>:<dump>, e.g., 0-1000000:1000";

        auto pos = rangeString.find_first_not_of("0123456789-:");
        if (pos != string::npos)
            errorAll(Error::SYNTAXERROR, message, rangeString);
        
        auto count = std::count(rangeString.begin(), rangeString.end(), '-');
        if (count != 1)
            errorAll(Error::SYNTAXERROR, message, rangeString);
        
        auto count = std::count(rangeString.begin(), rangeString.end(), ':');
        if (count != 1)
            errorAll(Error::SYNTAXERROR, message, rangeString);

        pos = rangeString.find_first_of(":-");
        if (pos == 0UL)
            errorAll(Error::SYNTAXERROR, message, rangeString);
        if (rangeString[pos] == ':')
            errorAll(Error::SYNTAXERROR, message, rangeString);
        auto hyphenPos = pos;
        pos = rangeString.find(':');

        auto range = StepRange(
            std::stoul(rangeString.substr(0, hyphenPos - 1)),
            std::stoul(rangeString.substr(hyphenPos, pos - hyphenPos - 1)),
            std::stoul(rangeString.substr(pos))
        );

        return range;
    }

    StepRange InputReader::parseRange(
        const std::string & rangeString,
        const uint64_t dumpStep)
    {
        const char message[] = "Invalid range syntax: %s\nMust be of form <init>-<end>, e.g., 0-1000000";

        auto pos = rangeString.find_first_not_of("0123456789-");
        if (pos != string::npos)
            errorAll(Error::SYNTAXERROR, message, rangeString);
        
        auto count = std::count(rangeString.begin(), rangeString.end(), '-');
        if (count != 1)
            errorAll(Error::SYNTAXERROR, message, rangeString);

        pos = rangeString.find('-');
        if (pos == 0UL || pos == rangeString.size() - 1UL)
            errorAll(Error::SYNTAXERROR, message, rangeString);

        auto range = StepRange(
            std::stoul(rangeString.substr(0, pos - 1)),
            std::stoul(rangeString.substr(pos)),
            dumpStep
        );

        return range;
    }


}