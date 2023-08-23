#include "readInput.hpp"

#include <sstream>

namespace fs = std::filesystem;
using std::string;
using std::vector;

namespace MDPAT
{
InputReader::InputReader(const string& inputFile) : 
    m_inputFile(inputFile)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &m_me);
    m_commandMap["msd"] = meanSquaredDisplacement;
}

InputReader::~InputReader() {}

void InputReader::runFile()
{
    if (!fs::is_regular_file(m_inputFile))
        errorAll(Error::IOERROR, "Input file %s does not exist.", m_inputFile.c_str());

    vector<string> lines;
    if (m_me == 0)
    {
        std::ifstream ifs(m_inputFile);
        if (ifs.fail())
            errorOne(Error::IOERROR, "Couldn't open input file: %s", m_inputFile.c_str());

        string line;
        while (!ifs.eof())
        {
            std::getline(ifs, line);
            lines.push_back(line);
        }
        ifs.close();
    }
    bcast(lines, 0, MPI_COMM_WORLD);

    for (const auto& line : lines)
    {
        const auto words = parseLine(line);
        if (words.size() > 1)
            executeCommand(words);
    }
}

vector<string> InputReader::parseLine(const string& line)
{
    vector<string> words;
    if (line.find_first_not_of(" \t\n") == string::npos)
        return words;

    string word;
    std::istringstream iss(line);
    iss >> word;

    if (m_commandMap.find(word.data()) == m_commandMap.end())
    {
        if (word == "traj") ;  // written this way because we may add more non-analysis commands
        else
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

void InputReader::executeCommand(const vector<string>& words)
{
    const string command = words[0];
    if (command == "trajectory")
    {
        trajCmd(words);
    }
    else if (m_commandMap.find(command.c_str()) != m_commandMap.end()) 
    {
        if (!m_trajectory.isLoaded)
            errorAll(Error::ARGUMENTERROR, "Command `%s` called without a loaded trajectory", command.c_str());
        const auto args = {words.begin()+1, words.end()};
        m_commandMap[command.c_str()](m_trajectory, args);
        // m_trajectory->reset();  // undo any permutation of the data?
    }
    else 
    {
        errorAll(Error::SYNTAXERROR, "Unknown command: %s", command.c_str());
    }
}

void InputReader::locateTrajFiles()
{
    auto subStartIdx = m_dumpfileString.find('%');
    if (subStartIdx == string::npos)
    {
        m_dumpfilePath = m_parentDir / m_dumpfileString;
        if (!std::filesystem::is_regular_file(m_dumpfilePath))
            errorAll(Error::IOERROR, "Could not locate dumpfile %s", m_dumpfilePath.c_str());
        return;
    }

    auto subEndIdx = m_dumpfileString.find_first_of("DIdi", subStartIdx+1);
    if (subEndIdx == string::npos)
        errorAll(Error::SYNTAXERROR, "Invalid dumpfile string: %s", m_dumpfileString.c_str());
    
    const char typeChar = m_dumpfileString[subEndIdx];

    string prefix = m_dumpfileString.substr(0, subStartIdx);
    string fmt = m_dumpfileString.substr(subStartIdx + 1, subEndIdx - subStartIdx - 1);
    string suffix = m_dumpfileString.substr(subEndIdx + 1);

    for (const char c : fmt)
        if (c < '0' || c > '9')
            errorAll(Error::SYNTAXERROR, "Invalid dumpfile string: %s", m_dumpfileString.c_str());
    
    uint32_t width = 0;
    bool fill = false;
    if (fmt.size() != 0)
    {
        fill = fmt[0] == '0';
        width = fill ? std::stoul(fmt.substr(1)) : std::stoul(fmt);
    }
    
    for (uint64_t step = m_stepRange.initStep; step <= m_stepRange.endStep; step += m_stepRange.dumpStep)
    {
        std::ostringstream oss;
        if (fill)
            oss.fill('0');
        oss << prefix;
        if (width)
            oss << std::setw(width);
        oss << step << suffix;

        fs::path filepath = m_parentDir / oss.str();
        if (!fs::is_regular_file(filepath))
            errorAll(Error::IOERROR, "Could not locate dumpfile %s", m_dumpfilePath.c_str());
        
        m_dumpfilePathsVec.push_back(filepath);
    }
}

void InputReader::trajCmd(const vector<string> &words)
{
    if (words.size() != 3)
        incorrectArgs(words[0], 2, words.size() - 1);

    fs::path tmp(words[1]);
    m_parentDir = tmp.parent_path();
    if (!fs::is_directory(m_parentDir))
        errorAll(Error::IOERROR, "Could not find directory for dumpfiles: %s", m_parentDir.c_str());

    m_dumpfileString = tmp.string();

    // TODO: make general for single-file dumpfiles and per-timestep dumpfiles
    // For single-file dumpfiles, `m_dumpfileString` will be a single file with no % signs.
    // These will have to be parsed step-by-step in `readTrajectories.cpp` with a check on the recorded timestep.
    // For per-timestep dumpfiles, `m_dumpfileString` will contain a % sign (substitution point)
    // and should be able to reference multiple files within a directory. 
    // This is probably easiest to do in `readTrajectories.cpp` completely, 
    // so that `m_dumpfileString` and `m_stepRange` will be smaller and easier to pass around.

    m_stepRange = parseRange(words[2]);

    locateTrajFiles();
    
    if (m_dumpfilePathsVec.size() != 0)
        m_trajectory.read(m_dumpfilePathsVec);
    else
        m_trajectory.read(m_dumpfilePath);
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
//TODO: Make these into other constructors for StepRange!
StepRange& InputReader::parseRange(const std::string& rangeString)
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

    StepRange stepRange;

    stepRange.initStep = std::stoul(rangeString.substr(0, hyphenPos - 1));
    stepRange.endStep = std::stoul(rangeString.substr(hyphenPos, pos - hyphenPos - 1));
    stepRange.dumpStep = std::stoul(rangeString.substr(pos));

    if (stepRange.endStep <= stepRange.initStep || stepRange.dumpStep == 0)
        errorAll(
            Error::ARGUMENTERROR,
            "Invalid range: %d-%d:%d\n"
            "End of range must be at least as large as beginning,"
            " and the step size must be greater than 0.",
            stepRange.initStep,
            stepRange.endStep,
            stepRange.dumpStep);
    
    return stepRange;
}

StepRange& InputReader::parseRange(const std::string& rangeString, const uint64_t dumpStep)
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

    StepRange stepRange;

    stepRange.initStep = std::stoul(rangeString.substr(0, pos - 1));
    stepRange.endStep = std::stoul(rangeString.substr(pos));
    stepRange.dumpStep = dumpStep;

    if (stepRange.endStep <= stepRange.initStep || stepRange.dumpStep == 0)
        errorAll(
            Error::ARGUMENTERROR,
            "Invalid range: %d-%d:%d\n"
            "End of range must be at least as large as beginning,"
            " and the step size must be greater than 0.",
            stepRange.initStep,
            stepRange.endStep,
            stepRange.dumpStep);
    
    return stepRange;
}
}