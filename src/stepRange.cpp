#include "stepRange.hpp"

#include <algorithm>
#include <string>
#include <sstream>

#include "error.hpp"

using std::string;

namespace MDPAT
{
StepRange::StepRange() : initStep(0UL), endStep(0UL), dumpStep(0UL), nSteps(0UL) {}
StepRange::StepRange(uint64_t iStep, uint64_t eStep, uint64_t dStep) :
    initStep(iStep), endStep(eStep), dumpStep(dStep)
{
    if (endStep <= initStep || dumpStep == 0)
        errorAll(
            Error::ARGUMENTERROR,
            "Invalid range: %d-%d:%d\n"
            "End of range must be at least as large as beginning,"
            " and the step size must be greater than 0.",
            initStep,
            endStep,
            dumpStep);

    if ((endStep - initStep) % dumpStep != 0)
        errorAll(
            Error::ARGUMENTERROR,
            "Invalid range: %d-%d:%d\n"
            "The difference between the last and first step must"
            " be divisible by the increment.",
            initStep,
            endStep,
            dumpStep);
    
    nSteps = (endStep - initStep) / dumpStep + 1UL;
}
StepRange::StepRange(const std::string& rangeString)
{
    const char message[] = "Invalid range syntax: %s\nMust be of form <init>-<end>:<dump>, e.g., 0-1000000:1000";

    if (rangeString.size() < 3)
        errorAll(Error::SYNTAXERROR, "Range string too small");

    if (rangeString[0] < '0' || rangeString[0] > '9')
        errorAll(Error::SYNTAXERROR, message, rangeString);

    std::ostringstream oss;
    bool init = false, end = false;
    for (const char c : rangeString)
    {
        if (c < '0' && c > '9' && c != '-' && c != ':')
            errorAll(Error::SYNTAXERROR, message, rangeString);
        
        if (c == '-')
        {
            std::string str = oss.str();
            if (str.size() == 0 || init || end)
                errorAll(Error::SYNTAXERROR, message, rangeString);

            initStep = std::stoull(str);
            init = true;

            std::ostringstream tmp;
            oss.swap(tmp);
        }
        else if (c == ':')
        {
            std::string str = oss.str();
            if (str.size() == 0 || !init || end)
                errorAll(Error::SYNTAXERROR, message, rangeString);

            endStep = std::stoull(str);
            end = true;

            std::ostringstream tmp;
            oss.swap(tmp);
        }
        else
        {
            oss << c;
        }
    }

    if (!init)
        errorAll(Error::SYNTAXERROR, message, rangeString);

    std::string str = oss.str();
    if (end)
    {
        dumpStep = std::stoull(str);
    }
    else
    {
        endStep = std::stoull(str);
        dumpStep = 1UL;
    }

    if (endStep <= initStep || dumpStep == 0)
        errorAll(
            Error::ARGUMENTERROR,
            "Invalid range: %d-%d:%d\n"
            "End of range must be at least as large as beginning,"
            " and the step size must be greater than 0.",
            initStep,
            endStep,
            dumpStep);

    if ((endStep - initStep) % dumpStep != 0)
        errorAll(
            Error::ARGUMENTERROR,
            "Invalid range: %d-%d:%d\n"
            "The difference between the last and first step must"
            " be divisible by the increment.",
            initStep,
            endStep,
            dumpStep);
    
    nSteps = (endStep - initStep) / dumpStep + 1UL;
}
}