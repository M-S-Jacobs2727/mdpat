#include "splitValues.hpp"

namespace MDPAT
{

    /*
     * Returns the first index and number of values for process me when split evenly
     * among nprocs processes.
     */
    std::pair<uint64_t, uint64_t> splitValues(uint64_t totalNumValues, int me, int nProcs)
    {
        uint64_t q = totalNumValues / nProcs;
        uint64_t r = totalNumValues % nProcs;
        uint64_t numValues = (me < r) ? (q + 1) : q;
        uint64_t firstIndex = (me < r) ? me * (q + 1) : r * (q + 1) + (me - r) * q;
        return std::pair<uint64_t, uint64_t>(firstIndex, numValues);
    }
}