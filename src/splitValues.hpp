#pragma once

#include <cstdint>
#include <cstdlib>
#include <utility>

namespace MDPAT
{
    /*
     * Returns the first index and number of values for process me when split evenly
     * among nprocs processes.
     */
    std::pair<uint64_t, uint64_t> splitValues(uint64_t totalNumValues, int me, int nProcs)
    {
        std::pair<uint64_t, uint64_t> pair;
        const auto div = std::div(totalNumValues, nProcs);

        pair.first = me * (div.quot + 1);
        pair.second = div.quot + 1;
        
        if (me >= div.rem)
        {
            pair.first += div.rem - me;
            pair.second -= 1;
        }
        return pair;
    }
}