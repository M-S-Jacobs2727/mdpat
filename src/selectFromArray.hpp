#pragma once

#include <algorithm>
#include <vector>

namespace MDPAT
{
    // The atomType should be in the typeCol column (0-indexed, default=1).
    template <typename T, typename T2>
    T2 selectAtomsByType(
        std::vector<T> &atoms,
        T2 nSteps,
        T2 nAtoms,
        int nCols,
        const std::vector<int> &selectedTypes,
        int typeCol = 1);

    template <typename T, typename T2>
    int selectColumns(
        std::vector<T> &atoms,
        T2 nSteps,
        T2 nAtoms,
        int nCols,
        std::vector<int> columns);
}