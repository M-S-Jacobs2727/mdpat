#include "selectFromArray.hpp"

using std::vector;

namespace MDPAT
{
// The atomType should be in the typeCol column (0-indexed, default=1).
template<typename T, typename T2>
T2 selectAtomsByType(vector<T> & atoms, T2 nSteps, T2 nAtoms, int nCols, const vector<int> & selectedTypes, int typeCol) {
    vector<int> atomList;
    atomList.reserve(nAtoms);

    if (selectedTypes.empty()) return -1;

    auto tFirst = selectedTypes.begin();
    auto tLast = selectedTypes.end();
    for (int i = 0; i < nAtoms; ++i) {
        if (
            tLast != std::find(
                tFirst,
                tLast,
                (int)atoms[i*nCols + typeCol]
            )
        ) {
            atomList.push_back(i);
        }
    }

    vector<T> outVec;
    outVec.reserve(nSteps * atomList.size() * nCols);
    for (auto it = atoms.begin(); it < atoms.end(); it += nCols*nAtoms)
        for (auto ind : atomList)
            outVec.insert(outVec.end(), it + ind*nCols, it + (ind+1)*nCols);
    
    atoms.swap(outVec);
    vector<T>().swap(outVec);

    return (T2)atomList.size();
}

template<typename T, typename T2>
int selectColumns(vector<T> & atoms, T2 nSteps, T2 nAtoms, int nCols, vector<int> columns) {
    vector<T> outVec;
    outVec.reserve(nSteps*nAtoms*columns.size());

    if (columns.empty()) return -1;

    for (auto it = atoms.begin(); it < atoms.end(); it += nCols)
        for (auto c : columns)
            outVec.push_back(*(it+c));

    return 0;
}
}