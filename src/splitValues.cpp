#include "splitValues.hpp"

template<typename T>
std::pair<T, T> splitValues(T totalNumValues, int me, int nProcs) {
    T q = totalNumValues / nProcs;
    T r = totalNumValues % nProcs;
    T numValues = (me<r) ? (q+1) : q;
    T firstIndex = (me<r) ? me*(q+1) : r*(q+1) + (me-r)*q;
    return std::pair<T, T>(firstIndex, numValues);
}
