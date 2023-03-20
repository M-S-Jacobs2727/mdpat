#pragma once

#include <utility>

template<typename T>
std::pair<T, T> splitValues(T totalNumValues, int me, int nProcs);
