#pragma once

#include <array>
#include <iostream>
#include <vector>

template<int N>
int getLinearIndices(const std::array<int, N> & indices, const std::array<int, N> & dimLengths, const std::array<int, N> & newDims);

template<typename T, int N>
std::vector<T> permuteDims(const std::vector<T> & oldvec, const std::array<int, N> & dimLengths, const std::array<int, N> & newDims);
