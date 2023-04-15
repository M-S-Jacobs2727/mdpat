#pragma once

#include <array>
#include <iostream>
#include <vector>

// #include <boost/multi_array.hpp>
#include <mpi.h>

#include "splitValues.hpp"

namespace MDPAT
{
template<int N>
size_t getLinearIndex(const std::array<int, N> & indices, 
                      const std::array<int, N> & dimLengths, 
                      const std::array<int, N> & newDims);

template<typename T, int N>
int permuteDims(std::vector<T> & vec,
                std::array<int, N> const & dimLengths,
                std::array<int, N> const & newDims
);

// template <>
// int permuteDims<int, 2>(std::vector<int> & vec,
//                         std::array<int, 2> const & dimLengths,
//                         std::array<int, 2> const & newDims);

// template<typename T, std::size_t N>
// int permuteDimsParallel(boost::multi_array<T, N> arr,
//                         const std::array<std::size_t, N> & newDims,
//                         std::size_t me,
//                         std::size_t nprocs,
//                         MPI_Comm comm);
}
