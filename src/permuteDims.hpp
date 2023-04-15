#pragma once

#include <array>
#include <cstdint>
#include <iostream>
#include <vector>

#include <mpi.h>

#include "splitValues.hpp"

namespace MDPAT
{
uint64_t getLinearIndex(const std::vector<uint64_t> & indices, 
                        const std::vector<uint64_t> & dimLengths, 
                        const std::vector<uint64_t> & newDims);

template<typename T>
int permuteDims(std::vector<T> & vec,
                std::vector<uint64_t> const & dimLengths,
                std::vector<uint64_t> const & newDims
);

// template<typename T, std::size_t N>
// int permuteDimsParallel(boost::multi_array<T, N> arr,
//                         const std::array<std::size_t, N> & newDims,
//                         std::size_t me,
//                         std::size_t nprocs,
//                         MPI_Comm comm);
}
