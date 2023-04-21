#pragma once

#include <cstdint>
#include <iostream>
#include <memory>
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
                const std::vector<uint64_t> & dimLengths,
                const std::vector<uint64_t> & newDims
);

template<typename T>
int permuteDimsParallel(std::vector<T> & vec,
                        const std::vector<uint64_t> & oldDimLengths,
                        const std::vector<uint64_t> & newDims,
                        const int me,
                        const int nprocs,
                        const MPI_Comm comm);
}
