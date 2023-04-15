#pragma once

#include <cstddef>
#include <cstdint>
#include <utility>

std::pair<uint64_t, uint64_t> splitValues(uint64_t totalNumValues, int me, int nProcs);
