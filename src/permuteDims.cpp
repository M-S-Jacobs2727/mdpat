#include "permuteDims.hpp"

template<int N>
int getLinearIndices(
    const std::array<int, N> & indices, 
    const std::array<int, N> & dimLengths, 
    const std::array<int, N> & newDims
) {
    int index = 0;
    int tmp;
    for (int i = 0; i < nDims; ++i) {
        tmp = indices[newDims[i]];
        for (int j = i+1; j < nDims; ++j) {
            tmp *= dimLengths[newDims[j]];
        }
        index += tmp;
    }
    return index;
}

template<typename T, int N>
std::vector<T> permuteDims(
    std::vector<T> & vec, 
    const std::array<int, N> & dimLengths, 
    const std::array<int, N> & newDims
) {
    int i, j;
    int prod = 1;
    int isSame = 1;
    for (i = 0; i < N; ++i)
        prod *= dimLengths[i];
    if (prod != vec.size()) {
        std::cerr << "ERROR: size of vector should equal" <<
            " accumulated product of `dimLengths`." << std::endl;
        return -2;
    }

    for (i = 0; i < N; ++i) {
        if (newDims[i] != i) {
            isSame = 0;
            break;
        }
    }
    if (isSame)
        return -1;

    std::array<int, N> indices = {};
    std::vector<T> tmpVec;
    tmpVec.resize(vec.size());

    int index;

    for (i = 0; i < prod; ++i) {
        index = getLinearIndices(indices, dimLengths, newDims);
        tmpVec[index] = vec[i];
        ++indices[N-1];
        for (j = N-1; j > 0; --j) {
            if (indices[j] < dimLengths[j])
                break;
            indices[j] = 0;
            ++indices[j-1];
        }
    }
    
    tmpVec.swap(vec);
    vector<T>().swap(tmpVec);

    return 0;
}
