#include "permuteDims.hpp"

using std::vector;
using std::size_t;

namespace MDPAT
{
uint64_t getLinearIndex(const vector<uint64_t> & indices, 
                        const vector<uint64_t> & dimLengths, 
                        const vector<uint64_t> & newDims)
{
    auto N = indices.size();
    uint64_t index = 0;
    for (int i = 0; i < N; ++i) {
        uint64_t tmp = indices[newDims[i]];
        for (int j = i+1; j < N; ++j) {
            tmp *= dimLengths[newDims[j]];
        }
        index += tmp;
    }
    return index;
}

// TODO: return std::array of new dimLengths
template<typename T>
int permuteDims(vector<T> & vec,
                const vector<uint64_t> & dimLengths,
                const vector<uint64_t> & newDims)
{
    int i, j;
    int prod = 1;
    int isSame = 1;
    auto const N = dimLengths.size();
    
    if (newDims.size() != N) {
        std::cerr << "ERROR: size of arguments `dimLengths` and `newDims` should be equal.\n";
        return -1;
    }
    
    for (i = 0; i < N; ++i)
        prod *= dimLengths[i];
    if (prod != vec.size()) {
        std::cerr << "ERROR: size of vector should equal" <<
            " accumulated product of `dimLengths`.\n";
        return -2;
    }

    for (i = 0; i < N; ++i) {
        if (newDims[i] != i) {
            isSame = 0;
            break;
        }
    }
    if (isSame)
        return -3;

    std::vector<uint64_t> indices(N, 0);
    std::vector<T> tmpVec;
    tmpVec.resize(vec.size());

    uint64_t index;

    for (i = 0; i < prod; ++i) {
        index = getLinearIndex(indices, dimLengths, newDims);
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
    std::vector<T>().swap(tmpVec);

    return 0;
}

/*  
 * Takes a 1-D vector on each proc, representative of an N-D array sorted in
 * row-major order according to `dimLengths`, and permutes those dimensions
 * according to `newDims` (a permutation of {0, 1, ..., N-1}).
 *
 * The vector on each proc is assumed to be of shape D1_p x D2 x ... x DN.
 * D1_p is not assumed to be identical on each proc, but D2-DN are.
 *
 * For example, if `dimLengths = {N/P, M, L}` where `P` is the number of
 * processes `nprocs`, and `newDims = {1, 2, 0}, then the resulting vector
 * will have size (M/P) * L * N, organized first by M, then L, then N.
 *
 * Assumes that the vectors are order according to the values of `me`
 * (i.e., proc 0 should hold the first N/P values, then proc 1
 * holds the next N/P values, etc.)
 */
template<typename T>
int permuteDimsParallel(vector<T> & vec,
                        vector<uint64_t> & dimLengths,
                        const vector<uint64_t> & newDims,
                        const int me,
                        const int nprocs,
                        const MPI_Comm comm)
{
    // TODO: make more efficient. Currently, this naively gathers all data to proc 0, proc 0 
    //     executes `permuteDims`, then the data is redistributed. This can be improved with
    //     a clever algorithm, but that's a heavier task for later. This works.

    if (dimLengths.size() <= 1 || newDims.size() <= 1)
        return -1;
    
    if (dimLengths.size() != newDims.size())
        return -2;

    int i, j;
    const auto numDims = newDims.size();
    
    // Check that the newDims are different from the old ones
    bool isSame = true;
    for (i = 0; i < numDims; ++i)
    {
        if (newDims[i] >= numDims)
            return -3;
        for (j = 0; j < i; ++j)
            if (newDims[i] == newDims[j])
                return -4;
        if (isSame && newDims[i] != i)
            isSame = false;
    }
    if (isSame)
        return -5;

    // Check that the dimLengths are the same along all but the first dimension
    auto procDimLengths = std::make_unique<uint64_t[]>(nprocs * numDims);
    MPI_Allgather(dimLengths.data(), numDims, MPI_UNSIGNED_LONG, procDimLengths.get(), numDims, MPI_UNSIGNED_LONG, comm);
    isSame = true;
    for (i = 1; i < nprocs; ++i)
        for (j = 1; j < numDims; ++j)
            if (isSame && procDimLengths[i*numDims+j] != procDimLengths[j])
                isSame = false;
    if (!isSame)
        return -6;

    // If we're not permuting the split dimension, then just permute alone without communication
    if (newDims[0] == 0)
    {
        permuteDims(vec, dimLengths, newDims);
        return 0;
    }

    // Gatherv all vectors to proc 0
    auto mySize = vec.size();
    // int totalSize = mySize;

    vector<uint64_t> totalDimLengths(numDims, 0);
    for (i = 0; i < numDims; ++i)
        totalDimLengths[i] = dimLengths[i];

    MPI_Allreduce(MPI_IN_PLACE, totalDimLengths.data(), 1, MPI_UNSIGNED_LONG, MPI_SUM, comm);
    
    auto recvcounts = std::make_unique<int[]>(nprocs);
    auto displs = std::make_unique<int[]>(nprocs);
    recvcounts[0] = mySize;

    uint64_t prodOtherDims = dimLengths[1];
    for (i = 2; i < numDims; ++i)
        prodOtherDims *= dimLengths[i];
    
    MPI_Datatype dtype = mpi_get_type<T>();

    if (me == 0)
    {
        recvcounts[0] = mySize;
        displs[0] = 0;
        auto totalSize = mySize;
        for (i = 1; i < nprocs; ++i)
        {
            recvcounts[i] = procDimLengths[i*numDims] * prodOtherDims;
            displs[i] = totalSize;
            totalSize += recvcounts[i];
        }

        vec.resize(totalSize);
        MPI_Gatherv(MPI_IN_PLACE, mySize, dtype, vec.data(), recvcounts.get(), displs.get(), dtype, 0, comm);
    }
    else
    {
        MPI_Gatherv(vec.data(), mySize, dtype, vec.data(), recvcounts.get(), displs.get(), dtype, 0, comm);
    }

    // permute on proc 0
    if (me == 0)
        permuteDims(vec, totalDimLengths, newDims);

    MPI_Barrier(comm);

    for (i = 0; i < numDims; ++i)
        dimLengths[i] = totalDimLengths[newDims[i]];

    auto sendcounts = std::make_unique<int[]>(nprocs);
    for (i = 0; i < nprocs; ++i)
    {
        auto [firstIndex, numValues] = splitValues(dimLengths[0], i, nprocs);
        sendcounts[i] = numValues * prodOtherDims;
        displs[i] = firstIndex * prodOtherDims;
    }
    dimLengths[0] = sendcounts[me] / prodOtherDims;

    // Scatterv new vectors
    if (me == 0)
    {
        MPI_Scatterv(vec.data(), sendcounts.get(), displs.get(), dtype, MPI_IN_PLACE, sendcounts[me], dtype, 0, comm);
        vec.resize(sendcounts[me]);
    }
    else
    {
        vec.resize(sendcounts[me]);
        MPI_Scatterv(vec.data(), sendcounts.get(), displs.get(), dtype, vec.data(), sendcounts[me], dtype, 0, comm);
    }

    return 0;
}

template int permuteDims(vector<float> & vec,
                         const vector<uint64_t> & dimLengths,
                         const vector<uint64_t> & newDims);
template int permuteDims(vector<double> & vec,
                         const vector<uint64_t> & dimLengths,
                         const vector<uint64_t> & newDims);
template int permuteDims(vector<uint64_t> & vec,
                         const vector<uint64_t> & dimLengths,
                         const vector<uint64_t> & newDims);
template int permuteDims(vector<int> & vec,
                         const vector<uint64_t> & dimLengths,
                         const vector<uint64_t> & newDims);
template int permuteDimsParallel(vector<float> & vec,
                                 vector<uint64_t> & dimLengths,
                                 const vector<uint64_t> & newDims,
                                 const int me,
                                 const int nprocs,
                                 const MPI_Comm comm);
template int permuteDimsParallel(vector<double> & vec,
                                 vector<uint64_t> & dimLengths,
                                 const vector<uint64_t> & newDims,
                                 const int me,
                                 const int nprocs,
                                 const MPI_Comm comm);
template int permuteDimsParallel(vector<uint64_t> & vec,
                                 vector<uint64_t> & dimLengths,
                                 const vector<uint64_t> & newDims,
                                 const int me,
                                 const int nprocs,
                                 const MPI_Comm comm);
template int permuteDimsParallel(vector<int> & vec,
                                 vector<uint64_t> & dimLengths,
                                 const vector<uint64_t> & newDims,
                                 const int me,
                                 const int nprocs,
                                 const MPI_Comm comm);

}