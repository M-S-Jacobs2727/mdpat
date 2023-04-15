#include "permuteDims.hpp"

using std::size_t;
using std::vector;

namespace MDPAT
{
uint64_t getLinearIndex(const vector<uint64_t> & indices, 
                        const vector<uint64_t> & dimLengths, 
                        const vector<uint64_t> & newDims)
{
    auto N = indices.size();
    uint64_t index = 0;
    uint64_t tmp;
    for (int i = 0; i < N; ++i) {
        tmp = indices[newDims[i]];
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
    auto N = dimLengths.size();
    
    if (newDims.size() != N) {
        std::cerr << "ERROR: size of arguments `dimLengths` and `newDims` should be equal.\n";
        return -1
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

    std::array<int, N> indices = {};
    std::vector<T> tmpVec;
    tmpVec.resize(vec.size());

    size_t index;

    for (i = 0; i < prod; ++i) {
        index = getLinearIndex<N>(indices, dimLengths, newDims);
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


/*  Takes a vector on each proc, sorted in row-major order according to
 *  `dimLengths`, with each proc having a contiguous subset in the first
 *  dimension, and reorganizes the vectors according to the first dimension
 *  of `newDims`. 
 * 
 *  For example, if `dimLengths = {N/P, M, L}` where `P` is the number of
 *  processes `nprocs`, and `newDims = {1, 2, 0}, then the resulting vector
 *  will have size (M/P) * L * N, organized first by M, then L, then N.
 * 
 *  Assumes that the vectors should be in order according to the values
 *  of `me`.
 */
// TODO: return std::array of new dimLengths
// template<typename T, size_t N>
// int permuteDimsParallel(std::vector<T> & vec,
//                         std::array<int, N> & dimLengths,
//                         const std::array<int, N> & newDims,
//                         const int me,
//                         const int nprocs,
//                         MPI_Comm comm)
// {
//     // TODO: add checks for dimLengths to be the same on all procs (except first element)
//     // TODO: make sure that the first dim is different (i.e., newDims[0] != 0)

//     // Gatherv all vectors to proc 0

//     int mySize = vec.size();
//     int totalSize = mySize;
    
//     auto recvcounts = std::make_shared<int[]>(nprocs);
//     auto displs = std::make_shared<int[]>(nprocs);
//     recvcounts[me] = mySize;

//     if (me == 0)
//     {
//         MPI_Reduce(recvcounts, MPI_IN_PLACE, nprocs, MPI_INT, MPI_SUM, 0, comm);
//         for (int i = 1; i < nprocs; ++i)
//         {
//             displs[i] = displs[i-1] + recvcounts[i-1];
//             totalSize += recvcounts[i];
//         }
//     }
//     else
//     {
//         MPI_Reduce(recvcounts, recvcounts, nprocs, MPI_INT, MPI_SUM, 0, comm);
//     }

//     if (me == 0)
//     {
//         vec.resize(totalSize);
//         MPI_Gatherv(vec.data(), mySize, MPI_FLOAT, MPI_IN_PLACE, recvcounts, displs, MPI_FLOAT, 0, comm);
//     }
//     else
//     {
//         MPI_Gatherv(vec.data(), mySize, MPI_FLOAT, vec.data(), recvcounts, displs, MPI_FLOAT, 0, comm);
//     }

//     // permute m_Dims on proc 0
//     if (me == 0) {
//         dimLengths[0] = totalSize;
//         for (int i = 1; i < N; ++i)
//             dimLengths[0] /= dimLengths[i];
//         permuteDims(vec, dimLengths, newDims);
//     }

//     // Scatterv new vectors
//     size_t dimcount = dimLengths[newDims[0]];
//     auto [firstIndex, numValues] = splitValues(dimcount, me, nprocs);
//     auto sendcounts = std::make_shared<int[]>(nprocs);

//     if (me == 0)
//     {
//         displs[0] = 0;
//         sendcounts[0] = numValues;
//         for (size_t i = 1; i < nprocs; ++i) 
//         {
//             auto [firstIndexTmp, numValuesTmp] = splitValues(dimcount, i, nprocs);
//             sendcounts[i] = numValuesTmp;
//             displs[i] = displs[i-1] + sendcounts[i-1];
//         }
//         MPI_Scatterv(vec.data(), sendcounts, displs, MPI_FLOAT, MPI_IN_PLACE, recvcount, MPI_FLOAT, 0, comm);
//         vec.resize(numValues);
//     }
//     else
//     {
//         vec.resize(numValues);
//         MPI_Scatterv(vec.data(), sendcounts, displs, MPI_FLOAT, vec.data(), numValues, MPI_FLOAT, 0, comm);
//     }

//     return 0;
// }

// template<typename T, size_t N>
// int permuteDimsParallel(boost::multi_array<T, N> arr,
//                         const std::array<size_t, N> & newDims,
//                         size_t me,
//                         size_t nprocs,
//                         MPI_Comm comm)
// {
//     // Gatherv all vectors to proc 0
//     size_t i;

//     int * shape = arr.shape();
//     auto allShapes = std::make_unique<int[]>(nprocs*N);
//     for (i = 0; i < N; ++i)
//         allShapes[me*N + i] = shape[i];
    
//     // auto recvcounts = std::make_unique<size_t[]>(nprocs);
//     // auto displs = std::make_unique<size_t[]>(nprocs);
//     // recvcounts[me] = myNumElements;

//     if (me == 0)
//     {
//         MPI_Reduce(allShapes, MPI_IN_PLACE, nprocs*N, MPI_INT, MPI_SUM, 0, comm);
//         // MPI_Reduce(recvcounts, MPI_IN_PLACE, nprocs, MPI_INT, MPI_SUM, 0, comm);
//         for (int i = 1; i < nprocs; ++i)
//         {
//             if ()
//             displs[i] = displs[i-1] + recvcounts[i-1];
//             totalSize += recvcounts[i];
//         }
//     }
//     else
//     {
//         MPI_Reduce(recvcounts, recvcounts, nprocs, MPI_INT, MPI_SUM, 0, comm);
//     }

//     if (me == 0)
//     {
//         arr.resize(totalSize);
//         MPI_Gatherv(&arr[0], myNumElements, MPI_FLOAT, MPI_IN_PLACE, recvcounts, displs, MPI_FLOAT, 0, comm);
//     }
//     else
//     {
//         MPI_Gatherv(&arr[0], myNumElements, MPI_FLOAT, &arr[0], recvcounts, displs, MPI_FLOAT, 0, comm);
//     }

//     // permute m_Dims on proc 0
//     if (me == 0) {
//         dimLengths[0] = totalSize;
//         for (int i = 1; i < N; ++i)
//             dimLengths[0] /= dimLengths[i];
//         permuteDims(vec, dimLengths, newDims);
//     } 

//     // Scatterv new vectors
//     auto [firstIndex, numValues] = splitValues(dimLengths[newDims[0]], me, nprocs);

//     if (me == 0)
//     {
//         MPI_Scatterv();
//         vec.resize(numValues);
//     }
//     else
//     {
//         vec.resize(numValues);
//         MPI_Scatterv();
//     }
// }

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

}